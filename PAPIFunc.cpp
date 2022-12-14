/* 
 * This file is part of the pebil project.
 * 
 * Copyright (c) 2010, University of California Regents
 * All rights reserved.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * PAPI Function Instrumentation
 * The instrumentation reads the env var FPAPI_HWC0, FPAPI_HWC1, ... up to 
 * FPAPI_HWC31, in that order,
 * but stops at the first that is not defined (no gaps allowed).
 * The env variables should be set to PAPI present events, e.g. PAPI_TOT_CYC.
 * Then, for each loop instrumented, the value of the counters specified is 
 * accumulated and reported on at the end of execution. The values are printed 
 * in a set_%d.meta_%d.ftpapiinst file.
 * "set" allows measuring compatible hw counters in batched without over-writing
 * existing data files.
 * There is no check that events are compatible, please use the 
 * papi_event_chooser to verify events compatibility.
 *
 * Usage example:
 *
 * pebil --tool FunctionTimer --app bench --lnc libpapifunc.so,/opt/papi/lib/
     libpapi.so --fbl functionblacklistfile
 *
 * export FPAPI_HWC0=PAPI_TOT_INS
 * export FPAPI_HWC1=PAPI_TOT_CYC
 * export FPAPI_HWC_SET_NUMBER=0  (0 if you want to call the set 0, 
    but can use any int)
 * ./bench.ftminst
 * Should see output files like: bench.set_0.meta_0.ftpapiinst
 *
 * The CPU clock frequency can either be hard-coded (see the CLOCK_RATE_HZ 
 * above) or it can be passed on as an env variable -- i.e., by defining the 
 * FPAPI_CPU_FREQ
 * env variable and setting it to the CPU frequency of the system in Hz. 
 * For CPU frequency scaling experiments, one need to set this variable to the 
 * nominal/base 
 * CPU frequency only; i.e., no need to set this to individual scaled CPU 
 * frequencies.
 */

#include <InstrumentationCommon.hpp>
#include <DataManager.hpp>
#include <DynamicInstrumentation.hpp>
#include <Metasim.hpp>
#include <ThreadedCommon.hpp>
#include <PAPIFunc.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <papi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>

using namespace std;

// mandel's CPU Frequency: you can either hard-code the frequency here
// or use the FPAPI_CPU_FREQ env var.
#define CLOCK_RATE_HZ 3200000000
static uint32_t timerCPUFreq = CLOCK_RATE_HZ;

DynamicInstrumentation* DynamicPoints = NULL;

static uint32_t hwcSetNumber = 0;

// Set with FPAPI_SHUTOFF - enables function shutoff
#define SHUTOFF_SET 0
static uint32_t shutoffCounters = SHUTOFF_SET;
// Set with FPAPI_ITERS - # of iterations allowed before averaging time 
// per visit to determine if instrumentation should be shutoff for a function
#define SHUTOFF_ITERS 100
static uint32_t shutoffIters = SHUTOFF_ITERS;
// Set with FPAPI_THRESHOLD - Average time per visit allowed for a function in
// microseconds
#define SHUTOFF_THRESHOLD 5000
static uint32_t timingThreshold = SHUTOFF_THRESHOLD;

inline uint64_t read_timestamp_counter() {
    unsigned low, high;
    __asm__ volatile ("rdtsc" : "=a" (low), "=d"(high));
    return ((unsigned long long)low | (((unsigned long long)high) << 32));
}

DataManager<FunctionPAPI*>* AllData = NULL;

FunctionPAPI* GenerateFunctionPAPI(FunctionPAPI* counters, uint32_t typ, 
  image_key_t iid, thread_key_t tid, image_key_t firstimage) {

    FunctionPAPI* retval;
    retval = new FunctionPAPI();
    retval->master = counters->master && typ == DataManagerType_Image;
    retval->application = counters->application;
    retval->extension = counters->extension;
    retval->functionCount = counters->functionCount;
    retval->functionNames = counters->functionNames;
    retval->functionHashes = counters->functionHashes;
    retval->tmpValues = new values_t[retval->functionCount];
    retval->accumValues = new values_t[retval->functionCount];
    retval->num = 0;
    retval->currentlyMeasuring = 0;
    retval->eventSet = PAPI_NULL;
    retval->functionTimerAccum = new uint64_t[retval->functionCount];
    retval->functionTimerLast = new uint64_t[retval->functionCount];
    retval->inFunctionP = new uint32_t[retval->functionCount];
    retval->functionEntryCounts = new uint64_t[retval->functionCount];
    retval->functionShutoff = new uint32_t[retval->functionCount];
    
    memset(retval->functionTimerAccum, 0, sizeof(uint64_t) * 
      retval->functionCount);
    memset(retval->functionTimerLast, 0, sizeof(uint64_t) *
      retval->functionCount);
    memset(retval->inFunctionP, 0, sizeof(uint32_t) * retval->functionCount);
    memset(retval->functionEntryCounts, 0, sizeof(uint64_t) * 
      retval->functionCount);
    memset(retval->functionShutoff, 0, sizeof(uint32_t) * 
      retval->functionCount);
    
    /* The CPU clock frequency can either be hard-coded (see the CLOCK_RATE_HZ 
     * above) or it can be passed on as an env variable -- i.e., by defining 
     * the FPAPI_CPU_FREQ env variable and setting it to the CPU 
     * frequency of the system in Hz. For CPU frequency scaling experiments, 
     * one need to set this variable to the nominal/base CPU frequency only; 
     * i.e., no need to set this to individual scaled CPU frequencies.
     */
    if (ReadEnvUint32("FPAPI_CPU_FREQ", &timerCPUFreq)) {
        inform << "Received custom TIMER_CPU_FREQ ***(in Hz)** from the user "
          ":: " << timerCPUFreq << endl;
    } else {
        timerCPUFreq = CLOCK_RATE_HZ;
    }
    
    /* Hardware counter measurements are usually done in multiple iterations; 
     * with each iteration measuring a set of compatible counters. To give a
     * name/id to different sets, use this env variable. By default, if the var
     * is not set, the set id is set to 0.
     */
    if (ReadEnvUint32("FPAPI_HWC_SET_NUMBER", &hwcSetNumber)) {
        inform << "Received hardware counter set number from user :: " << 
          hwcSetNumber << endl;
    }

    // Check if function shutoff is set
    if (!ReadEnvUint32("FPAPI_SHUTOFF", &shutoffCounters)) {
        shutoffCounters = SHUTOFF_SET;
    }

    // If function shutoff is set, check for shutoff configuration vars
    if (shutoffCounters) {
        if (!ReadEnvUint32("FPAPI_ITERS", &shutoffIters)) {
            shutoffIters = SHUTOFF_ITERS;
        }
        if (!ReadEnvUint32("FPAPI_THRESHOLD", &timingThreshold)) {
            timingThreshold = SHUTOFF_THRESHOLD;
        }
    }
    
    return retval;
}
#undef CLOCK_RATE_HZ
#undef SHUTOFF_SET
#undef SHUTOFF_ITERS
#undef SHUTOFF_THRESHOLD

void DeleteFunctionPAPI(FunctionPAPI* counters){
    delete counters->functionTimerAccum;
    delete counters->functionTimerLast;
    delete counters->inFunctionP;
    delete counters->functionEntryCounts;
    delete counters->functionShutoff;
}

uint64_t ReferenceFunctionPAPI(FunctionPAPI* counters){
    return (uint64_t)counters;
}

extern "C"
{

  // function entry instrumentation
  int32_t function_entry(uint32_t funcIndex, image_key_t* key) {
    thread_key_t tid = pthread_self();
  
      FunctionPAPI* counters = AllData->GetData(*key, pthread_self());
      assert(counters != NULL);
      assert(counters->functionTimerLast != NULL);
  
      int eventSet = counters->eventSet;
      
      // if its a recursive function call increment the recursion depth and 
      // entry count and return
      if (counters->inFunctionP[funcIndex] != 0) {
          counters->inFunctionP[funcIndex]++;
          counters->functionEntryCounts[funcIndex]++;
          return 0;
      }
  
      // See if there are any active functions
      // Read the counters and update all active functions
      if (counters->currentlyMeasuring != 0) {
          // No error checking to minimize overhead
          PAPI_read(eventSet, counters->tmpValues[funcIndex]);
      
          // Uncomment error checking for debugging purposes 
          //if (PAPI_read(eventSet, counters->tmpValues[funcIndex]) != PAPI_OK){
          //    fprintf(stderr, "Error reading the values!\n");
          //    exit(1);
          //}
  
          // associate the counter values to all the active functions
          for (std::set<int>::iterator it = counters->activeFunctions.begin();
            it != counters->activeFunctions.end(); ++it) {
              for (int i = 0; i < counters->num; i++) {
                  counters->accumValues[*it][i] += 
                    counters->tmpValues[funcIndex][i];
              }
          }
      }
    
      // increment the currentlyMeasuring count
      counters->currentlyMeasuring++;
  
      // insert this function to the actively measuring set
      counters->activeFunctions.insert(funcIndex);
  
      //initialize PAPI for each thread (if not already).
      if (!counters->num) { 
          int retval = PAPI_library_init(PAPI_VER_CURRENT);
          if (retval != PAPI_VER_CURRENT) {
              fprintf(stderr, "PAPI library init error!\n");
              exit(1);
          }
  
          // Create the Event Set 
          if (PAPI_create_eventset(&eventSet) != PAPI_OK) {
              fprintf(stderr, "PAPI create event set error!\n");
              exit(1);
          }
        
          int eventCode = counters->eventCode;
          // Check each environment variable
          while (counters->num < MAX_HWC) {
              char hwc_var[MAX_HWC];
              sprintf(hwc_var, "FPAPI_HWC%d", counters->num);
              char* hwc_name = getenv(hwc_var);
              // If env var is set, add it to the events
              if (hwc_name) {
                  retval = PAPI_event_name_to_code(hwc_name, counters->events + 
                    counters->num);
                  if (retval != PAPI_OK) {
                      fprintf(stderr, "Unable to determine code for hwc %d: %s,"
                        " %d\n", counters->num, hwc_name, retval);
                  } else {
                      fprintf(stderr, "Thread 0x%llx parsed counter %s\n", tid, 
                        hwc_name);
                  }
                  if (PAPI_add_event(eventSet, 
                    *(counters->events+counters->num)) != PAPI_OK) {
                      fprintf(stderr, "PAPI add event error! Either the "
                        "specified counter(s) not available or compatible.\n");
                      exit(1);
                  }
                  ++counters->num;
              // If env var is not set, end while loop
              } else {
                  // Set PAPI_TOT_CYC as a default
                  if (counters->num == 0) {
                      fprintf(stderr, "No counters defined in the env. Adding "
                        "PAPI_TOT_CYC as default. \n");
                      PAPI_add_event(eventSet, PAPI_TOT_CYC);
                      *(counters->events+counters->num) = PAPI_TOT_CYC;
                      ++counters->num;
                  } else {
                      break;
                  }
              }
              fprintf(stderr, "Parsed %d counters for thread: 0x%llx\n", 
                counters->num, tid);
          }
      }

      counters->inFunctionP[funcIndex]++;
      counters->functionEntryCounts[funcIndex]++;
      
  
      // if this is the first entry, start the measurements
      if (counters->papiMeasurementsStarted == 0) {
          // No error checking to minimize overhead
          PAPI_start(eventSet);

          // uncomment for debugging purposes
          //if (PAPI_start(eventSet) != PAPI_OK) {
          //    fprintf(stderr, "Error in PAPI start!\n");
          //    exit(1);
          //}

          // indicate that the measurements have started 
          counters->papiMeasurementsStarted = 1;
      } else {
          // else reset the counters (again doing it without the check)
            PAPI_reset(eventSet);
      }
      counters->eventSet = eventSet;
      counters->functionTimerLast[funcIndex] = read_timestamp_counter();
      return 0;
  }

  //function exit instrumentation
  int32_t function_exit(uint32_t funcIndex, image_key_t* key) {
    
      uint64_t now = read_timestamp_counter();
      thread_key_t tid = pthread_self();
      FunctionPAPI* counters = AllData->GetData(*key, pthread_self());
      int eventSet = counters->eventSet;
      PAPI_read(eventSet, counters->tmpValues[funcIndex]);
      uint32_t recDepth = counters->inFunctionP[funcIndex];
      
      if (recDepth == 0) {
          if (GetTaskId() == 0) {
              warn << "Thread " << AllData->GetThreadSequence(tid) << " Leaving"
                " never entered function " << funcIndex << ":" << 
                counters->functionNames[funcIndex] << ENDL;
              print_backtrace();
          }
          counters->inFunctionP[funcIndex] = 0;
          return 0; 
      } else if (recDepth < 0) {
          if (GetTaskId() == 0) {
              warn << "Negative call depth for " << 
                counters->functionNames[funcIndex] << ENDL;
          }
          counters->inFunctionP[funcIndex] = 0;
          return 0;
      }
      
      --recDepth;
      if (recDepth == 0) {
          uint64_t last = counters->functionTimerLast[funcIndex];
        
          counters->functionTimerAccum[funcIndex] += now - last;
          counters->functionTimerLast[funcIndex] = now;
          for (std::set<int>::iterator it=counters->activeFunctions.begin();
            it!=counters->activeFunctions.end(); ++it) {
  
              for (int i = 0; i < counters->num; i++) {
                  counters->accumValues[*it][i] += 
                    counters->tmpValues[funcIndex][i];
              }
          }
     
          counters->currentlyMeasuring--;
          // remove this function from the active function set.
          counters->activeFunctions.erase(funcIndex);
      
          // if there are active functions remaining, we need to reset the 
          // counters
          if (counters->currentlyMeasuring != 0) {
              PAPI_reset(eventSet);
          }
      }
    
      counters->inFunctionP[funcIndex] = recDepth;

      // Shutoff counters for this function?
      if (shutoffCounters) {
          if (counters->functionEntryCounts[funcIndex] % shutoffIters == 0) {
              double timePerVisit = ((double)counters->functionTimerAccum[
                funcIndex]) / ((double)counters->functionEntryCounts[funcIndex])
                / timerCPUFreq;

              if (timePerVisit < (((double)timingThreshold) / 1000000.0)) {
                  uint64_t imageSeq = AllData->GetImageSequence(*key);
                  AllData->WriteLock();
                  uint64_t this_key = GENERATE_UNIQUE_KEY(funcIndex, imageSeq,
                    PointType_functionExit);
                  uint64_t corresponding_entry_key = GENERATE_UNIQUE_KEY(
                    funcIndex, imageSeq, PointType_functionEntry);
  
                  set<uint64_t> inits;
                  inits.insert(this_key);
                  inits.insert(corresponding_entry_key);
                  DynamicPoints->SetDynamicPoints(inits, false);
                  counters->functionShutoff[funcIndex] = 1;
                  AllData->UnLock();
              }
          }
      }
    
      return 0;
  }
  
  static pthread_mutex_t dynamic_init_mutex = PTHREAD_MUTEX_INITIALIZER;  
  void* tool_dynamic_init(uint64_t* count, DynamicInst** dyn, bool* 
    isThreadedModeFlag) {
      pthread_mutex_lock(&dynamic_init_mutex);
      if (DynamicPoints == NULL) {
          DynamicPoints = new DynamicInstrumentation();
      }
      DynamicPoints->InitializeDynamicInstrumentation(count, dyn,
        isThreadedModeFlag);
      pthread_mutex_unlock(&dynamic_init_mutex);
      return NULL;
  }
  
  void* tool_mpi_init() {
      return NULL;
  }
  
  void* tool_thread_init(thread_key_t tid) {
      if (AllData) {
          if (DynamicPoints->IsThreadedMode()) {
              AllData->AddThread(tid);
          }
      } else {
          ErrorExit("Calling PEBIL thread initialization library for thread " 
            << hex << tid << " but no images have been initialized.", 
            MetasimError_NoThread);
      }
      return NULL;
  }
  
  void* tool_thread_fini(thread_key_t tid) {
      return NULL;
  }

  static pthread_mutex_t image_init_mutex = PTHREAD_MUTEX_INITIALIZER;  
  void* tool_image_init(void* args, image_key_t* key, ThreadData* td) {
  
      pthread_mutex_lock(&image_init_mutex);
      FunctionPAPI* counters = (FunctionPAPI*)args;
    
      if (AllData == NULL) {
          AllData = new DataManager<FunctionPAPI*>(GenerateFunctionPAPI, 
            DeleteFunctionPAPI, ReferenceFunctionPAPI);
      }
    
      AllData->AddImage(counters, td, *key);
    
      counters = AllData->GetData(*key, pthread_self());
    
      if (PAPI_num_hwctrs() < PAPI_OK) {
          fprintf(stderr, "PAPI initialization failed");
          return NULL;
      }

      set<uint64_t> inits;
      inits.insert(GENERATE_KEY(*key, PointType_inits));
      DynamicPoints->SetDynamicPoints(inits, false);
 
      pthread_mutex_unlock(&image_init_mutex);
    
      return NULL;
  }
  
  void* tool_image_fini(image_key_t* key) {
      image_key_t iid = *key;

      static bool finalized = false;
      if (finalized)
          return NULL;

      finalized = true;

      if (DynamicPoints != NULL) {
          delete DynamicPoints;
      }
  
      if (AllData == NULL) {
          ErrorExit("data manager does not exist. no images were intialized", 
            MetasimError_NoImage);
          return NULL;
      }
  
      FunctionPAPI* counters = AllData->GetData(iid, pthread_self());
      if (counters == NULL) {
          ErrorExit("Cannot retrieve image data using key " << dec << (*key), 
            MetasimError_NoImage);
          return NULL;
      }
  
      if (!counters->master) {
          printf("Image is not master, skipping\n");
          return NULL;
      }
  
      char outFileName[1024];
      sprintf(outFileName, "%s.set_%0d.meta_%0d.%s", counters->application, 
        hwcSetNumber, GetTaskId(), "ftpapiinst");
  
      FILE* outFile = fopen(outFileName, "w");
      if (!outFile) {
          cerr << "error: cannot open output file %s" << outFileName << ENDL;
          exit(-1);
      }
      
      //print title of PAPI event names
      char EventName[512];
      int i, j, retval;
      char** units;
      units = new char*[MAX_HWC];
      PAPI_event_info_t evinfo;
    
      for (i = 0; i < counters-> num; i++) {
          units[i] = new char[PAPI_MIN_STR_LEN];
          PAPI_event_code_to_name(*(counters->events + i), EventName);
          retval = PAPI_get_event_info(*(counters->events + i), &evinfo);
          if (retval != PAPI_OK) {
              cerr << "\n\t Error getting event info\n";
              exit(-1);
          }
          strncpy(units[i], evinfo.units, PAPI_MIN_STR_LEN);
          fprintf(outFile, "%s ", EventName);
      }
    
      fprintf(outFile, "\n\n");
  
      //print function name:
      //then next line is Thread: threadhex Time: time hwc values
      //to match lppapiinst style
      for (set<image_key_t>::iterator iit = AllData->allimages.begin(); iit != 
        AllData->allimages.end(); ++iit) {
          FunctionPAPI* imageData = AllData->GetData(*iit, pthread_self());
          //uint64_t imgHash = *iit;
          char** functionNames = imageData->functionNames;
          uint64_t functionCount = imageData->functionCount;
  
          for (uint64_t funcIndex = 0; funcIndex < functionCount; ++funcIndex) {
              char* fname = functionNames[funcIndex];
              fprintf(outFile, "%s:\n", fname);
              for (set<thread_key_t>::iterator tit = 
                AllData->allthreads.begin(); tit != AllData->allthreads.end(); 
                ++tit) {
                  FunctionPAPI* icounters = AllData->GetData(*iit, *tit);
                  
                  fprintf(outFile, "\tThread: %d\tTime: %f\tEntries: "
                    "%lld\tHash: 0x%llx\t", 
                    AllData->GetThreadSequence(*tit), (double)(icounters->functionTimerAccum[funcIndex]) / 
                    timerCPUFreq, icounters->functionEntryCounts[funcIndex], 
                    icounters->functionHashes[funcIndex]);
                  // Add * if function was shutoff
                  if (icounters->functionShutoff[funcIndex] == 1) {
                      fprintf(outFile, "*\t");
                  }
  
                  int hwc;
                  double scaledValue;
                  for (hwc = 0; hwc < icounters->num; ++hwc) {
                      PAPI_event_code_to_name(*(icounters->events + hwc), 
                        EventName);
                      retval = PAPI_get_event_info(*(icounters->events + hwc), 
                        &evinfo);
                      if (retval != PAPI_OK) {
                          cerr << "\n\t Error getting event info\n";
                          exit(-1);
                      }
                      if (strstr(units[hwc],"nJ")) {
                          scaledValue = (double)(icounters->accumValues[
                            funcIndex][hwc]/(1.0e9));
                          double watts = scaledValue / 
                            ((double)(icounters->functionTimerAccum[funcIndex])
                            / timerCPUFreq);
                          fprintf(outFile, "%.4f ", watts);
                      } else {
                          fprintf(outFile, "%lld ", 
                            icounters->accumValues[funcIndex][hwc]);
                      }
                  }
                  fprintf(outFile, "\n");
              }
          }
      }
      
      fflush(outFile);
      fclose(outFile);
      return NULL;
  }
};

// helpers borrowed from Simulation.cpp

bool ParsePositiveInt32(string token, uint32_t* value) {
    return ParseInt32(token, value, 1);
}

// returns true on success... allows things to continue on failure if desired
bool ParseInt32(string token, uint32_t* value, uint32_t min) {
    int32_t val;
    uint32_t mult = 1;
    bool ErrorFree = true;

    istringstream stream(token);
    if (stream >> val) {
        if (!stream.eof()) {
            char c;
            stream.get(c);

            c = ToLowerCase(c);
            if (c == 'k') {
	              mult = KILO;
            } else if (c == 'm') {
	              mult = MEGA;
            } else if (c == 'g') {
	              mult = GIGA;
            } else {
	              ErrorFree = false;
            }

            if (!stream.eof()) {
	              stream.get(c);

	              c = ToLowerCase(c);
	              if (c != 'b') {
	                  ErrorFree = false;
	              }
            }
        }
    }

    if (val < min){
      ErrorFree = false;
    }

    (*value) = (val * mult);
    return ErrorFree;
}

// returns true on success... allows things to continue on failure if desired
bool ParsePositiveInt32Hex(string token, uint32_t* value) {
    int32_t val;
    bool ErrorFree = true;

    istringstream stream(token);

    char c1, c2;
    stream.get(c1);
    if (!stream.eof()) {
        stream.get(c2);
    }

    if (c1 != '0' || c2 != 'x') {
        stream.putback(c1);
        stream.putback(c2);
    }

    stringstream ss;
    ss << hex << token;
    if (ss >> val)  {
    }

    if (val <= 0) {
        ErrorFree = false;
    }

    (*value) = val;
    return ErrorFree;
}


char ToLowerCase(char c) {
    if (c < 'a') {
        c += ('a' - 'A');
    }
    return c;
}

bool ReadEnvUint32(string name, uint32_t* var) {
    char* e = getenv(name.c_str());
    if (e == NULL) {
      return false;
      inform << "unable to find " << name << " in environment" << ENDL;
    }
    string s (e);
    if (!ParseInt32(s, var, 0)) {
      return false;
      inform << "unable to parse " << name << " in environment" << ENDL;
    }
    return true;
}

