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
 * Time spent in each function
 *
 * file per rank
 * function: total
 *   - per thread time
 * Timer
 */

#include <InstrumentationCommon.hpp>
#include <DataManager.hpp>
#include <DynamicInstrumentation.hpp>
#include <Metasim.hpp>
#include <ThreadedCommon.hpp>
#include <TimerFunctions.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>

using namespace std;

DataManager<FunctionTimers*>* AllData = NULL;

DynamicInstrumentation* DynamicPoints = NULL;

// by default, do not shut off function timing instrumentation.
// please set FTIMER_SHUTOFF to something other than zero to enable
//  function timer shutoff.
static uint32_t shutoffFunctionTimers=0;
// by default, the function timer tool allows 100 invocations of the
//  function and then averages the time per visit to determine 
//  whether to shut off timing measurement. please set FTIMER_ITERS
//  to control the number of iterations.
static uint32_t shutoffIters=100;
// by default, the function timer shuts of the timing for function if
//  the per visit time is less than 5000 microseconds.
// please set FTIMER_THRESHOLD env variable to control the number of
// microseconds per visit.
static uint32_t timingThreshold=5000;
static uint64_t timerCPUFreq=3200000000;
// HPE EPYC: note that if the env variable is not defined, we default to 
//    what is defined here:
#define CLOCK_RATE_HZ 3200000000
// Clark
//#define CLOCK_RATE_HZ 2 600 079 000

//#define CLOCK_RATE_HZ 2800000000
//#define CLOCK_RATE_HZ 3326000000
inline uint64_t read_timestamp_counter(){
    unsigned low, high;
    __asm__ volatile ("rdtsc" : "=a" (low), "=d"(high));
    return ((unsigned long long)low | (((unsigned long long)high) << 32));
}

static double diffTime(struct timeval t1, struct timeval t2)
{
    struct timeval diff;
    if(t2.tv_usec < t1.tv_usec) {
            diff.tv_usec = 1000000 + t2.tv_usec - t1.tv_usec;
            diff.tv_sec = t2.tv_sec - t1.tv_sec - 1;
    } else {
            diff.tv_usec = t2.tv_usec - t1.tv_usec;
            diff.tv_sec = t2.tv_sec - t1.tv_sec;
    }

    double time = (double)diff.tv_sec + (diff.tv_usec / 1000000.0);
    return time;
}

/*
 * When a new image is added, called once per existing thread
 * When a new thread is added, called once per loaded image
 *
 * timers: some pre-existing data
 * typ: ThreadTyp when called via AddThread
 *      ImageTyp when called via AddImage
 * iid: image the new data will be for
 * tid: thread the new data will be for
 * firstimage: key of first image created
 * 
 */
FunctionTimers* GenerateFunctionTimers(FunctionTimers* timers, uint32_t typ, image_key_t iid, thread_key_t tid, image_key_t firstimage) {

    FunctionTimers* retval;
    retval = new FunctionTimers();

    retval->master = timers->master && typ == DataManagerType_Image;
    //retval->master = timers->master && typ == AllData->ImageType;
    retval->application = timers->application;
    retval->extension = timers->extension;
    retval->functionCount = timers->functionCount;
    retval->functionNames = timers->functionNames;
    retval->functionHashes = timers->functionHashes;
    retval->functionTimerAccum = new uint64_t[retval->functionCount];
    retval->functionTimerLast = new uint64_t[retval->functionCount];
    retval->inFunction = new uint32_t[retval->functionCount];
    retval->functionEntryCounts = new uint64_t[retval->functionCount];
    retval->functionShutoff = new uint32_t[retval->functionCount]; 

    memset(retval->functionTimerAccum, 0, sizeof(*retval->functionTimerAccum) * retval->functionCount);
    memset(retval->functionTimerLast, 0, sizeof(*retval->functionTimerLast) * retval->functionCount);
    memset(retval->inFunction, 0, sizeof(*retval->inFunction) * retval->functionCount);
    memset(retval->functionEntryCounts, 0, sizeof(*retval->functionEntryCounts) * retval->functionCount);
    memset(retval->functionShutoff, 0, sizeof(*retval->functionShutoff) * retval->functionCount);

    retval->appTimeStart = timers->appTimeStart;
    retval->appTimeOfDayStart = timers->appTimeOfDayStart;
    retval->sanitize = timers->sanitize;

    // read in key environment variables
    if (!ReadEnvUint32("FTIMER_SHUTOFF", &shutoffFunctionTimers)){
        shutoffFunctionTimers=0;
    }


    // see if the FTIMER_CPU_FREQ env var is defined
    char * ftimeCPU = getenv("FTIMER_CPU_FREQ");
    if (ftimeCPU != NULL) {

        std::stringstream strStream;
        strStream << ftimeCPU;
        strStream >> timerCPUFreq;
        inform << "Got custom FTIMER_CPU_FREQ ***(in Hz)** from the user :: " 
          << timerCPUFreq << endl;
    } else {
        inform << "***Using the default CPU clock rate to calculate timings****"
          << CLOCK_RATE_HZ << endl;
        timerCPUFreq=CLOCK_RATE_HZ;
    }

    if(shutoffFunctionTimers) {
        if (!ReadEnvUint32("FTIMER_ITERS", &shutoffIters)){
            shutoffIters=100;
        }

        if (!ReadEnvUint32("FTIMER_THRESHOLD", &timingThreshold)){
            timingThreshold=5000;
        }

        //warn << "Dynamic Turning off of Function Timers is enabled." << ENDL;
        //warn << "Number of iterations to consider before averaging time per visit: " << shutoffIters << ENDL;
        //warn << "Timer per visit threshold is at: " << timingThreshold << " micro-seconds. " << ENDL;

    }

    return retval;
}

void DeleteFunctionTimers(FunctionTimers* timers){
    delete timers->functionTimerAccum;
    delete timers->functionTimerLast;
    delete timers->inFunction;
    delete timers->functionEntryCounts;
    delete timers->functionShutoff;
}

uint64_t ReferenceFunctionTimers(FunctionTimers* timers){
    return (uint64_t)timers;
}

extern "C"
{

    // start timer
    int32_t function_entry(uint32_t funcIndex, image_key_t* key) {
        thread_key_t tid = pthread_self();

        FunctionTimers* timers = AllData->GetData(*key, pthread_self());
        assert(timers != NULL);
        assert(timers->functionTimerLast != NULL);

        if(timers->inFunction[funcIndex] == 0){
            timers->functionEntryCounts[funcIndex]++;
            timers->functionTimerLast[funcIndex] = read_timestamp_counter();

            if(GetTaskId() == 0) {
                    //warn << "Thread " << AllData->GetThreadSequence(tid) << " Entering function " << funcIndex << ":" << timers->functionNames[funcIndex] << ENDL;
                    //print_backtrace();
            }

        } else if(GetTaskId() == 0) {
            //warn << "Thread " << AllData->GetThreadSequence(tid) << " Re-entering function " << timers->functionNames[funcIndex] << ENDL;
        }
        ++timers->inFunction[funcIndex];

        return 0;
    }

    // end timer
    int32_t function_exit(uint32_t funcIndex, image_key_t* key) {
        thread_key_t tid = pthread_self();
        uint64_t last, now;
        FunctionTimers* timers = AllData->GetData(*key, pthread_self());

        int32_t recDepth = timers->inFunction[funcIndex];
        if(recDepth == 0) {
            if(GetTaskId() == 0) {
                warn << "Thread " << AllData->GetThreadSequence(tid) << 
                  " Leaving never entered function " << funcIndex << ":" << 
                  timers->functionNames[funcIndex] << ENDL;
                print_backtrace();
            }
            timers->inFunction[funcIndex] = 0;
            return 0; 

        } else if(recDepth < 0) {
            if(GetTaskId() == 0) warn << "Negative call depth for " <<                        timers->functionNames[funcIndex] << ENDL;
            timers->inFunction[funcIndex] = 0;
            return 0;
        }

        --recDepth;
        if(recDepth == 0) {
            last = timers->functionTimerLast[funcIndex];
            now = read_timestamp_counter();
            timers->functionTimerAccum[funcIndex] += now - last;
            timers->functionTimerLast[funcIndex] = now;

            if(GetTaskId() == 0) {
                //warn << "Leaving function " << timers->functionNames[funcIndex] << ENDL;
            }
        }
        timers->inFunction[funcIndex] = recDepth;

        if(shutoffFunctionTimers) {
            if (timers->functionEntryCounts[funcIndex] % shutoffIters == 0){
                // time per visit is total t
                double timeInFunction = timers->functionTimerAccum[funcIndex];
                double numVisits = (double)timers->functionEntryCounts[
                  funcIndex];
                double timePerVisit= timeInFunction / numVisits / timerCPUFreq;

                if(timePerVisit < (((double)timingThreshold)/1000000.0)) {
                    uint64_t imageSeq = AllData->GetImageSequence(*key);
                    AllData->WriteLock();
                    uint64_t this_key = GENERATE_UNIQUE_KEY(funcIndex, imageSeq,
                      PointType_functionExit);
                    uint64_t corresponding_entry_key = GENERATE_UNIQUE_KEY(
                      funcIndex, imageSeq, PointType_functionEntry);

                    //warn << "Shutting off timing for function " << timers->functionNames[funcIndex] << "; time per visit averaged over " << timers->functionEntryCounts[funcIndex] << " entries is " << timePerVisit << "s; specified cut-off threshold is " << (((double)timingThreshold)/1000000.0) << "s." << ENDL;
                    set<uint64_t> inits;
                    inits.insert(this_key);
                    inits.insert(corresponding_entry_key);
                    DynamicPoints->SetDynamicPoints(inits, false); 
                    timers->functionShutoff[funcIndex] = 1;
                    AllData->UnLock();
                }
            }
        }
        return 0;
    }


    // Just after MPI_Init is called
    void* tool_mpi_init() {
        return NULL;
    }

    // Entry function for threads
    void* tool_thread_init(thread_key_t tid) {
        if (AllData){
            if(DynamicPoints->IsThreadedMode())
                AllData->AddThread(tid);
        } else {
            ErrorExit("Calling PEBIL thread initialization library for thread "
              << hex << tid << " but no images have been initialized.", 
              MetasimError_NoThread);
        }
        return NULL;
    }

    // Optionally? called on thread join/exit?
    void* tool_thread_fini(thread_key_t tid) {
        return NULL;
    }

    // Create mutex to ensure that Dynamics is initialized exactly once
    static pthread_mutex_t dynamic_init_mutex = PTHREAD_MUTEX_INITIALIZER;
    // initialize dynamic instrumentation
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

    // Create mutex to ensure that each image is initialized exactly once
    static pthread_mutex_t image_init_mutex = PTHREAD_MUTEX_INITIALIZER;
    // Called when new image is loaded
    void* tool_image_init(void* args, image_key_t* key, ThreadData* td) {

        pthread_mutex_lock(&image_init_mutex);
        FunctionTimers* timers = (FunctionTimers*)args;

        // If this is the first image, set up a data manager
        if (AllData == NULL){
            AllData = new DataManager<FunctionTimers*>(GenerateFunctionTimers, 
              DeleteFunctionTimers, ReferenceFunctionTimers);
        }

        // Check if added already
        if (AllData->allimages.count(*key) != 0) {
            pthread_mutex_unlock(&image_init_mutex);
            return NULL;
        }

        // image time
        timers->appTimeStart = read_timestamp_counter();
        gettimeofday(&timers->appTimeOfDayStart, NULL);

        // Add this image
        AllData->AddImage(timers, td, *key);

        // Remove this instrumentation
        // Must be done after the image is added, or threads may get to the 
        // instrumentation before the image is initialized
        set<uint64_t> inits;
        inits.insert(GENERATE_KEY(*key, PointType_inits));
        DynamicPoints->SetDynamicPoints(inits, false);

        pthread_mutex_unlock(&image_init_mutex);
        return NULL;
    }

    // 
    void* tool_image_fini(image_key_t* key) {

        image_key_t iid = *key;

        // Only print one file with data from all images
        static bool finalized = false;
        if (finalized)
            return NULL;

        finalized = true;

        if (DynamicPoints != NULL) {
            delete DynamicPoints;
        }

        if (AllData == NULL){
            ErrorExit("data manager does not exist. no images were intialized",
              MetasimError_NoImage);
            return NULL;
        }

        FunctionTimers* timers = AllData->GetData(iid, pthread_self());
        if (timers == NULL){
            ErrorExit("Cannot retrieve image data using key " << dec << (*key),
              MetasimError_NoImage);
            return NULL;
        }

        if (!timers->master){
            printf("Image is not master, skipping\n");
            return NULL;
        }


        uint64_t appTimeEnd = read_timestamp_counter();
        struct timeval tvEnd;
        gettimeofday(&tvEnd, NULL);

        char outFileName[1024];
        sprintf(outFileName, "%s.meta_%0d.%s", timers->application, GetTaskId(),
          timers->extension);

        FILE* outFile = fopen(outFileName, "w");
        if (!outFile){
            cerr << "error: cannot open output file %s" << outFileName << ENDL;
            exit(-1);
        }

        fprintf(outFile, "App timestamp time: %lld %lld %f\n", 
          timers->appTimeStart, appTimeEnd, (double)(appTimeEnd - timers->appTimeStart) / timerCPUFreq);
        fprintf(outFile, "App timeofday time: %lld %lld %f\n", 
          timers->appTimeOfDayStart.tv_sec, tvEnd.tv_sec, diffTime(timers->appTimeOfDayStart, tvEnd));
        // for each image
        //   for each function
        //     for each thread
        //       print time
        for (set<image_key_t>::iterator iit = AllData->allimages.begin(); iit != AllData->allimages.end(); ++iit) {
            FunctionTimers* imageData = AllData->GetData(*iit, pthread_self());

            char** functionNames = imageData->functionNames;
            uint64_t functionCount = imageData->functionCount;
            for (uint64_t funcIndex = 0; funcIndex < functionCount; ++funcIndex)            {
                if (timers->sanitize){
                    fprintf(outFile, "\n0x%llx:\t", timers->functionHashes[funcIndex]);
                    for (set<thread_key_t>::iterator tit = 
                      AllData->allthreads.begin(); tit != 
                      AllData->allthreads.end(); ++tit) {
                        FunctionTimers* timers = AllData->GetData(*iit, *tit);

                        if(timers->functionShutoff[funcIndex]==1) {
                            fprintf(outFile, "\tThread: %d\tTime: %f\tEntries: "
                              "%lld\tImage: %d\t*", AllData->GetThreadSequence(*tit), (double)(timers->
                              functionTimerAccum[funcIndex]) / timerCPUFreq, 
                              timers->functionEntryCounts[funcIndex],
                              AllData->GetImageSequence(*iit));
                        } else {
                            fprintf(outFile, "\tThread: %d\tTime: %f\tEntries: "
                              "%lld\tImage: %d\t", AllData->GetThreadSequence(*tit), (double)(timers->
                              functionTimerAccum[funcIndex]) / timerCPUFreq, 
                              timers->functionEntryCounts[funcIndex],
                              AllData->GetImageSequence(*iit));
                        }
                    }
                } else {
                    char* fname; 
                    fname = functionNames[funcIndex];//ELIZABETH REPLACE
                    fprintf(outFile, "\n%s:\t", fname);
                    for (set<thread_key_t>::iterator tit = 
                      AllData->allthreads.begin(); tit != 
                      AllData->allthreads.end(); ++tit) {
                        FunctionTimers* timers = AllData->GetData(*iit, *tit);

                        if(timers->functionShutoff[funcIndex]==1) {
                            fprintf(outFile, "\tThread: %d\tTime: %f\tEntries: "
                              "%lld\tHash: 0x%llx\tImage: %d*\t", AllData->GetThreadSequence(*tit), (double)(timers->
                              functionTimerAccum[funcIndex]) / timerCPUFreq, 
                              timers->functionEntryCounts[funcIndex], timers->
                              functionHashes[funcIndex],
                              AllData->GetImageSequence(*iit));
                        } else {
                            fprintf(outFile, "\tThread: %d\tTime: %f\tEntries: "
                              "%lld\tHash: 0x%llx\tImage: %d\t", AllData->GetThreadSequence(*tit), (double)(timers->
                              functionTimerAccum[funcIndex]) / timerCPUFreq, 
                              timers->functionEntryCounts[funcIndex], timers->
                              functionHashes[funcIndex],
                              AllData->GetImageSequence(*iit));
                        }
                    }
                }
            }
        }
        fflush(outFile);
        fclose(outFile);
        return NULL;
    }
};

// helpers borrowed from Simulation.cpp

bool ParsePositiveInt32(string token, uint32_t* value){
    return ParseInt32(token, value, 1);
}

// returns true on success... allows things to continue on failure if desired
bool ParseInt32(string token, uint32_t* value, uint32_t min){
    int32_t val;
    uint32_t mult = 1;
    bool ErrorFree = true;

    istringstream stream(token);
    if (stream >> val){
        if (!stream.eof()){
            char c;
            stream.get(c);

            c = ToLowerCase(c);
            if (c == 'k'){
                mult = KILO;
            } else if (c == 'm'){
                mult = MEGA;
            } else if (c == 'g'){
                mult = GIGA;
            } else {
                ErrorFree = false;
            }

            if (!stream.eof()){
                stream.get(c);

                c = ToLowerCase(c);
                if (c != 'b'){
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
bool ParsePositiveInt32Hex(string token, uint32_t* value){
    int32_t val;
    bool ErrorFree = true;

    istringstream stream(token);

    char c1, c2;
    stream.get(c1);
    if (!stream.eof()){
        stream.get(c2);
    }

    if (c1 != '0' || c2 != 'x'){
        stream.putback(c1);
        stream.putback(c2);        
    }

    stringstream ss;
    ss << hex << token;
    if (ss >> val){
    }

    if (val <= 0){
        ErrorFree = false;
    }

    (*value) = val;
    return ErrorFree;
}


char ToLowerCase(char c){
    if (c < 'a'){
        c += ('a' - 'A');
    }
    return c;
}

bool ReadEnvUint32(string name, uint32_t* var){
    char* e = getenv(name.c_str());
    if (e == NULL){
        return false;
        inform << "unable to find " << name << " in environment" << ENDL;
    }
    string s (e);
    if (!ParseInt32(s, var, 0)){
        return false;
        inform << "unable to parse " << name << " in environment" << ENDL;
    }
    return true;
}

