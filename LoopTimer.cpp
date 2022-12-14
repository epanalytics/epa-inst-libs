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
#include <LoopTimer.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <cstring>

using namespace std;

// HPE EPYC: note that if the env variable is not defined, we default to                                            
//    what is defined here:
static uint32_t timerCPUFreq=2200000000;

DynamicInstrumentation* DynamicPoints = NULL;

DataManager<LoopTimers*>* AllData = NULL;

// HPE Epyc
#define CLOCK_RATE_HZ 2200000000
// Clark
//#define CLOCK_RATE_HZ 2600079000

// Xeon Phi Max Rate
//#define CLOCK_RATE_HZ 1333332000
inline uint64_t read_timestamp_counter(){
    unsigned low, high;
    __asm__ volatile ("rdtsc" : "=a" (low), "=d"(high));
    return ((unsigned long long)low | (((unsigned long long)high) << 32));
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
LoopTimers* GenerateLoopTimers(LoopTimers* timers, uint32_t typ, image_key_t iid, thread_key_t tid, image_key_t firstimage) {

    LoopTimers* retval;
    retval = new LoopTimers();

    retval->master = timers->master && typ == DataManagerType_Image;
    retval->application = timers->application;
    retval->extension = timers->extension;
    retval->loopCount = timers->loopCount;
    retval->loopHashes = timers->loopHashes;
    retval->loopTimerAccum = new uint64_t[retval->loopCount];
    retval->loopTimerLast = new uint64_t[retval->loopCount];
    retval->entryCounts = new uint64_t[retval->loopCount];

    memset(retval->loopTimerAccum, 0, sizeof(uint64_t) * retval->loopCount);
    memset(retval->loopTimerLast, 0, sizeof(uint64_t) * retval->loopCount);
    memset(retval->entryCounts, 0, sizeof(uint64_t) * retval->loopCount);

    if (ReadEnvUint32("TIMER_CPU_FREQ", &timerCPUFreq)) {
        inform << "Got custom TIMER_CPU_FREQ ***(in MHz)** from the user :: " << timerCPUFreq << endl;
        // convert timerCPUFreq from MHz to Hz
        timerCPUFreq=timerCPUFreq*1000;
    } else {
        timerCPUFreq=CLOCK_RATE_HZ;
    }


    return retval;
}

void DeleteLoopTimers(LoopTimers* timers){
    delete timers->loopTimerAccum;
    delete timers->loopTimerLast;
    delete timers->entryCounts;
}

uint64_t ReferenceLoopTimers(LoopTimers* timers){
    return (uint64_t)timers;
}

extern "C"
{

    // start timer
    int32_t loop_entry(uint32_t loopIndex, image_key_t* key) {
        thread_key_t tid = pthread_self();

        LoopTimers* timers = AllData->GetData(*key, pthread_self());
        assert(timers != NULL);
        assert(timers->loopTimerLast != NULL);
        
        timers->loopTimerLast[loopIndex] = read_timestamp_counter();
        timers->entryCounts[loopIndex]++;
        return 0;
    }

    // end timer
    int32_t loop_exit(uint32_t loopIndex, image_key_t* key) {
        thread_key_t tid = pthread_self();

        LoopTimers* timers = AllData->GetData(*key, pthread_self());
        uint64_t last = timers->loopTimerLast[loopIndex];
        uint64_t now = read_timestamp_counter();
        timers->loopTimerAccum[loopIndex] += now - last;
        //timers->loopTimerLast = now;

        return 0;
    }

    // initialize dynamic instrumentation
    static pthread_mutex_t dynamic_init_mutex = PTHREAD_MUTEX_INITIALIZER;
    void* tool_dynamic_init(uint64_t* count, DynamicInst** dyn, bool* 
      isThreadedModeFlag) {
        pthread_mutex_lock(&dynamic_init_mutex);
        if (DynamicPoints == NULL)
            DynamicPoints = new DynamicInstrumentation();

        DynamicPoints->InitializeDynamicInstrumentation(count, dyn,
          isThreadedModeFlag);
        pthread_mutex_unlock(&dynamic_init_mutex);
        return NULL;
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
        ErrorExit("Calling PEBIL thread initialization library for thread " << hex << tid << " but no images have been initialized.", MetasimError_NoThread);
        }
        return NULL;
    }

    // Optionally? called on thread join/exit?
    void* tool_thread_fini(thread_key_t tid) {
        return NULL;
    }

    // Called when new image is loaded
    static pthread_mutex_t image_init_mutex = PTHREAD_MUTEX_INITIALIZER;
    void* tool_image_init(void* args, image_key_t* key, ThreadData* td) {
        pthread_mutex_lock(&image_init_mutex);
        LoopTimers* timers = (LoopTimers*)args;

        // If this is the first image, set up a data manager
        if (AllData == NULL){
            init_signal_handlers();
            AllData = new DataManager<LoopTimers*>(GenerateLoopTimers, DeleteLoopTimers, ReferenceLoopTimers);
        }

        // Add this image
        AllData->AddImage(timers, td, *key);

        // Remove this instrumentation
        set<uint64_t> inits;
        inits.insert(GENERATE_KEY(*key, PointType_inits));
        DynamicPoints->SetDynamicPoints(inits, false);

        pthread_mutex_unlock(&image_init_mutex);
        return NULL;
    }

    // 
    void* tool_image_fini(image_key_t* key) {
        image_key_t iid = *key;

        static bool finalized = false;
        if (finalized)
            return NULL;

        finalized = true;

        if (DynamicPoints != NULL) {
            delete DynamicPoints;
        }

        if (AllData == NULL){
            ErrorExit("data manager does not exist. no images were intialized", MetasimError_NoImage);
            return NULL;
        }

        LoopTimers* timers = AllData->GetData(iid, pthread_self());
        if (timers == NULL){
            ErrorExit("Cannot retrieve image data using key " << dec << (*key), MetasimError_NoImage);
            return NULL;
        }

        if (!timers->master){
            printf("Image is not master, skipping\n");
            return NULL;
        }

        char outFileName[1024];
        sprintf(outFileName, "%s.meta_%0d.%s", timers->application, GetTaskId(), timers->extension);
        FILE* outFile = fopen(outFileName, "w");
        if (!outFile){
            cerr << "error: cannot open output file %s" << outFileName << ENDL;
            exit(-1);
        }

        // for each image
        //   for each loop
        //     for each thread
        //       print time
        for (set<image_key_t>::iterator iit = AllData->allimages.begin(); iit != AllData->allimages.end(); ++iit) {
            LoopTimers* imageData = AllData->GetData(*iit, pthread_self());

            uint64_t imgHash = *iit;
            uint64_t* loopHashes = imageData->loopHashes;
            uint64_t loopCount = imageData->loopCount;

            for (uint64_t loopIndex = 0; loopIndex < loopCount; ++loopIndex){
                uint64_t loopHash = loopHashes[loopIndex];
                fprintf(outFile, "0x%llx:0x%llx:\n", imgHash, loopHash);
                for (set<thread_key_t>::iterator tit = AllData->allthreads.begin(); tit != AllData->allthreads.end(); ++tit) {
                    LoopTimers* timers = AllData->GetData(*iit, *tit);
                    fprintf(outFile, "\tThread: 0x%llx\tTime: %f\tEntries: %lld\n", *tit, (double)(timers->loopTimerAccum[loopIndex]) / timerCPUFreq, timers->entryCounts[loopIndex]);
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

