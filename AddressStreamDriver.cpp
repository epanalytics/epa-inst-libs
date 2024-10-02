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

#include <InstrumentationCommon.hpp>
#include <DataManager.hpp>
#include <DynamicInstrumentation.hpp>
#include <Metasim.hpp>
#include <ThreadedCommon.hpp>
#include <AddressStreamBase.hpp>
#include <AddressStreamDriver.hpp>

#include <AddressRange.hpp>
#include <CacheSimulation.hpp>
#include <ReuseDistanceASI.hpp>
#include <ScatterGatherLength.hpp>
#include <SpatialLocality.hpp>

#ifdef HAS_ARIEL_FRONTEND
#include <ArielFrontend.hpp>
#endif

#ifdef HAS_EPA_TOOLS
#include <DataCentricAddressRange.hpp>
#include <DataCentricSpatialLocality.hpp>
#include <DataCentricCacheSimulation.hpp>
#include <DataCentricReuseDistance.hpp>
#include <EntropyRange.hpp>
#include <DataCentricEntropyRange.hpp>
#include <PrefetchSimulation.hpp>
#include <SpatialLocalityPerMemOp.hpp>
#endif

#ifdef HAS_DATA_STRUCTURE_MODULE
#include <DataStructureModule.hpp>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string.h>
#include <assert.h>

using namespace std;

// Define directives to keep #ifdefs out of code
#ifdef HAS_ARIEL_FRONTEND
  #define GENERATE_ARIEL_TOOL new ArielFrontendTool()
#else
  #define GENERATE_ARIEL_TOOL 0
#endif

#ifdef HAS_EPA_TOOLS
  #define GENERATE_PREFETCH_TOOL new PrefetchSimulationTool()
  #define GENERATE_SPATIAL_MEMOP_TOOL new SpatialLocalityPerMemOpTool()
#else
  #define GENERATE_PREFETCH_TOOL 0
  #define GENERATE_SPATIAL_MEMOP_TOOL 0
#endif

#ifdef HAS_DATA_STRUCTURE_MODULE
  #define ENTER_TOOL(m) m->EnterTool()
  #define EXIT_TOOL(m, b) if(runDataCentric) m->ExitTool(b)
  #define GENERATE_DATA_TOOL(m) new m()
  #define GENERATE_MODULE(m) m = new DataStructureModule()
  #define GET_DATA_STRUCTURE_ID(m, a, l) m->GetDataStructureID(a, l)
  #define GET_NUM_DATA_STRUCTURES(m, parser) m->GetNumberOfDataStructures(parser)
  #define DELETE_MODULE(m) delete m
  #define PAUSE_MODULE(m) if(runDataCentric) m->PauseMemoryWrappers()
  #define PRINT_DATA_STRUCTURE_REPORT(m, s) if(runDataCentric) \
    m->PrintDataStructureReport(s)
  #define READLOCK(m, l) if(runDataCentric) m->ReadLock(l)
  #define REGISTER_TOOL(m) if (runDataCentric) m->RegisterThreadInDynamicTool()
  #define UNPAUSE_MODULE(m) if(runDataCentric) m->UnpauseMemoryWrappers()
  #define UNLOCK(m, l) if(runDataCentric) m->UnLock(l)
  #define WRITELOCK(m, l) if(runDataCentric) m->WriteLock(l)
#else
  #define ENTER_TOOL(m) false
  #define EXIT_TOOL(m, b) 0
  #define GENERATE_DATA_TOOL(m) 0
  #define GENERATE_MODULE(m) 0
  #define GET_DATA_STRUCTURE_ID(m, a, l) 0
  #define GET_NUM_DATA_STRUCTURES(m, parser) 0
  #define DELETE_MODULE(m) 0
  #define PAUSE_MODULE(m) 0
  #define PRINT_DATA_STRUCTURE_REPORT(m, s) 0
  #define READLOCK(m, l) 0
  #define REGISTER_TOOL(m) 0
  #define UNPAUSE_MODULE(m) 0
  #define UNLOCK(m, l) 0
  #define WRITELOCK(m, l) 0
#endif

// Default Constructor
AddressStreamDriver::AddressStreamDriver() {

    // Only run Cache Simulation by default
    runAddressRange = false;
    runArielFrontend = false;
    runCacheSimulation = true;
    runEntropyRange = false;
    runHardwarePrefetching = false;
    runReuseDistance = false;
    runScatterLength = false;
    runSpatialLocality = false;
    runSpatialLocalityPerMemOp = false;

    // Only run code-centric by default
    runCodeCentric = true;
    runDataCentric = false;

    // Create the vector to store the tools
    tools = new vector<AddressStreamTool*>();
    numCodeCentricTools = 0;

    numMemoryHandlers = 0;
    numCodeCentricMemoryHandlers = 0;

    // Create a parser for parsing
    parser = new StringParser();

    variableNameFile = "";

    GENERATE_MODULE(dataStructureModule);
}

AddressStreamDriver::~AddressStreamDriver() {
    if (sampler != NULL)
        delete sampler;
    if (liveMemoryAccessInstPointKeys != NULL)
        delete liveMemoryAccessInstPointKeys;
    if (dynamicPoints != NULL) {
        delete dynamicPoints;
    }
    if (parser != NULL)
        delete parser;
    for (vector<AddressStreamTool*>::iterator it = tools->begin(); it !=
      tools->end(); it++) {
          delete (*it);
    }
    tools->clear();
    delete tools; 
    delete fastData;

    DELETE_MODULE(dataStructureModule);
}

bool AddressStreamDriver::BuiltWithArielFrontend() {
    #ifdef HAS_ARIEL_FRONTEND
    return true;
    #else
    return false;
    #endif
}

bool AddressStreamDriver::BuiltWithDataStructureModule() {
    #ifdef HAS_DATA_STRUCTURE_MODULE
    return true;
    #else
    return false;
    #endif
}

bool AddressStreamDriver::BuiltWithEPATools() {
    #ifdef HAS_EPA_TOOLS
    return true;
    #else
    return false;
    #endif
}

void AddressStreamDriver::CreateFastData(uint64_t capacity) {
    assert(fastData == NULL);
    fastData = new FastData<AddressStreamStats*, BufferEntry*>(GetBufferIds,
      allData, capacity);
    assert(fastData);
}

void AddressStreamDriver::CreateSamplingMethod() {
    if (sampler != NULL) 
       delete sampler;
    uint32_t sampleMax;
    uint32_t sampleOn;
    uint32_t sampleOff;
    if (!(parser->ReadEnvUint32("METASIM_SAMPLE_MAX", &sampleMax))){
        sampleMax = DEFAULT_SAMPLE_MAX;
    }
    if (!(parser->ReadEnvUint32("METASIM_SAMPLE_OFF", &sampleOff))){
        sampleOff = DEFAULT_SAMPLE_OFF;
    }
    if (!(parser->ReadEnvUint32("METASIM_SAMPLE_ON", &sampleOn))){
        sampleOn = DEFAULT_SAMPLE_ON;
    }

    sampler = new SamplingMethod(sampleMax, sampleOn, sampleOff);
    sampler->Print();
}

void AddressStreamDriver::DeleteAllData() {
    delete allData;
}

// Return if the tool needed to be entered
// Pass to exit tool
bool AddressStreamDriver::EnterTool() {
    if (runDataCentric)
        return ENTER_TOOL(dataStructureModule);
    else
        return false;
}

// Pass value from enter tool
void AddressStreamDriver::ExitTool(bool needToExit) {
    EXIT_TOOL(dataStructureModule, needToExit);
}

bool AddressStreamDriver::HasLiveInstrumentationPoints(bool lock) {
    // if there are keys, then still live
    sampler->ReadLock(lock);
    bool stillLive = !(liveMemoryAccessInstPointKeys->empty());
    sampler->UnLock(lock);
    return stillLive;
}

// Should only be called once per image (only one thread should call it)
void* AddressStreamDriver::FinalizeImage(image_key_t* key) {
    image_key_t iid = *key;

    allData->SetTimer(iid, 1);
    SAVE_STREAM_FLAGS(cout);

#ifdef MPI_INIT_REQUIRED
    if (!IsMpiValid()){
        warn << "Process " << dec << getpid() << " did not execute "
          << "MPI_Init, will not print execution count files" << ENDL;
        RESTORE_STREAM_FLAGS(cout);
        return NULL;
    }
#endif

    if (allData == NULL){
        ErrorExit("data manager does not exist. no images were "
          "initialized", MetasimError_NoImage);
        return NULL;
    }

    AddressStreamStats* stats = (AddressStreamStats*)allData->GetData(iid,
      pthread_self());
    if (stats == NULL){
        ErrorExit("Cannot retreive image data using key " << dec << (*key),
          MetasimError_NoImage);
        return NULL;
    }

    // only print stats when the master image exits
    if (!stats->Master){
        RESTORE_STREAM_FLAGS(cout);
        return NULL;
    }

    // clear all threads' buffers
    for (set<thread_key_t>::iterator it = allData->allthreads.begin(); 
      it != allData->allthreads.end(); it++) {
        ProcessThreadBuffer(iid, (*it));
    }

    AddressStreamStats* statss = allData->GetData(iid, pthread_self());
    string fileName = "";
    fileName.append(statss->Application);
    PRINT_DATA_STRUCTURE_REPORT(dataStructureModule, fileName);
    
    // Create the reports 
    for (vector<AddressStreamTool*>::iterator it = tools->begin(); it !=
      tools->end(); it++) {
          AddressStreamTool* currentTool = (*it);
          currentTool->FinalizeTool(allData, sampler);
    }

    if(HasLiveInstrumentationPoints())
        ShutOffInstrumentationInAllBlocks();
    
    double t = (allData->GetTimer(*key, 1) - allData->GetTimer(*key, 0));
    inform << "CXXX Total Execution time for instrumented application " 
      << t << ENDL;
    // TODO Is this right?
    double m = (double)(GetNumMemoryHandlers() * sampler->GetAccessCount());
    inform << "CXXX - Address Stream Library - Memops simulated per "
      << "second: " << (m/t) << ENDL;
    RESTORE_STREAM_FLAGS(cout);
    return NULL;
}

void AddressStreamDriver::FinalizeThread(thread_key_t tid) {
    allData->FinishThread(tid);
}

// Look for a file that has variable names and locations
// This will get passed onto the data structure module
void AddressStreamDriver::GetAndSetVariableNameFile() {
    char* fileName = parser->GetEnv("METASIM_VAR_NAME_FILE");

    if (fileName != NULL) {
        variableNameFile = string(fileName);
    }
}

// Should only be called once per driver
void AddressStreamDriver::InitializeAddressStreamDriver(
  DataManager<AddressStreamStats*>* d) {

    // Initialize AllData
    allData = d;

    // Initialize Sampler
    CreateSamplingMethod();

    // Set up the tools!
    SetUpTools();

    // Set up the data structure module -- Must be done after SetUpTools
    // Otherwise runDataCentric will not be set
    if (runDataCentric) {
        GetAndSetVariableNameFile();
        SetUpDataStructureModule();
    }

}

// Initialize the Instrumentation Points that the sampler needs to turn off
// Requires sampler and allData!
void AddressStreamDriver::InitializeKeys() {
    if (liveMemoryAccessInstPointKeys == NULL)
        liveMemoryAccessInstPointKeys = new set<uint64_t>();
    assert(liveMemoryAccessInstPointKeys != NULL);

    // Get all the instrumentation points that put memory addresses in the 
    // buffer (PointType_bufferfill) so the sampler can turn them on/off
    set<uint64_t> keys;
    dynamicPoints->GetAllDynamicKeys(keys);
    sampler->WriteLock();
    for (set<uint64_t>::iterator it = keys.begin(); it != keys.end(); it++) {
        uint64_t k = (*it);
        if (GET_TYPE(k) == PointType_bufferfill && 
          allData->allimages.count(k) == 0){
            liveMemoryAccessInstPointKeys->insert(k);
        }
    }
    sampler->UnLock();

  
    // Disable them if sampling is turned off
    if (sampler->GetSamplingFrequency() == 0){
        inform << "Disabling all simulation-related instrumentation"
          " because METASIM_SAMPLE_ON is set to 0" << ENDL;
        ShutOffInstrumentationInAllBlocks();
    }

    // If EPA_SLICER_START_OFF is set then turn instrumentation off
    uint32_t startOff = 0;
    (void) parser->ReadEnvUint32("EPA_SLICER_START_OFF", &startOff);
    if (startOff != 0) {
        sampler->WriteLock();
        //SetDynamicPoints(false);
        dynamicPoints->SetDynamicPoints(*liveMemoryAccessInstPointKeys, false);
        sampler->UnLock();
    }

}

// Meant to only be called once per image (thus only one thread should 
// ever call this)
void* AddressStreamDriver::InitializeNewImage(image_key_t* iid, 
  AddressStreamStats* stats, ThreadData* threadData){

    // If already added, just return
    if (allData->allimages.count(*iid) != 0) {
        return NULL;
    }

    // Add image to allData
    allData->AddImage(stats, threadData, *iid);

    // If fastData not created yet, create it
    if(fastData == NULL) {
        CreateFastData(BUFFER_CAPACITY(stats));
    }

    // Add image to fastData
    fastData->AddImage();
    
    // Set image and thread id in stats
    stats->imageid = *iid;
    stats->threadid = allData->GenerateThreadKey();

    // Initialize instrumentation point keys for the new image
    InitializeKeys();

    // Start the application timer for this image
    allData->SetTimer(*iid, 0);

    // Remove initialization instrumentation points for this image
    dynamicPoints->SetDynamicPoint(GENERATE_KEY(*iid, PointType_inits), false);

    // TODO: Do we want to do this for every tool or just Ariel?
    fprintf(stderr, "ACC: Do openmp in InitializeNewImage\n");
    volatile int x = 0;
    #pragma omp parallel
    {
        #pragma omp critical
        {
            x += 1;
        }
    }

    return NULL;
}

void* AddressStreamDriver::InitializeNewThread(thread_key_t tid){
    RegisterThreadInDynamicTool();
    bool entered = EnterTool();
    SAVE_STREAM_FLAGS(cout);
    if (allData){
        if(dynamicPoints->IsThreadedMode()) {
            assert(fastData);
            fastData->Lock();
            allData->WriteLock();
            allData->AddThread(tid, false);
            fastData->AddThread(tid, false);
            allData->UnLock();
            fastData->UnLock();
        }
        InitializeSuspendHandler();
    } else {
        ErrorExit("Calling PEBIL thread initialization library for thread " 
          << hex << tid << " but no images have been initialized.", 
          MetasimError_NoThread);
    }

    RESTORE_STREAM_FLAGS(cout);
    ExitTool(entered);
    return NULL;
}

// Not thread-safe: a write lock must be held before using
void AddressStreamDriver::InitializeStatsWithNewHandlers(AddressStreamStats* 
  stats) {
    assert(GetNumMemoryHandlers() > 0);
    stats->Handlers = new MemoryStreamHandler*[GetNumMemoryHandlers()];
    bzero(stats->Handlers, sizeof(MemoryStreamHandler*) * 
      GetNumMemoryHandlers());

    for (vector<AddressStreamTool*>::iterator it = tools->begin(); it !=
      tools->end(); it++) {
          AddressStreamTool* currentTool = (*it);
          currentTool->AddNewHandlers(stats);
    }
}

// Not thread-safe: a write lock must be held before using
void AddressStreamDriver::InitializeStatsWithNewStreamStats(AddressStreamStats*
  stats) {
    assert(GetNumMemoryHandlers() > 0);
    // Create a StreamStats object for each test/memory handler
    stats->Stats = new StreamStats*[GetNumMemoryHandlers()];
    bzero(stats->Stats, sizeof(StreamStats*) * GetNumMemoryHandlers());

    uint32_t originalAllocCount = stats->AllocCount;

    uint32_t toolIndex = 0;
    for (vector<AddressStreamTool*>::iterator it = tools->begin(); it !=
      tools->end(); it++) {
          // For Data-Centric tools, set AllocCount to number of data
          // structures if no METASIM_DS_SIZE is set, else use the set
          // METASIM_DS_SIZE
        if (toolIndex == numCodeCentricTools) {
            stats->AllocCount = GET_NUM_DATA_STRUCTURES(dataStructureModule, parser);
        }
        toolIndex++;
        AddressStreamTool* currentTool = (*it);
        currentTool->AddNewStreamStats(stats);
    }

    // Reset AllocCount
    stats->AllocCount = originalAllocCount;
}

void AddressStreamDriver::PauseApplicationWrappers() {
    PAUSE_MODULE(dataStructureModule);
}

// Process all the addresses in each threads buffers. This goes against the 
// assumption that thread A can only access its own buffer and therefore it 
// can update it without locks. 
// We "get around this" by grabbing a read lock in "ProcessThreadBuffer"
// and grabbing the corresponding WriteLock here. Additionally, we stop 
// all threads so that they cannot add to their buffers. We grab all other 
// locks that are required for processing the buffer so that another thread 
// is not suspended while holding a required lock
void AddressStreamDriver::ProcessAllBuffers(ProcessBuffersExtra extra) {

    //Suspend all threads
    bool entered = EnterTool();
    // Get Data Structure Module FIRST, since process_buffer will 
    // take it first too
    WriteLockDSM();
    fastData->Lock();
    allData->WriteLock();
    sampler->WriteLock();
    SuspendAllThreads(allData->livethreads.size(),
      allData->livethreads.begin(), allData->livethreads.end());

    // Go through each image and thread and process their buffers
    for (set<image_key_t>::iterator iit = allData->allimages.begin();
      iit != allData->allimages.end(); iit++) {
        for (set<thread_key_t>::iterator it = allData->allthreads.begin(); 
          it != allData->allthreads.end(); it++) {
            ProcessThreadBuffer((*iit), (*it), false);
        }
    }

    // Do we need to turn instrumentation on/off after processing?
    if (extra == ProcessBuffersExtra_setDynamicOn) {
        dynamicPoints->SetDynamicPoints(*liveMemoryAccessInstPointKeys, true);
    } else if (extra == ProcessBuffersExtra_setDynamicOff) {
        dynamicPoints->SetDynamicPoints(*liveMemoryAccessInstPointKeys, false);
    }

    // resume all threads
    ResumeAllThreads();
    sampler->UnLock();
    allData->UnLock();
    fastData->UnLock();
    UnLockDSM();
    ExitTool(entered);
}

//void AddressStreamDriver::ProcessAllBuffers() {
//    for (set<image_key_t>::iterator iit = allData->allimages.begin();
//      iit != allData->allimages.end(); iit++) {
//        for (set<thread_key_t>::iterator it = allData->allthreads.begin();
//          it != allData->allthreads.end(); it++) {
//            ProcessThreadBuffer((*iit), (*it));
//        }
//    }
//}

// Thread-safe function
// Returns number of elements skipped
uint64_t AddressStreamDriver::ProcessBufferForEachHandler(image_key_t iid, 
  thread_key_t tid, uint32_t numElementsInBuffer, bool lock) {

    uint64_t numSkipped = 0;
    AddressStreamStats** faststats = fastData->GetBufferStats(tid, lock);
    assert(faststats != NULL);
    uint32_t elementIndex = 0; 
    for (elementIndex = 0; elementIndex < numElementsInBuffer; 
      elementIndex++){

        debug(assert(elementIndex >= 0));
        debug(assert(elementIndex < numElementsInBuffer));
        debug(assert(faststats[elementIndex]));
        debug(assert(faststats[elementIndex]->Stats));

        AddressStreamStats* stats = faststats[elementIndex];
        // If stats is null, then this buffer entry was empty.
        // This is possible when you have multiple threads, and sampling was 
        // turned on/off in the middle
        if (stats == NULL) {
            numSkipped++;
            continue;
        }
        assert(stats != NULL);
        uint64_t maxNumAddresses = stats->maxNumAddresses;

        BufferEntry* reference = BUFFER_ENTRY(stats, elementIndex);
        if (reference->imageid == 0){
            debug(assert(AllData->CountThreads() > 1));
            continue;
        }
        uint64_t memSeq = reference->memseq;
        uint64_t dataCentricSeq = reference->memseq;
        bool ldstFlag = reference->loadstoreflag;
        bool memvecFlag = false; 
        // for single memory entry, length is one
        uint64_t length = 1;
        if (reference->type == MEM_ENTRY) {
            if (reference->address != 0) { 
                stats->addressesForProcessing[0] = reference->address;
                if (runDataCentric)
                    dataCentricSeq = GET_DATA_STRUCTURE_ID(dataStructureModule, 
                      reference->address, false);
            } else {
                inform << "found address 0, skipping\n";
            }
        // end of if memory entry 
        } else if (reference->type == VECTOR_ENTRY ) {
            uint64_t currAddr;
            uint16_t mask = (reference->vectorAddress).mask;
            // for vec entry, length is determined by the mask.
            length = 0;
            memvecFlag = true;
            uint32_t loopCheck = (reference->vectorAddress).numIndices;
            // if this is false, we won't have space to store all of the
            // addresses in the addresses array.
            assert(loopCheck <= maxNumAddresses);
            for (int i = 0; i < loopCheck; i++) {
                if (mask % 2 == 1) {
                    currAddr = (reference->vectorAddress).base
                      + (reference->vectorAddress).indexVector[i]
                      * (reference->vectorAddress).scale;
                    //we start at 0 for length and increment when there
                    //is an address we are accessing so we can use that
                    //to keep track of where we are in the array as well
                    //as its final length
                    stats->addressesForProcessing[length] = currAddr;
                    length++;
                }// mask check 
                mask = (mask >> 1);
            }// for num of indices

            if (runDataCentric) {
                dataCentricSeq = GET_DATA_STRUCTURE_ID(dataStructureModule, 
                  stats->addressesForProcessing[0], false);
                // Check if we have addresses from different data structures --
                // If so, we're gonna need to refactor
                for (int i = 1; i < length; i++) {
                    if (dataCentricSeq != GET_DATA_STRUCTURE_ID(
                      dataStructureModule, stats->addressesForProcessing[i],
                      false)) {

                        fprintf(stderr, "WARNING: Multiple data structures in "
                          "a vector...data will be a little off. The fix will "
                          "require a small refactor.\n");
                    }
                }
            }
        // If EPAX_VECTOR_ENTRY: a masked contiguous vector memop
        // Currently either takes the form of
        // ld1d z1.d, p0/z, [x0, x1, LSL #3] or  (scalar plus scalar)
        // ld1d z0.d, p0/z, [x0, #1, mul vl]     (scalar plus immediate)
        // Were x0 is the base address and x1 is an index (optionally shifted
        // so index can go by 1 instead of by datatype size).
        // Immediates are offsets that are multiplied by the vector length.
        //
        // We are given the base address, the maximum amount accessed (i.e.,
        // before predication), and the number of addresses accessed. From this,
        // we can calculate which addresses are accessed. Then, we use the
        // given predicate register value to calculate which addresses were
        // actually loaded/stored.
#ifdef EPAX_INST_TOOL
        } else if (reference->type == EPAX_VECTOR_ENTRY) {
            // The first address accessed (x0 in examples)
            uint64_t memAddress = reference->epaxVectorAddress.memAddress;
            // Amount accessed in bytes
            uint64_t access = reference->epaxVectorAddress.sizeOfAccess / 8;
            // Number of elements accessed
            uint16_t numElems = reference->epaxVectorAddress.numElements;
            // Value of the predicate register (note: stored as bytes)
            uint8_t* predReg = reference->epaxVectorAddress.predReg;
            // Size of *accessed* datatype in bytes
            uint16_t memElemSize = access / numElems;
            // SVE vector length in bytes
            uint32_t vecLen = stats->SVEVectorLength / 8;
            // Size of Z register datatype in bytes (not the same as accessed
            // datatype! There could be a datatype conversion)
            uint16_t regElemSize = vecLen / numElems;

            // for epax_vector_entry, length is determined by the mask.
            length = 0;
            // For each memory address accessed, calculate the address.
            // Then, check the predicate to see if the address was
            // loaded/stored. If so, add it to the addressesForProcessing.
            for (int elemNum = 0; elemNum < numElems; elemNum++) {
                // Calculate the address accessed in memory
                uint64_t curAddress = memAddress + (elemNum * memElemSize);

                // Figure out which is the corresponding bit in the predicate
                // register. If the Z register datatype is 1 byte, then each
                // bit in the predicate register corresponds to an element in
                // the Z register. If the datatype is 2 bytes, then it is every
                // other bit. 4 bytes - every 4 bits. 8 bytes - every 8 bits.
                //
                // Pred reg is stored as 8 bits to an element. First figure out
                // which element of the pred reg we want, and then which bit
                // in that element is the corresponding one.
                uint16_t predRegElemToCheck = (elemNum * regElemSize) / 8;
                uint8_t predRegElem = predReg[predRegElemToCheck];
                uint8_t bitToCheck = (elemNum * regElemSize) % 8;

                // If corresponding bit is 1, then we access this address
                bool isOn = (predRegElem & (1 << bitToCheck)) != 0;
                if (isOn) {
                    stats->addressesForProcessing[length] = curAddress;
                    length++;
                }
            } // For each memory address accessed
            memvecFlag = false;

            if (runDataCentric) {
                dataCentricSeq = GET_DATA_STRUCTURE_ID(dataStructureModule, 
                  stats->addressesForProcessing[0], false);
                // Check if we have addresses from different data structures --
                // This should NOT happen since this a contiguous mem op
                for (int i = 1; i < length; i++) {
                    if (dataCentricSeq != GET_DATA_STRUCTURE_ID(
                      dataStructureModule, stats->addressesForProcessing[i],
                      false)) {

                        fprintf(stderr, "WARNING: Multiple data structures in "
                          "a vector...data will be a little off. This probably "
                          "should not be happening with this memory type.\n");
                    }
                }
            }
        // If EPAX_INDIRECT_ENTRY: a masked indirect vector memop (i.e.,
        // scatter and gathers)
        // Currently either takes the form of
        // ld1d z0.d, p0/z, [z1.d, #N]   (vector plus immediate)
        //   where z1 is a vector of base addresses and #N is a fixed offset
        // ld1d z0.d, p0/z, [x0, z1.d]   (scalar plus vector)
        //  where x0 is a base address and z1 is an index (optionally extended
        //  and/or shifted)
        //
        // For the vector plus immediate case, we are given the base vector
        // (z1 in the example) and the immediate to add to each value in the
        // base vector.
        // For the scalar plus scalar case, we are given the base address (x0
        // in the example), the index vector (z1 in the example), and whether
        // to sign extend/shift the index.
        // We are also given the number if addresses accessed (before
        // predication), and value of the predicate register so we can
        // determine which addresses were actually loaded/stored.
        } else if (reference->type == EPAX_INDIRECT_ENTRY) {
            // The base address, if it has one (otherwise it is 0)
            uint64_t baseAddress = reference->epaxIndirectAddress.baseAddress;
            // Does the index need to be extended?
            uint8_t doesExtension =
              reference->epaxIndirectAddress.doesExtension;
            // If the index is extended, is it signed or unsigned
            uint8_t signedExtension =
              reference->epaxIndirectAddress.signedExtend;
            // Amount to shift the index (after any extension) -- 0 is no shift
            uint8_t shiftAmount = reference->epaxIndirectAddress.shiftAmount;
            // Amount to add to base address (in vector plus immediate case)
            // 0 if no immediate
            uint8_t immediate = reference->epaxIndirectAddress.immediate;
            // The number of elements that will be accessed
            uint16_t numElements = reference->epaxIndirectAddress.numElements;
            // Value of the predicate register (note: stored as bytes)
            uint8_t* predReg = reference->epaxIndirectAddress.predReg;
            // Value of the Z register (either the base addresses or the
            // indices) Note: Stored as little endian!
            uint8_t* baseVector = reference->epaxIndirectAddress.baseVector;
            // SVE vector length in bits
            uint32_t vecLen = stats->SVEVectorLength;
            // Size of Z register datatype in bits (z0 or z1, they should be
            // the same datatype)
            uint32_t elemSize = vecLen / numElements;

            // Calculate the addresses that could have been accessed and store
            // them in a vector (addresses):
            //    * Convert the Z register (either baseAddresses or indices)
            //      from little endian (bytes) to big endian (64-bit ints)
            //      so we can use the value
            //    * Do a signed extension, if required
            //    * Do a shift, if required
            //    * Add it to the immediate, if required
            //    * Add it to the base address, if required
            std::vector<uint64_t> addresses;
            for (size_t i = 0 ; i < numElements; i++) {
                uint64_t valToPush;
                if (elemSize == 8) {
                    valToPush = (uint64_t) baseVector[i];
                    // if we do a sign extension, back fill with 1s if we have
                    // a 1 in the most significant bit
                    if (doesExtension == 1 && signedExtension == 1) {
                        uint8_t bitToExtend = valToPush & 0x80; // signBit
                        if (bitToExtend != 0) { // fill with 1s
                            valToPush | 0xffffffffffffff00;
                        }
                    }
                } else if (elemSize == 16) {
                    valToPush = baseVector[i*2];
                    valToPush |= (((uint64_t) baseVector[(i*2)+1]) << 8);
                    // If the address is supposed to be sign extended
                    if (doesExtension == 1 && signedExtension == 1) {
                        uint16_t bitToExtend = valToPush & 0x8000;
                        if (bitToExtend != 0) { // fill with 1s
                            valToPush | 0xffffffffffff0000;
                        }
                    }
                } else if (elemSize == 32) {
                    valToPush = baseVector[i*4];
                    valToPush |= (((uint64_t) baseVector[(i*4)+1]) << 8);
                    valToPush |= (((uint64_t) baseVector[(i*4)+2]) << 16);
                    valToPush |= (((uint64_t) baseVector[(i*4)+3]) << 24);
                    // If the address is supposed to be sign extended
                    if (doesExtension == 1 && signedExtension == 1) {
                        uint32_t bitToExtend = valToPush & 0x80000000;
                        if (bitToExtend != 0) { // fill with 1s
                            valToPush | 0xffffffff00000000;
                        }
                    }
                } else if (elemSize == 64) {
                    valToPush = baseVector[i*8];
                    valToPush |= (((uint64_t) baseVector[(i*8)+1]) << 8);
                    valToPush |= (((uint64_t) baseVector[(i*8)+2]) << 16);
                    valToPush |= (((uint64_t) baseVector[(i*8)+3]) << 24);
                    valToPush |= (((uint64_t) baseVector[(i*8)+4]) << 32);
                    valToPush |= (((uint64_t) baseVector[(i*8)+5]) << 40);
                    valToPush |= (((uint64_t) baseVector[(i*8)+6]) << 48);
                    valToPush |= (((uint64_t) baseVector[(i*8)+7]) << 56);
                    // can't sign extend 64 bits
                } else {
                    assert(false && "Don't know how to read the Z register" 
                      " with the given element size");
                }

                // Do shift: lsl shift amount (if 0, still correct)
                valToPush = valToPush << shiftAmount;
                // Add optional immediate (if 0, still correct)
                // Note: The immediate is offset * mbytes in ARM documentation
                valToPush = valToPush + immediate;
                // Add to the base address (if 0, still correct)
                valToPush = valToPush + baseAddress;
                // Push the final address!
                addresses.push_back(valToPush);
            }

            // TODO remove when nolonger necessary?
            assert(addresses.size() == numElements);

            // for epax_indirect_entry, length is determined by the mask.
            length = 0;
            // For each memory address accessed, check the predicate to see if
            // the address was loaded/stored. If so, add it to the
            // addressesForProcessing.
            for (int elemNum = 0; elemNum < numElements; elemNum++) {
                // Get the address accessed
                uint64_t curAddress = addresses[elemNum];

                // Figure out which is the corresponding bit in the predicate
                // register. If the Z register datatype is 1 byte, then each
                // bit in the predicate register corresponds to an element in
                // the Z register. If the datatype is 2 bytes, then it is every
                // other bit. 4 bytes - every 4 bits. 8 bytes - every 8 bits.
                //
                // Pred reg is stored as 8 bits to an element. First figure out
                // which element of the pred reg we want, and then which bit
                // in that element is the corresponding one.
                uint16_t predRegElemToCheck = (elemNum * elemSize) / 64;
                uint8_t predRegElem = predReg[predRegElemToCheck];
                uint8_t bitToCheck = (elemNum * (elemSize / 8)) % 8;

                // If corresponding bit is 1, then we access this address
                bool isOn = (predRegElem & (1 << bitToCheck)) != 0;
                if (isOn) {
                    stats->addressesForProcessing[length] = curAddress;
                    length++;
                }
            } // For each memory address accessed
            
            memvecFlag = true;

            if (runDataCentric) {
                dataCentricSeq = GET_DATA_STRUCTURE_ID(dataStructureModule, 
                  stats->addressesForProcessing[0], false);
                // Check if we have addresses from different data structures --
                // If so, we're gonna need to refactor
                for (int i = 1; i < length; i++) {
                    if (dataCentricSeq != GET_DATA_STRUCTURE_ID(
                      dataStructureModule, stats->addressesForProcessing[i],
                      false)) {

                        fprintf(stderr, "WARNING: Multiple data structures in "
                          "a vector...data will be a little off. The fix will "
                          "require a small refactor.\n");
                    }
                }
            }
#endif // EPAX_INST_TOOL
        } // end of epax indirect address

        debug(assert(length <= maxNumAddresses));

        // Process for each memory handler
        for (uint32_t handlerIndex = 0; handlerIndex < GetNumMemoryHandlers(); 
          handlerIndex++) {
            MemoryStreamHandler* handler = stats->Handlers[handlerIndex];
            StreamStats* ss = stats->Stats[handlerIndex];

            // If this is the first data-centric handler, then change the 
            // memop ID to the data structure ID
            if (handlerIndex == numCodeCentricMemoryHandlers) {
                memSeq = dataCentricSeq; 
            }

            if (handlerIndex >= numCodeCentricMemoryHandlers) {
                ss->SetIsCodeCentric(false);
            }

            if (reference->type == INSN_COUNT) {
                handler->ProcessInstructions((void*)ss, memSeq,
                  reference->address);
            } else {
                // maxNumAddresses is the allocated size of the array when it
                // was created, the length is the number of actual elements used
                (void) handler->Process((void*)ss, memSeq, ldstFlag,
                  stats->addressesForProcessing, length, memvecFlag);
            }
        }// for number of handlers

        // 0 out addresses array to prevent passing stale data
        memset(stats->addressesForProcessing, 0, sizeof(uint64_t) *
          maxNumAddresses);
    }// for elements in the buffer

    return numSkipped;
}

// Thread-safe
// ProcessThreadBuffer
// Input: image ID and thread ID
// Return: none
// Side effects:
//   * The buffer should always be reset (bring "current" back to the front)
//   * If all address collection is shut off, then don't do anything else
//   * If sampling is "on" (collecting addresses), process the buffer through
//     each handlers Process function. Then, check if we should shut off any
//     address collection
//   * If sampling is "off", tell the handlers how many addresses were not 
//     collected
//   * Switch sampling on/off depending on sampler settings
void* AddressStreamDriver::ProcessThreadBuffer(image_key_t iid, thread_key_t 
  tid, bool suspend) {

    // If we don't need to suspend, then we already have locks
    bool lock = suspend;

#define DONE_WITH_BUFFER(...) BUFFER_CURRENT(stats) = 0;  return NULL;

    // Prevent another thread from executing this code for this thread's 
    // buffer at the same time as this thread. This currently can only 
    // happen during a data-centric, coming from ProcessAllBuffers.
    // ProcessAllBuffers takes the data structure module read lock, so 
    // getting a write lock would prevent it from processing this buffer.
    // We grab the data structure module read lock so that other threads can 
    // process their own buffers concurrently and so that we don't change 
    // anything for a code-centric-only run
    ReadLockDSM(lock);

    // Check if we are sampling
    // Thread-safe: Sampling method protected with lock
    bool isSampling;
    isSampling = sampler->CurrentlySampling(lock);

    assert(iid);
    if (allData == NULL){
        ErrorExit("data manager does not exist. no images were initialized",
          MetasimError_NoImage);
        return NULL;
    }

    // Buffer is shared between all images
    debug(inform << "Getting data for image " << hex << iid << " thread " 
      << tid << ENDL);

    // Thread-safe call
    AddressStreamStats* stats = (AddressStreamStats*)allData->GetData(iid, 
      tid, lock);

    // Thread-safe: Each thread has its own stats
    if (stats == NULL){
        ErrorExit("Cannot retreive image data using key " << dec << iid, 
          MetasimError_NoImage);
        return NULL;
    }

    // Thread-safe: Each thread has its own stats
    uint64_t numElements = BUFFER_CURRENT(stats);
    uint64_t capacity = BUFFER_CAPACITY(stats);

    // Thread-safe call
    uint32_t threadSeq = allData->GetThreadSequence(tid, lock);

    debug(inform << "Thread " << hex << tid << TAB << "Image " << hex 
      << iid << TAB << "Counter " << dec << numElements << TAB 
      << "Capacity " << dec << capacity << TAB << "Total " << dec 
      << sampler->GetAccessCount() << ENDL);

    // If there is no more instrumentation, return
    // Thread-Safe call
    if (!HasLiveInstrumentationPoints(lock)){
        UnLockDSM(lock);
        DONE_WITH_BUFFER();
    }

    if (isSampling){
        // Refresh FastStats so it can be used
        // Thread-safe call
        BufferEntry* buffer = &(stats->Buffer[1]);
        fastData->Refresh(buffer, numElements, tid, lock);

        // Process the buffer for each memory handler
        // Thread-safe call
        uint64_t numSkipped = ProcessBufferForEachHandler(iid, tid, 
          numElements, lock);
        if (numSkipped > 0) {
            for (uint32_t i = 0; i < GetNumMemoryHandlers(); i++) {
                MemoryStreamHandler* m = stats->Handlers[i];
                m->SkipAddresses(numSkipped);
            }
        }

        // Shut off any instrumentation if sample max is hit
        // Thread-safe: Calls thread-safe functions
        ShutOffInstrumentationInMaxedGroups(iid, tid, suspend);

    // if not sampling            
    } else {
        // Let each handler know that addresses were skipped
        // Thread-safe calls since each thread has its own stats/handlers
        for (uint32_t i = 0; i < GetNumMemoryHandlers(); i++) {
            MemoryStreamHandler* m = stats->Handlers[i];
            m->SkipAddresses(numElements);
        }
    }

    // Turn sampling on/off
    // Sampler is thread-safe
    if (sampler->SwitchesMode(numElements, lock)){
        if (suspend) {
            allData->ReadLock();
            // We are modifiying dynamic points. Use the sampler write 
            // lock to protect this action
            sampler->WriteLock();
            SuspendAllThreads(allData->livethreads.size(),
              allData->livethreads.begin(), allData->livethreads.end());
        }
        dynamicPoints->SetDynamicPoints(*liveMemoryAccessInstPointKeys,
          !(isSampling));
        if (suspend) {
            ResumeAllThreads();
            sampler->UnLock();
            allData->UnLock();
        }
    }

    // Thread-safe
    sampler->IncrementAccessCount(numElements, lock);

    // Wipe the buffer before exitting to prevent use of stale addresses later
    // on. Start with element 1, since the 0 element has metadata
    memset(&(stats->Buffer[1]), 0, sizeof(BufferEntry) * capacity);

    UnLockDSM(lock);
    DONE_WITH_BUFFER();
}

void AddressStreamDriver::ReadLockDSM(bool lock) {
    READLOCK(dataStructureModule, lock);
}

void AddressStreamDriver::RegisterThreadInDynamicTool() {
    REGISTER_TOOL(dataStructureModule);
}

void AddressStreamDriver::SetUpDataStructureModule() {
#ifdef HAS_DATA_STRUCTURE_MODULE
    int32_t stackDepth;
    bool setDepth = parser->ReadEnvInt32("METASIM_UNWIND_DEPTH", &stackDepth);
    dataStructureModule->WriteLock();
    dataStructureModule->CreateContainer();
    dataStructureModule->SetDriver(this);
    dataStructureModule->SetVariableNameFile(variableNameFile);
    dataStructureModule->ParseVariableFile();
    dataStructureModule->CreateDynamicTool();
    if (setDepth)
        dataStructureModule->SetStackDepth(stackDepth);
    dataStructureModule->UnLock();
#endif
}

void AddressStreamDriver::SetUpTools() {
    // Check for which tools to use
    uint32_t doAddressRange;
    uint32_t doArielFrontend;
    uint32_t doCacheSimulation;
    uint32_t doEntropyRange;
    uint32_t doHardwarePrefetching;
    uint32_t doReuseDistance;
    uint32_t doScatterGatherLength;
    uint32_t doSpatialLocality;
    uint32_t doSpatialLocalityPerMemOp;
    if (parser->ReadEnvUint32("METASIM_ADDRESS_RANGE", &doAddressRange)){
        runAddressRange = (doAddressRange == 0) ? false : true;
    }
    if (parser->ReadEnvUint32("METASIM_ARIEL_FRONTEND", &doArielFrontend)){
        runArielFrontend = (doArielFrontend == 0) ? false : true;
    }
    if (parser->ReadEnvUint32("METASIM_CACHE_SIMULATION", &doCacheSimulation)){
        runCacheSimulation = (doCacheSimulation == 0) ? false : true;
    }
    if (parser->ReadEnvUint32("METASIM_ENTROPY_RANGE", &doEntropyRange)){
        runEntropyRange = (doEntropyRange == 0) ? false : true;
    }
    if (parser->ReadEnvUint32("METASIM_HWPF_SIMULATION", 
      &doHardwarePrefetching)){
        runHardwarePrefetching = (doHardwarePrefetching == 0) ? false : true;
    }
    if (parser->ReadEnvUint32("METASIM_REUSE_DISTANCE", &doReuseDistance)){
        runReuseDistance = (doReuseDistance == 0) ? false : true;
    }
    if (parser->ReadEnvUint32("METASIM_SG_LENGTH", &doScatterGatherLength)){
        runScatterLength = (doScatterGatherLength == 0) ? false : true;
    }
    if (parser->ReadEnvUint32("METASIM_SPATIAL_LOCALITY", &doSpatialLocality)){
        runSpatialLocality = (doSpatialLocality == 0) ? false : true;
    }
    if (parser->ReadEnvUint32("METASIM_SPATIAL_LOCALITY_MEMOP", 
      &doSpatialLocalityPerMemOp)){
        runSpatialLocalityPerMemOp = (doSpatialLocalityPerMemOp == 0) ? false :
          true;
    }

    // Check for which types of tools to use
    uint32_t doCodeCentric;
    uint32_t doDataCentric;
    if (parser->ReadEnvUint32("METASIM_CODE_CENTRIC", &doCodeCentric)){
        runCodeCentric = (doCodeCentric == 0) ? false : true;
    }
    if (parser->ReadEnvUint32("METASIM_DATA_CENTRIC", &doDataCentric)){
        runDataCentric = (doDataCentric == 0) ? false : true;
    }

    // Check that one is set
    if (!(runCodeCentric || runDataCentric)) {
        DISPLAY_ERROR << "Neither Data Centric not Code Centric run is set. "
          << "Set one. Exitting." << ENDL;
        exit(0);
    }

    // First add code-centric tools.
    if (runAddressRange && runCodeCentric) {
        tools->push_back(new AddressRangeTool());
    }

    if (runArielFrontend && runCodeCentric) {
        if (BuiltWithArielFrontend()) {
            tools->push_back(GENERATE_ARIEL_TOOL);
        } else {
            DISPLAY_ERROR << "No ariel frontend library included. "
              << "Unset Ariel frontend library tool. Exitting." << ENDL;
            exit(0);
        }
    }

    if (runCacheSimulation && runCodeCentric) {
        tools->push_back(new CacheSimulationTool());
    }

    if (runEntropyRange && runCodeCentric) {
        tools->push_back(new EntropyRangeTool());
    }

    if (runHardwarePrefetching) {
        if (BuiltWithEPATools()) {
            tools->push_back(GENERATE_PREFETCH_TOOL);
        } else {
            DISPLAY_ERROR << "No hardware prefetching library linked. "
              << "Unset Hardware prefetching library tool. Exitting." << ENDL;
            exit(0);
        }
    }

    if (runReuseDistance && runCodeCentric) {
        tools->push_back(new ReuseDistanceTool());
    }

    if (runScatterLength) {
        tools->push_back(new ScatterGatherLengthTool());
    }

    if (runSpatialLocality && runCodeCentric) {
        tools->push_back(new SpatialLocalityTool());
    }

    if (runSpatialLocalityPerMemOp) {
        if (BuiltWithEPATools()) {
            tools->push_back(GENERATE_SPATIAL_MEMOP_TOOL);
        } else {
            DISPLAY_ERROR << "No spatial locality per memop library linked. "
              << "Unset Spatial locality per memop library tool. Exitting." 
              << ENDL;
            exit(0);
        }
    }

    numCodeCentricTools = tools->size();

    // THEN add the data-centric tools.
    if (!BuiltWithDataStructureModule() && runDataCentric) {
        DISPLAY_ERROR << "No data structure module linked. "
          << "Unset Data Centric Libraries. Exitting." << ENDL;
        exit(0);
    }

    if (runAddressRange && runDataCentric) {
        tools->push_back(GENERATE_DATA_TOOL(DataCentricAddressRangeTool));
    }

    if (runSpatialLocality && runDataCentric) {
        tools->push_back(GENERATE_DATA_TOOL(DataCentricSpatialLocalityTool));
    }

    if (runCacheSimulation && runDataCentric) {
        tools->push_back(GENERATE_DATA_TOOL(DataCentricCacheSimulationTool));
    }

    if (runEntropyRange && runDataCentric) {
        tools->push_back(GENERATE_DATA_TOOL(DataCentricEntropyRangeTool));
    }

    if (runReuseDistance && runDataCentric) {
        tools->push_back(GENERATE_DATA_TOOL(DataCentricReuseDistanceTool));
    }

    uint32_t toolIndex = 0;
    for (vector<AddressStreamTool*>::iterator it = tools->begin(); it != 
      tools->end(); it++) {
        AddressStreamTool* currentTool = (*it);
        uint32_t handlersAdded = currentTool->CreateHandlers(
          GetNumMemoryHandlers(), parser);
        assert(handlersAdded > 0);
        numMemoryHandlers += handlersAdded;
        if (toolIndex < numCodeCentricTools)
            numCodeCentricMemoryHandlers += handlersAdded;
        toolIndex++;
    }
}

void AddressStreamDriver::ShutOffInstrumentationInAllBlocks() {
    // map of imageSequences -> set of blocks to shut off
    map<uint64_t, set<uint64_t>> allBlocks;
    for (set<uint64_t>::iterator it = liveMemoryAccessInstPointKeys->begin();
      it != liveMemoryAccessInstPointKeys->end(); it++) {
        uint32_t blockID = GET_BLOCKID(*it);
        uint32_t imageSequence = GET_IMAGEID(*it);
        (allBlocks[imageSequence]).insert(blockID);
//        allBlocks.insert(GET_BLOCKID(*it));
    }
    for (map<uint64_t, set<uint64_t>>::iterator it = allBlocks.begin();
      it != allBlocks.end(); it++) {
        uint32_t imageSequence = (*it).first;
        image_key_t imageID = allData->GetImageId(imageSequence);
        ShutOffInstrumentationInBlocks(allBlocks[imageSequence], imageID);

    }
}

// Not thread-safe! For performance, thread suspension should happen outside 
// this function
void AddressStreamDriver::ShutOffInstrumentationInBlock(uint64_t blockID, 
  uint64_t imageSequence) {

    // Note: unique keys generated in InitializeDynamicInstrumentation
    set<uint64_t> keysToRemove;
    uint64_t kcheck = GENERATE_UNIQUE_KEY(blockID, imageSequence, 
      PointType_buffercheck);
    uint64_t kinc = GENERATE_UNIQUE_KEY(blockID, imageSequence, 
      PointType_bufferinc);
    uint64_t kfill = GENERATE_UNIQUE_KEY(blockID, imageSequence, 
      PointType_bufferfill);

    // If this key is not active, then done
    if (liveMemoryAccessInstPointKeys->count(kfill) == 0){
        return;
    }

    // Otherwise, remove the instrumentation for this block
    keysToRemove.insert(kcheck);
    keysToRemove.insert(kinc);
    keysToRemove.insert(kfill);
    
    dynamicPoints->SetDynamicPoints(keysToRemove, false);
    liveMemoryAccessInstPointKeys->erase(kfill);

}

void AddressStreamDriver::ShutOffInstrumentationInBlocks(set<uint64_t>& blocks,
  image_key_t iid, bool suspend) {
    // Make sure only one thread is executing this code
    if (suspend) {
        allData->ReadLock();
        sampler->WriteLock();
        SuspendAllThreads(allData->livethreads.size(),
          allData->livethreads.begin(), allData->livethreads.end());
    }

    uint64_t imageSequence = (uint32_t)allData->GetImageSequence(iid, false);
    
    for (set<uint64_t>::iterator it = blocks.begin(); it != blocks.end(); 
      it++) {
        uint64_t blockID = *it;
        ShutOffInstrumentationInBlock(blockID, imageSequence);
    }

    if (suspend) {
        ResumeAllThreads();
        sampler->UnLock();
        allData->UnLock();
    }
}

void AddressStreamDriver::ShutOffInstrumentationInMaxedGroups(image_key_t iid, 
  thread_key_t tid, bool suspend) {

    bool lock = suspend;

    // Thread-safe call
    AddressStreamStats* stats = (AddressStreamStats*)allData->GetData(iid, 
      tid, lock);

    // Make sure group counters are up to date
    for(uint32_t i = 0; i < (stats->BlockCount); i++) {
        uint32_t idx = i;
        if (stats->Types[i] == CounterType_instruction) {
            idx = stats->Counters[i];
        }
        uint64_t blocksGroupId = stats->GroupIds[i]; 
        uint64_t blockCount = stats->Counters[idx];
        if(stats->GroupCounters[blocksGroupId] < blockCount) {
            stats->GroupCounters[blocksGroupId] = blockCount;
        }
    }

    // Can't combine this with above because a later block could cause 
    // group to exceed max
    set<uint64_t> blocksToRemove;
    if (suspend) {
        allData->ReadLock();
        sampler->WriteLock();
        SuspendAllThreads(allData->livethreads.size(),
          allData->livethreads.begin(), allData->livethreads.end());
    }
    
    for (set<uint64_t>::iterator it = liveMemoryAccessInstPointKeys->begin();
      it != liveMemoryAccessInstPointKeys->end(); it++) {
        uint64_t blockID = GET_BLOCKID(*it);
        image_key_t imageID = allData->GetImageId(GET_IMAGEID(*it), false);
        // Only shut off for current image
        if (imageID != iid)
            continue;
        // If max count is reached, we will remove this block
        uint64_t blocksGroupId = stats->GroupIds[blockID]; 
        if (sampler->ExceedsAccessLimit(stats->GroupCounters[blocksGroupId], 
          false))
            blocksToRemove.insert(blockID);            
    }

    // Only call this if there are blocks to remove since it will suspend 
    // threads
    if (blocksToRemove.size() > 0)
        ShutOffInstrumentationInBlocks(blocksToRemove, iid, false);

    if (suspend) {
        ResumeAllThreads();
        sampler->UnLock();
        allData->UnLock();
    }
}

void AddressStreamDriver::UnpauseApplicationWrappers() {
    UNPAUSE_MODULE(dataStructureModule);
}

void AddressStreamDriver::UnLockDSM(bool lock) {
    UNLOCK(dataStructureModule, lock);
}

void AddressStreamDriver::WriteLockDSM(bool lock) {
    WRITELOCK(dataStructureModule, lock);
}

// For testing
AddressStreamTool* AddressStreamDriver::GetTool(uint32_t index) {
    return tools->at(index);
}

void AddressStreamDriver::SetParser(StringParser* p) {
    if (parser != NULL)
        delete parser;
    parser = p;
}

void AddressStreamDriver::SetSampler(SamplingMethod* s) {
    if (sampler != NULL)
        delete sampler;
    sampler = s;
}

void GetBufferIds(BufferEntry* b, image_key_t* i){
    *i = b->imageid;
}
