#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <InstrumentationCommon.h>

#define USING_MPI_WRAPPERS
#define USING_CSTD_WRAPPERS

FILE* logfile;
int64_t timerstart;
int64_t timerstop;
char* appName;
int32_t* currentSiteIndex;
char** fileNames;
int32_t* lineNumbers;

inline int starttimer(){
    timerstart = readtsc();
}

inline int stoptimer(){
    timerstop = readtsc();
}

inline int64_t gettimer(){
    return timerstop - timerstart;
}

#define printtimer(__stream) PRINT_INSTR(__stream, "timer value (in cycles): %lld", gettimer())

// do any initialization here
// NOTE: on static-linked binaries, calling any functions from here will cause some problems
int32_t initwrapper(int32_t* indexLoc, char** fNames, int32_t* lNum){
    // at each call site we will put the index of the originating point in this location
    currentSiteIndex = indexLoc;

    fileNames = fNames;
    lineNumbers = lNum;

    logfile = stdout;
}

// do any cleanup here
int32_t finishwrapper(){
}

#ifdef USING_CSTD_WRAPPERS
#include <IOWrappers.c>
#endif // USING_CSTD_WRAPPERS

#ifdef USING_MPI_WRAPPERS
#include <MPIWrappers.c>
#endif // USING_MPI
