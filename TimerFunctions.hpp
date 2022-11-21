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

#ifndef _TimerFunctions_hpp_
#define _TimerFunctions_hpp_
#include <sys/time.h>
#include <string>
using namespace std;


typedef struct FunctionTimers_s {
    bool master;
    char* application;
    char* extension;
    uint64_t functionCount;
    char** functionNames;
    uint64_t* functionHashes;
    uint64_t* functionTimerAccum;
    uint64_t* functionEntryCounts;
    uint32_t* functionShutoff;
    uint64_t* functionTimerLast;
    
    uint32_t* inFunction;
    uint64_t appTimeStart;
    struct timeval appTimeOfDayStart;
    bool sanitize = false;
} FunctionTimers;

#endif
#define KILO (1024)
#define MEGA (KILO*KILO)
#define GIGA (MEGA*KILO)
static char ToLowerCase(char c);
static bool ParsePositiveInt32(string token, uint32_t* value);
static bool ParseInt32(string token, uint32_t* value, uint32_t min);
static bool ParsePositiveInt32Hex(string token, uint32_t* value);
static bool ReadEnvUint32(string name, uint32_t* var);