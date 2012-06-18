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

#ifndef _CacheSim_h_
#define _CacheSim_h_

#include <stdint.h>

#ifdef PER_SET_RECENT
  #define MOSTRECENT(__setIdx) mostRecent[__setIdx]
#else
  #define MOSTRECENT(__setIdx) mostRecent
#endif

#define __MAX_LEVEL                      3
#define __MAX_LINEAR_SEARCH_ASSOC       64
#define __MAX_STRING_SIZE             1024


typedef uint32_t     Attribute_t;  
typedef uint64_t     Counter_t;    
#ifdef METASIM_32_BIT_LIB
  typedef uint32_t     Address_t;
#else
  typedef uint64_t     Address_t;
#endif

typedef enum {
    repl_lru = 0,
    repl_ran,
    repl_dir,
    repl_lru_vc,
    Total_ReplacementPolicy
} ReplacementPolicy;

typedef enum {
    inclusive_cache,
    victim_cache,
    prediction_cache
} CacheType;

#define IS_REPL_POLICY_VC(__policy)  (__policy == repl_lru_vc)
#define IS_REPL_POLICY_LRU(__policy) ((__policy == repl_lru) || (__policy == repl_lru_vc))
#define IS_REPL_POLICY_RAN(__policy) (__policy == repl_ran)
#define IS_REPL_POLICY_DIR(__policy) (__policy == repl_dir)


#define CACHE_LINE_INDEX(__addr,__bits) (__addr >> __bits)
#define __L1_CACHE_LEVEL 0

typedef enum {
    number_of_sets = 0,
    set_associativity,
    line_size_in_bits,
    replacement_policy,
    cache_type,
    Total_CacheAttributes
} CacheAttributes;

typedef enum {
    cache_hit = 0,
    cache_miss,
    Total_AccessStatus
} AccessStatus;
    
typedef struct {
    Address_t  key;
    uint16_t   value;
    uint16_t   next;
} HashEntry;

typedef struct {
    Attribute_t attributes[Total_CacheAttributes];
    Counter_t   hitMissCounters[Total_AccessStatus];
#ifdef PER_SET_RECENT
    uint16_t*   mostRecent;
#else
    uint16_t    mostRecent;
#endif
    Address_t*  content;
    HashEntry*  highAssocHash;
// Prediction Cache info, may want to move this elsewhere
    Address_t trainer;
    Address_t stream;
    uint8_t   nFetched;
    uint8_t   fetchDistance;
} Cache;

typedef struct {
    uint32_t  index;
    uint8_t   levelCount;
    Cache     levels[__MAX_LEVEL]; 
    uint8_t   isVictimCacheHierarchy;
} MemoryHierarchy;

extern void initCaches(MemoryHierarchy* systems, uint32_t systemCount);
extern void processVictimCache(Address_t toFind, Address_t* victim, uint8_t* nVictims, AccessStatus* status, Cache* cache);
extern void processInclusiveCache(Address_t toFind, Address_t* victim, uint8_t* nvictims, AccessStatus* status, Cache* cache);

#endif

