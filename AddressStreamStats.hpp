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

#ifndef _AddressStreamStats_hpp_
#define _AddressStreamStats_hpp_

//TODO see if below was actually needed
//#include <Metasim.hpp>
//#define debug(...) __VA_ARGS__
#define debug(...)

#define LOAD 1
#define STORE 0

#define NORMAL 0
#define SWPF   1

typedef uint64_t image_key_t;
typedef pthread_t thread_key_t;

enum EntryType: uint8_t {
  MEM_ENTRY = 0,
  VECTOR_ENTRY,
  EPAX_VECTOR_ENTRY,
  EPAX_INDIRECT_ENTRY,
  EntryType_Total
};

// pebil only
struct VectorAddress { 
    uint32_t indexVector[16];
    uint8_t  scale;
    uint64_t base;
    uint64_t mask;
    uint32_t  numIndices;
};

// EPAX-only. For use with SVE memops that do a contiguous memory access.
// The addressing modes are called:
// 1. Scalar plus immediate
// 2. Scalar plus scalar
// Since this is a contiguous memory access, you only need the first address
// accessed and how much it accesses. Use the number of elements to process
// each accessed address. Use the number of elements to determine which mask 
// bits to use.
struct EPAXVectorAddress {
    uint64_t memAddress;    // First Address in contiguous mem access
    uint64_t sizeOfAccess;  // Size of mem access in bits
    uint16_t numElements;   // Number of addresses accessed
    // The value of the predicate register
    // Allocate maximum length of predicate register == SVE VL / 64
    // (1 bit for each byte of the SVE Vector Length)
    // Max SVE VL == 2048; Max P Reg length == 2048 / 64 == 32
    uint8_t predReg[32];
};

// EPAX-only. For use with SVE memops that do a scatter/gather memory access.
// The addressing modes are called:
// 1. Scalar plus vector
// 2. Vector plus immediate
// These require calculating the address based on the ARM documentation and
// the given values. Scalar + vector uses a base address, may sign extend and/or
// shift the given index value and then add the values together. Vector +
// immediate does not use a base address and instead adds an immediate to each
// value in the given vector
struct EPAXIndirectAddress {
    uint64_t baseAddress;   // scalar + vec: Value of Xn OR 0
    bool doesExtension;     // scalar + vec: do a signed or unsigned extension
    bool signedExtend;      // scalar + vec: signedExtend or unsignedExtend
    uint8_t shiftAmount;    // scalar + vec: How much to shift index
    uint64_t immediate;     // vector + imm: offset * mbytes OR 0
                            // (see ARM documentation)
    uint16_t numElements;   // Number of addresses accessed
    // The value of the predicate register
    // Allocate maximum length of predicate register in bytes == SVE VL / 64
    // (1 bit for each byte of the SVE Vector Length)
    // Max SVE VL == 2048; Max P Reg length == 2048 / 64 == 32
    uint8_t predReg[32];
    // The value of the Z register in the memory operand (base or index vector)
    // Allocate maximum length of Z register in bytes == SVE VL / 8
    // Max SVE VL == 2048; Max Z Reg length == 2048 / 8  == 256
    uint8_t baseVector[256];
};

typedef struct BufferEntry_s {
    enum EntryType  type;
    uint8_t         swprefetchflag;  // Is a software prefetch op
    uint8_t         loadstoreflag;   // Dirty Caching
    uint64_t        imageid;         // Multi-image
    uint64_t        memseq;          // identifies memop in image
    union {
        uint64_t address;        // value simulated
        struct VectorAddress vectorAddress;
#ifdef  EPAX_INST_TOOL
        struct EPAXVectorAddress epaxVectorAddress;
        struct EPAXIndirectAddress epaxIndirectAddress;
#endif
    };
    //uint64_t    threadid;        // Error-checking
} BufferEntry;
#define __buf_current  vectorAddress.base
#define __buf_oldPosition  vectorAddress.mask
#define __buf_capacity memseq

class StreamStats;
class MemoryStreamHandler;
class ReuseDistance;

typedef struct AddressStreamStats_s {
    // memory buffer
    BufferEntry* Buffer;

    // metadata
    thread_key_t threadid;
    image_key_t imageid;
    bool FirstImage;    // Set to true if image is first image
    bool Initialized;   // Set to false when created by thread
    bool PerInstruction;
    bool LoopInclusion; // when terminating sampling for a block,
                        // do this for all blocks within the loop
                        // Note: includes all other blocks in the loop
    bool Master;        // Master image?
    uint32_t SVEVectorLength;  // Used only by EPAX
    uint32_t Phase;
    uint32_t AllocCount;
    uint32_t BlockCount;
    uint32_t GroupCount;
    uint32_t MemopCount;
    char* Application;
    char* Extension;

    // per-memop data
    uint64_t* BlockIds;   // Indices into per-block data, like counter

    // per-block data
    CounterTypes* Types; // If Counter is a count or index to a count
    uint64_t* Counters;
    uint32_t* MemopsPerBlock;
    // True for blocks where dynamic memops cannot be determinied at runtime
    // (e.g. memops with masks)
    bool* HasNonDeterministicMemop;
    char** Files;
    uint32_t* Lines;
    char** Functions;
    uint64_t* Hashes;
    uint64_t* Addresses;
    uint64_t* GroupIds;
    StreamStats** Stats; // indexed by handler
    MemoryStreamHandler** Handlers;
    ReuseDistance** RHandlers;

    // per-group data
    uint64_t* GroupCounters;

    // run data
    bool* UsedSlicerPtr;    // Pointer to global in Library
    uint64_t maxNumAddresses;
    uint64_t* addressesForProcessing;

} AddressStreamStats;

#define BUFFER_ENTRY(__stats, __n) (&(__stats->Buffer[__n+1]))
#define BUFFER_CAPACITY(__stats) (__stats->Buffer[0].__buf_capacity)
#define BUFFER_CURRENT(__stats) (__stats->Buffer[0].__buf_current)


#endif // _AddressStreamStats_hpp_
