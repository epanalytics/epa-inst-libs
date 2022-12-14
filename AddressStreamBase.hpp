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

#ifndef _AddressStreamBase_hpp_
#define _AddressStreamBase_hpp_

#include <vector>
#include <string>
#include <AddressStreamStats.hpp>

template <class T> class DataManager;
class MemoryStreamHandler;
class SamplingMethod;
class StringParser;

#define KILO (1024)
#define MEGA (KILO*KILO)
#define GIGA (MEGA*KILO)

class AddressStreamTool {
  protected:
    std::vector<MemoryStreamHandler*> handlers;
    int32_t indexInStats = -1;
  public:
    AddressStreamTool() : indexInStats(-1) {};
    virtual ~AddressStreamTool();
    virtual void AddNewHandlers(AddressStreamStats* stats) = 0;
    virtual void AddNewStreamStats(AddressStreamStats* stats) = 0;
    virtual uint32_t CreateHandlers(uint32_t, StringParser*) = 0;
    virtual void FinalizeTool(DataManager<AddressStreamStats*>*, 
      SamplingMethod*) = 0;
};

class StreamStats {
  public:
    virtual uint64_t GetAccessCount(uint32_t memid) = 0;
    virtual bool Verify() = 0;
};

// Note: User required to check if limit is hit
class SamplingMethod {
  private:
    uint32_t AccessLimit;
    uint32_t SampleOn;
    uint32_t SampleOff;
    uint64_t AccessCount;

    // Samplers can be shared by threads -- use a r/w lock so that 
    // multiple threads can read at once
    pthread_rwlock_t sampling_rwlock;
    pthread_rwlockattr_t sampling_rwlock_attr;

    bool CurrentlySampling(uint64_t count);

  public:
    SamplingMethod(uint32_t limit, uint32_t on, uint32_t off);
    virtual ~SamplingMethod();

    virtual bool CurrentlySampling();
    virtual bool ExceedsAccessLimit(uint64_t count);
    virtual uint64_t GetAccessCount() { return AccessCount; }
    uint64_t GetAccessLimit() { return AccessLimit; }
    virtual double GetSamplingFrequency();
    uint32_t GetSampleOn() { return SampleOn; }
    uint32_t GetSampleOff() { return SampleOff; }
    void IncrementAccessCount(uint64_t count);
    virtual bool SwitchesMode(uint64_t count);
    void Print();

    bool ReadLock();
    bool UnLock();
    bool WriteLock();
};

// DFP and other interesting memory things extend this class.
class MemoryStreamHandler {
  protected:
    pthread_mutex_t mlock;
  public:
    MemoryStreamHandler();
    virtual ~MemoryStreamHandler();

    virtual void Print(std::ofstream& f) = 0;
    virtual uint32_t Process(void* stats, BufferEntry* access) = 0;
    // Number of addresses that appeared but aren't processed
    virtual void SkipAddresses(uint32_t numToSkip) {};
    virtual bool Verify() = 0;
    bool Lock();
    bool UnLock();
    bool TryLock();

};

// Common string functions
class StringParser {
  public:
    StringParser() {};
    virtual ~StringParser() {};

	virtual char* GetEnv(const char* var);
    virtual bool IsEmptyComment(std::string str);
    virtual bool ParseInt32(std::string token, int32_t* value, int32_t min);
    virtual bool ParsePositiveInt32(std::string token, uint32_t* value);
    virtual bool ReadEnvUint32(std::string name, uint32_t* var);
    virtual char ToLowerCase(char c);
};

class Randomizer {
  public:
    Randomizer() {};
    virtual ~Randomizer() {}    

    virtual uint32_t RandomInt(uint32_t max);
};

class EasyHash {
  private:
    pebil_map_type<uint32_t, uint32_t>* internal_map; //key to running count

  public:
    EasyHash();
    ~EasyHash();
    bool contains(uint32_t key); //checks if key is currently in map
    void add(uint32_t key, uint32_t valueToAdd);//takes the key, and adds to the running count
                                           //if no key is present, creates and sets the running
                                           //count to the passed in valueToAdd
    uint32_t get(uint32_t key); //returns running count returns 0 if key is not found
};

class NestedHash {
  private:
    pebil_map_type<uint32_t, EasyHash*>* internal_hash; //set to EasyHash that has line and
                                                        //running count

  public:
    NestedHash();
    ~NestedHash();
    bool contains(uint32_t set, uint32_t line); //checks if the set and line combo is in the map
    //vvv Takes the set and line combination and adds to the running count,
    //if no key is found, creates the key and sets valueToAdd as the starting running count
    void put(uint32_t set, uint32_t line, uint32_t valueToAdd);
    uint32_t get(uint32_t set, uint32_t line);//returns running count for set and line
                                              //returns 0 if key is not found
};

#endif /* _AddressStreamBase_hpp_ */

