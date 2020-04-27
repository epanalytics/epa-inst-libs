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

#ifndef _CacheSimulation_hpp_
#define _CacheSimulation_hpp_

#include <AddressStreamBase.hpp>
#include <string>
#include <fstream>

#define DEFAULT_CACHE_FILE "instcode/CacheDescriptions.txt"

#define INVALID_CACHE_LEVEL (0xffffffff)

class CacheStructureHandler;
class CacheStats;

enum CacheLevelType {
    CacheLevelType_Undefined,
    CacheLevelType_InclusiveLowassoc,
    CacheLevelType_NonInclusiveLowassoc,
    CacheLevelType_InclusiveHighassoc,
    CacheLevelType_Total
};

enum ReplacementPolicy {
    ReplacementPolicy_Undefined,
    ReplacementPolicy_trulru,
    ReplacementPolicy_nmru,
    ReplacementPolicy_random,
    ReplacementPolicy_direct,
    ReplacementPolicy_Total
};

static const char* ReplacementPolicyNames[ReplacementPolicy_Total] = {
    "undefined",
    "truelru",
    "nmru",
    "random",
    "direct"
};

struct EvictionInfo {
    uint64_t addr;
    uint32_t level;
    uint32_t setid;
    uint32_t lineid;
    bool coldMiss;
};

struct LevelStats {
    uint64_t hitCount;
    uint64_t missCount;
    uint64_t loadCount; 
    uint64_t storeCount;
};

class MainMemory {
public:
    MainMemory(uint32_t setSize, uint32_t numOfLines, uint32_t lineSize); 
    MainMemory(MainMemory& mem);
    ~MainMemory();

    //Below HashMaps are used to keep count of loads and store to and from main memory
    NestedHash* writeOutsMap = nullptr; //keyed by set and line
    NestedHash* readInsMap = nullptr; //keyed by set and line

    EasyHash* dirOutsMap = nullptr; //if using a direct map last level cache,
    EasyHash* dirInsMap = nullptr; //should hopefully cut down even more memory usage

    uint32_t GetLoads(); //loops through all of readIns and gets a total sum
    uint32_t GetStores(); //loops through all of writeOuts and gets a total sum

    uint32_t numOfSets; //based on Last Level Cache
    uint32_t numOfLinesInSet;
    uint32_t sizeOfLine;
};

class CacheSimulationTool : public AddressStreamTool {
  public:
    CacheSimulationTool() : AddressStreamTool() {}
    virtual void AddNewHandlers(AddressStreamStats* stats);
    virtual void AddNewStreamStats(AddressStreamStats* stats);
    void CacheSimulationFileName(AddressStreamStats* stats, std::string& oFile);
    virtual uint32_t CreateHandlers(uint32_t index, StringParser* parser);
    virtual void FinalizeTool(DataManager<AddressStreamStats*>* AllData,
      SamplingMethod* Sampler);
    std::string GetCacheDescriptionFileName() { return CacheDescriptionFile; }
    void GetAndSetCacheDescriptionFile(StringParser* parser);
    virtual void HandleEnvVariables(StringParser* parser);
    bool IsLoadStoreLogging() { return loadStoreLogging; }
    void LogFileName(AddressStreamStats* stats, std::string& oFile);
    virtual uint32_t ReadCacheDescription(std::istream& cacheStream, 
      StringParser* parser);

    void SetCacheDescriptionFileName(std::string f) {CacheDescriptionFile = f;}

  protected:
    std::string CacheDescriptionFile;
    uint32_t MinimumHighAssociativity = 256;
    bool loadStoreLogging = false;
    uint32_t DirtyCacheHandling = 0;
    void PrintApplicationHeader(std::ofstream& file, 
      DataManager<AddressStreamStats*>* AllData, SamplingMethod* Sampler, 
      uint64_t totalMemop, uint64_t sampledCount);
    void PrintPerBlockCacheSimData(std::ofstream& file,
      DataManager<AddressStreamStats*>* AllData);
    void PrintSysidInfo(std::ofstream& file, CacheStats* c, std::set<image_key_t>::iterator iit);
    void PrintThreadidInfo(std::ofstream& file, thread_key_t thread, 
      DataManager<AddressStreamStats*>* AllData);
};


class CacheStats : public StreamStats {
public:
    uint32_t LevelCount;
    uint32_t SysId;
    LevelStats** Stats; // indexed by [memid][level]
    MainMemory** mainMemoryStats; // indexed by [memid]
    uint32_t Capacity;
    CacheStats(uint32_t lvl, uint32_t sysid, uint32_t capacity, uint32_t 
      hybridCache);
    ~CacheStats();
    void InitMainMemoryStats(CacheStructureHandler* handler);

    bool HasMemId(uint32_t memid);
    void ExtendCapacity(uint32_t newSize);
    void NewMem(uint32_t memid);

    void Hit(uint32_t memid, uint32_t lvl);

    void Miss(uint32_t memid, uint32_t lvl);

    void Hit(uint32_t memid, uint32_t lvl, uint32_t cnt);

    void Miss(uint32_t memid, uint32_t lvl, uint32_t cnt);

    void Load(uint32_t memid,uint32_t lvl);
    void Load(uint32_t memid, uint32_t lvl, uint32_t cnt);
   
    void Store(uint32_t memid,uint32_t lvl);
    void Store(uint32_t memid, uint32_t lvl, uint32_t cnt);
 
    uint64_t GetLoads(uint32_t memid, uint32_t lvl);
    uint64_t GetLoads(uint32_t lvl);
    
    uint64_t GetStores(uint32_t memid, uint32_t lvl);
    uint64_t GetStores(uint32_t lvl);    

    static float GetHitRate(LevelStats* stats);
    static float GetHitRate(uint64_t hits, uint64_t misses);

    uint64_t GetHits(uint32_t memid, uint32_t lvl);

    uint64_t GetHits(uint32_t lvl);
    
    uint64_t GetMisses(uint32_t memid, uint32_t lvl);

    uint64_t GetMisses(uint32_t lvl);

    LevelStats* GetLevelStats(uint32_t memid, uint32_t lvl);
    uint64_t GetAccessCount(uint32_t memid);
    float GetHitRate(uint32_t memid, uint32_t lvl);
    float GetCumulativeHitRate(uint32_t memid, uint32_t lvl);

    bool Verify();
};

#define USES_MARKERS(__pol) (__pol == ReplacementPolicy_nmru)
#define CacheLevel_Init_Interface uint32_t lvl, uint32_t sizeInBytes, uint32_t assoc, uint32_t lineSz, ReplacementPolicy pol, uint32_t loadStore, uint32_t dirtyCache
#define CacheLevel_Init_Arguments lvl, sizeInBytes, assoc, lineSz, pol, loadStore, dirtyCache

struct history {
    uint32_t prev;
    uint32_t next;
};

class CacheLevel {
protected:

    CacheLevelType type;

    uint32_t level;
    uint32_t levelCount;
    uint32_t size;
    uint32_t associativity;
    uint32_t linesize;
    ReplacementPolicy replpolicy;

    uint32_t countsets;
    uint32_t linesizeBits;

    uint64_t** contents = nullptr;
    bool**  dirtystatus = nullptr;
    uint32_t* recentlyUsed = nullptr;
    history** historyUsed = nullptr;
    bool toEvict;

	bool loadStoreLogging;
	uint32_t dirtyCacheHandling;

public:
    std::vector<uint64_t>* toEvictAddresses = nullptr;
    CacheLevel();
    virtual ~CacheLevel();

    uint32_t GetLevelCount() { return levelCount;}
    uint32_t SetLevelCount(uint32_t inpLevelCount) { return levelCount = 
      inpLevelCount; }
    virtual CacheLevelType GetType() { return type; }
    ReplacementPolicy GetReplacementPolicy() { return replpolicy; }
    virtual uint32_t GetLevel() { return level; }
    virtual uint32_t GetSizeInBytes() { return size; }
    virtual uint32_t GetAssociativity() { return associativity; }
    virtual uint32_t GetSetCount() { return countsets; }
    virtual uint32_t GetLineSize() { return linesize; }
    virtual bool GetLoadStoreLog() { return IsLoadStoreLogging(); }
    virtual uint32_t GetDirtyCacheHandle() { return dirtyCacheHandling; }
    uint64_t CountColdMisses();
    virtual bool IsLoadStoreLogging() { return loadStoreLogging; }

    void Print(std::ofstream& f, uint32_t sysid);

    // re-implemented by Exclusive/InclusiveCacheLevel
    virtual bool IsExclusive() { return false; }
    virtual uint32_t Process(CacheStats* stats, uint32_t memid, uint64_t addr, 
      uint64_t loadstoreflag, bool* anyEvict, EvictionInfo* info);
    virtual uint32_t EvictProcess(CacheStats* stats, uint32_t memid, uint64_t 
      addr, uint64_t loadstoreflag, void* info);    

    virtual void EvictDirty(CacheStats* stats, CacheLevel** levels, uint32_t 
      memid, void* info); // void* info is needed since eventually 'Process' 
      // needs to be called! 
    virtual bool GetEvictStatus();

    std::vector<uint64_t>* passEvictAddresses() { return toEvictAddresses;}

protected:

    uint64_t GetStorage(uint64_t addr);
    uint64_t GetAddress(uint64_t store);
    uint32_t GetSet(uint64_t addr);
    uint32_t LineToReplace(uint32_t setid);
    bool MultipleLines(uint64_t addr, uint32_t width);

    void MarkUsed(uint32_t setid, uint32_t lineid,uint64_t loadstoreflag);

    // re-implemented by HighlyAssociativeCacheLevel
    virtual bool Search(uint64_t addr, uint32_t* set, uint32_t* lineInSet);
    virtual uint64_t Replace(uint64_t addr, uint32_t setid, uint32_t lineid,uint64_t loadstoreflag);

    virtual const char* TypeString() = 0;
    virtual void Init (CacheLevel_Init_Interface);
    
   // Both store and lineid is being sent since while calling these methods we do not make distinction as whether the object belongs to CacheLevel or HighlyAssociateCacheLevel
    virtual void SetDirty(uint32_t setid, uint32_t lineid,uint64_t store);
    virtual void ResetDirty(uint32_t setid, uint32_t lineid,uint64_t store);
    virtual bool GetDirtyStatus(uint32_t setid, uint32_t lineid,uint64_t store);
};

class InclusiveCacheLevel : public virtual CacheLevel {
public:
    InclusiveCacheLevel() {}
    ~InclusiveCacheLevel() {}

    virtual void Init (CacheLevel_Init_Interface){
        CacheLevel::Init(CacheLevel_Init_Arguments);
        type = CacheLevelType_InclusiveLowassoc;
    }
    bool IsExclusive() { return false; }
    virtual const char* TypeString() { return "inclusive"; }
};

class NonInclusiveCacheLevel : public virtual CacheLevel {
public:
    uint32_t FirstExclusive;
    uint32_t LastExclusive;

    NonInclusiveCacheLevel() {}
    virtual void Init(CacheLevel_Init_Interface){
        CacheLevel::Init(CacheLevel_Init_Arguments);
        type = CacheLevelType_NonInclusiveLowassoc;
    }
    bool IsExclusive() { return false; }
    uint32_t Process(CacheStats* stats, uint32_t memid, uint64_t addr, uint64_t
      loadstoreflag, bool* anyEvict, EvictionInfo* info);
    virtual const char* TypeString() { return "exclusive"; }
};

class HighlyAssociativeCacheLevel : public virtual CacheLevel {
protected:
    pebil_map_type <uint64_t, uint32_t>** fastcontents;
    pebil_map_type <uint64_t, bool>** fastcontentsdirty;
public:
    HighlyAssociativeCacheLevel() {}
    ~HighlyAssociativeCacheLevel();

    bool Search(uint64_t addr, uint32_t* set, uint32_t* lineInSet);
    uint64_t Replace(uint64_t addr, uint32_t setid, uint32_t lineid,uint64_t loadstoreflag);
    virtual void Init (CacheLevel_Init_Interface);

   // Both store and lineid is being sent since while calling these methods we do not make distinction as whether the object belongs to CacheLevel or HighlyAssociateCacheLevel
   void SetDirty(uint32_t setid, uint32_t lineid,uint64_t store);
   void ResetDirty(uint32_t setid, uint32_t lineid,uint64_t store);
   bool GetDirtyStatus(uint32_t setid, uint32_t lineid,uint64_t store);  
};

class HighlyAssociativeInclusiveCacheLevel : public InclusiveCacheLevel, public HighlyAssociativeCacheLevel {
public:
    HighlyAssociativeInclusiveCacheLevel() {}
    virtual void Init (CacheLevel_Init_Interface){
        InclusiveCacheLevel::Init(CacheLevel_Init_Arguments);
        HighlyAssociativeCacheLevel::Init(CacheLevel_Init_Arguments);
        type = CacheLevelType_InclusiveHighassoc;
    }
    const char* TypeString() { return "inclusive_H"; }
};

class CacheStructureHandler : public MemoryStreamHandler {
  protected: 

    uint32_t sysId;
    uint32_t levelCount;

    CacheLevel** levels;
    //std::string description;

    uint32_t MinimumHighAssociativity = 256;
    bool loadStoreLogging = false;  // Record loads/stores with cache activity
    uint32_t DirtyCacheHandling = 0;

    bool isInitialized = false;

    uint64_t hits;
    uint64_t misses;
    //uint64_t AddressRangesCount;
    StringParser* parser;
    std::vector<uint64_t>* toEvictAddresses;
    uint32_t processAddress(void* stats, uint64_t address, uint64_t memseq, 
      uint8_t loadstoreflag);

  public:      
    // note that this doesn't contain any stats gathering code. that is done at
    // the thread level and is therefore done in ThreadData

    CacheStructureHandler();
    CacheStructureHandler(CacheStructureHandler& h);
    ~CacheStructureHandler();
    virtual bool Init(std::string desc, uint32_t MinimumHighAssociativity, bool 
      doLoadStore, uint32_t DirtyCacheHandling);

    CacheLevel* GetCacheLevel(uint32_t lvl) { return levels[lvl]; } //0-indexed
    uint32_t GetNumberOfCacheLevels() { return levelCount; }
    uint32_t GetSysId() { return sysId; }

    bool IsInitialized() { return isInitialized; }
    bool IsLoadStoreLogging() { return loadStoreLogging; }

    CacheLevel* ParseCacheLevelTokens(std::stringstream& tokenizer, 
      uint32_t levelId, uint32_t* firstExcl);
    void Print(std::ofstream& f);
    uint32_t Process(void* stats, BufferEntry* access);
    void SetParser(StringParser* p) { parser = p; }
    bool Verify();

    uint64_t GetHits(){return hits;}
    uint64_t GetMisses(){ return misses;} 

    //bool CheckRange(CacheStats* stats,uint64_t addr,uint64_t loadstoreflag,uint32_t memid); //, uint32_t* set, uint32_t* lineInSet);    
    //void ExtractAddresses();

    /* For testing only */
    void SetCacheLevel(uint32_t lvl, CacheLevel* c) { levels[lvl] = c; }
    void SetLoadStoreLogging(bool b) { loadStoreLogging = b; }
    void SetNumberOfCacheLevels(uint32_t n) { levelCount = n; }
};


#endif /* _CacheSimulation_hpp_ */
