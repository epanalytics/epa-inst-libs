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

#ifndef _ArielFrontend_hpp_
#define _ArielFrontend_hpp_

#include <sst/core/sst_config.h>


#include <AddressStreamBase.hpp>
#include <string>

namespace SST {
  namespace ArielComponent {
    class ArielTunnel;
  }

  namespace Core {
    namespace Interprocess {
      template <typename TunnelType> class SHMChild;
    }
  }
}

class ArielFrontendHandler;

class ArielFrontendTool : public AddressStreamTool {
  protected:
    std::string shmemName;
    uint32_t traceRank; // Ariel picks a rank to collect data from
    ArielFrontendHandler* tunnelCreator;
    std::vector<ArielFrontendHandler*> tunnelUsers;
    bool usesMPI;
  public:
    ArielFrontendTool() : AddressStreamTool(), shmemName(""), traceRank(0),
      tunnelCreator(NULL), usesMPI(false) {}
    virtual void AddNewHandlers(AddressStreamStats* stats);
    virtual void AddNewStreamStats(AddressStreamStats* stats);
    virtual uint32_t CreateHandlers(uint32_t index, StringParser* parser);
    virtual void FinalizeTool(DataManager<AddressStreamStats*>* AllData,
      SamplingMethod* Sampler);
    virtual void NotifyArielOutputStats(AddressStreamStats* stats);
    virtual void NotifyDoneMPIInit();
};

class ArielStats : public StreamStats {
private:
    uint64_t threadId;
    bool* isDP;
    bool* isFP;
    uint32_t* sizeInBytes;
    uint64_t* instPtr;

public:

    ArielStats(uint32_t threadSeq);
    virtual ~ArielStats();

    uint64_t GetAccessCount(uint32_t memid) { return 0; }
    uint32_t GetSize(uint32_t memid) { return sizeInBytes[memid]; }
    uint32_t GetThread() { return threadId; }
    uint32_t GetInstPtr(uint32_t memid) { return instPtr[memid]; }

    bool IsDP(uint32_t memseq) { return isDP[memseq]; }
    bool IsFP(uint32_t memseq) { return isFP[memseq]; }

    void SetIsDP(bool* newIsDP) { isDP = newIsDP; }
    void SetIsFP(bool* newIsFP) { isFP = newIsFP; }
    void SetSize(uint32_t* newSize) { sizeInBytes = newSize; }
    void SetInstPtr(uint64_t* newPtr) { instPtr = newPtr; }

    bool Verify();
};

class ArielFrontendHandler : public MemoryStreamHandler {
private:
    std::string shmemName;
    uint32_t traceRank; // Ariel picks a rank to collect data from
    SST::Core::Interprocess::SHMChild<SST::ArielComponent::ArielTunnel>*
      tunnelmgr;
    SST::ArielComponent::ArielTunnel* tunnel;
public:
    ArielFrontendHandler(std::string n);
    ~ArielFrontendHandler();

    void FinalizeTunnel();
    void InitializeTunnel();
    void InitializeTunnel(ArielFrontendHandler& h);
    void OutputStats(uint32_t threadID);
    void Print(std::ofstream& f);
    uint32_t Process(void* stats, uint64_t memSeq, bool ldstFlag,
      uint64_t* addresses, uint64_t length, bool memvecFlag);
    void ProcessInstructions(void* stats, uint64_t memSeq, uint64_t numInsns);
    void SetTraceRank(uint32_t r) { traceRank = r; }
    bool Verify() { return true; }
};


#endif /* _ArielFrontend_hpp_ */

