/* 
 * This file is part of the pebil project.
 * 
 * Copyright (c) 2010, University of California Regents
 * All rights reserved.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version
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
#include <Metasim.hpp>
#include <ArielFrontend.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cassert>

#include <sst/core/interprocess/shmchild.h>
#include "ariel_shmem.h"

using namespace SST::ArielComponent;
using namespace std;

void ArielFrontendTool::AddNewHandlers(AddressStreamStats* stats) {
    ArielFrontendHandler* oldHandler = (ArielFrontendHandler*)(handlers[0]);
    ArielFrontendHandler* newHandler = new ArielFrontendHandler(*oldHandler);
    stats->Handlers[indexInStats] = newHandler;
}

void ArielFrontendTool::AddNewStreamStats(AddressStreamStats* stats) {
    stats->Stats[indexInStats] = new ArielStats(stats->ThreadSeq);
    ((ArielStats*)(stats->Stats[indexInStats]))->SetIsDP(stats->IsDP);
    ((ArielStats*)(stats->Stats[indexInStats]))->SetIsFP(stats->IsFP);
    ((ArielStats*)(stats->Stats[indexInStats]))->SetSize(stats->SizeInBytes);
}

uint32_t ArielFrontendTool::CreateHandlers(uint32_t index, StringParser* parser) {
    indexInStats = index;
    char* e = parser->GetEnv("METASIM_SST_SHMEM");
    if (e == NULL) {
        ErrorExit("Please set METASIM_SST_SHMEM", MetasimError_Env);
    }

    ShmemName = (string)e;
    handlers.push_back(new ArielFrontendHandler(ShmemName));
    return handlers.size();
}

void ArielFrontendTool::FinalizeTool(DataManager<AddressStreamStats*>* AllData,
  SamplingMethod* Sampler) {
    
    AddressStreamStats* stats = AllData->GetData(AllData->GetFirstImage(),
      pthread_self());

}

ArielStats::ArielStats(uint32_t threadSeq) : threadId(threadSeq) {
}

ArielStats::~ArielStats(){
}

bool ArielStats::Verify(){
    return true;
}

//SST::Core::Interprocess::SHMChild<ArielTunnel> * tunnelmgr;
ArielFrontendHandler::ArielFrontendHandler(std::string n) : ShmemName(n) {
    tunnelmgr = new SST::Core::Interprocess::SHMChild<ArielTunnel>(ShmemName);
    tunnel = tunnelmgr->getTunnel();
}
ArielFrontendHandler::~ArielFrontendHandler() {
    if (tunnel != NULL) {
        ArielCommand ac;
        ac.command = ARIEL_PERFORM_EXIT;
        ac.instPtr = (uint64_t) 0;
        tunnel->writeMessage(0, ac);
        
        delete tunnel;
    }
    tunnel = NULL;
}

void ArielFrontendHandler::Print(ofstream& f){
    f << "ArielFrontendHandler" << ENDL;
}

uint32_t ArielFrontendHandler::Process(void* stats, uint64_t memSeq, 
  bool ldstFlag, uint64_t* addresses, uint64_t length, bool memvecFlag) {
    ArielStats* s = (ArielStats*)stats;
    ArielCommand ac;

    if (length <= 0)
        return 0;

    // Send Start instruction
    ac.command = ARIEL_START_INSTRUCTION;
    ac.instPtr = memSeq;
    ac.inst.instClass = ARIEL_INST_UNKNOWN;
    if (s->IsFP(memSeq)) {
        if (s->IsDP(memSeq))
            ac.inst.instClass = ARIEL_INST_DP_FP;
        else
            ac.inst.instClass = ARIEL_INST_SP_FP;
    }
    ac.inst.simdElemCount = length;
    tunnel->writeMessage(s->GetThread(), ac);

    for(int i = 0; i < length; i++) {
        uint64_t addr = addresses[i];
        if (addr != 0) {

            // if load
            if (ldstFlag)
                ac.command = ARIEL_PERFORM_READ;
            else
                ac.command = ARIEL_PERFORM_WRITE;

            ac.instPtr = memSeq;
            ac.inst.addr = addr;
            ac.inst.size = s->GetSize(memSeq);

            //if (ldstFlag)
            //fprintf(stderr, "ACC: ARIEL_PERFORM_READ: %d, %#lx\n", ac.instPtr,
            //  ac.inst.addr);
            //else
            //fprintf(stderr, "ACC: ARIEL_PERFORM_WRITE: %d, %#lx\n", ac.instPtr,
            //  ac.inst.addr);
            tunnel->writeMessage(s->GetThread(), ac);
        }
    }

    ac.command = ARIEL_END_INSTRUCTION;
    ac.instPtr = memSeq;
    tunnel->writeMessage(s->GetThread(), ac);
    return 0;

}
                
