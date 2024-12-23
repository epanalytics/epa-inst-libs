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
#include <cstring>
#include <cassert>

#include <sst/core/interprocess/shmchild.h>
#include "ariel_shmem.h"

using namespace SST::ArielComponent;
using namespace std;

void ArielFrontendTool::AddNewHandlers(AddressStreamStats* stats) {
    ArielFrontendHandler* oldHandler = (ArielFrontendHandler*)(handlers[0]);
    ArielFrontendHandler* newHandler = new ArielFrontendHandler(*oldHandler);
    tunnelUsers.push_back(newHandler);
    stats->Handlers[indexInStats] = newHandler;
}

void ArielFrontendTool::AddNewStreamStats(AddressStreamStats* stats) {
    stats->Stats[indexInStats] = new ArielStats(stats->ThreadSeq);
    ((ArielStats*)(stats->Stats[indexInStats]))->SetIsDP(stats->IsDP);
    ((ArielStats*)(stats->Stats[indexInStats]))->SetIsFP(stats->IsFP);
    ((ArielStats*)(stats->Stats[indexInStats]))->SetSize(stats->SizeInBytes);
    ((ArielStats*)(stats->Stats[indexInStats]))->SetInstPtr(stats->Addresses);
}

uint32_t ArielFrontendTool::CreateHandlers(uint32_t index, StringParser* parser) {
    indexInStats = index;
    char* e = parser->GetEnv("METASIM_SST_SHMEM");
    if (e == NULL) {
        ErrorExit("Please set METASIM_SST_SHMEM", MetasimError_Env);
    }
    shmemName = (string)e;

    uint32_t mpiUsage;
    if (parser->ReadEnvUint32("METASIM_SST_USE_MPI", &mpiUsage))
        usesMPI = (mpiUsage != 0);
    else
        usesMPI = false;

    uint32_t rankToTrace;
    if (parser->ReadEnvUint32("METASIM_SST_TRACE_RANK", &rankToTrace))
        traceRank = rankToTrace;
    else
        traceRank = 0;

    tunnelCreator = new ArielFrontendHandler(shmemName);
    tunnelCreator->SetTraceRank(traceRank);
    // If this is MPI, initialize the tunnel *after* MPI_Init gets called
    // (see NotifyDoneMPIInit()). If this is not MPI, initialize the tunnel
    // now. The user must set METASIM_SST_USE_MPI to have it work with MPI.
    // This is okay, since Ariel requires a special mpi-launcher for MPI
    // apps so we will know ahead of time if we are using MPI
    if (!usesMPI)
        tunnelCreator->InitializeTunnel();
    handlers.push_back(tunnelCreator);
    return handlers.size();
}

void ArielFrontendTool::FinalizeTool(DataManager<AddressStreamStats*>* AllData,
  SamplingMethod* Sampler) {

    // Currently only passing data from rank 0 // TODO --> Generalize
    if (GetTaskId() == traceRank)
        tunnelCreator->FinalizeTunnel();
}

// Use this to Initialize the tunnel for MPI applications.
// There can only be one initialized tunnel, so we give it to the rank that
// is being traced for Ariel
void ArielFrontendTool::NotifyDoneMPIInit() {
    if (!usesMPI) {
        warn << "METASIM_SST_USE_MPI was not set but MPI_Init was found"
          << ENDL;
        return;
    }

    if (GetTaskId() != traceRank)
        return;

    // There can only be one initialized tunnel, so threads need to share it
    // Have our "tunnel creator" initialize the tunnel and then other threads
    // "use" the tunnel by copying over its pointer
    tunnelCreator->InitializeTunnel();
    for (auto uItr = tunnelUsers.begin(); uItr != tunnelUsers.end(); uItr++)
        (*uItr)->InitializeTunnel(*tunnelCreator);
}

ArielStats::ArielStats(uint32_t threadSeq) : threadId(threadSeq) {
}

ArielStats::~ArielStats(){
}

bool ArielStats::Verify(){
    return true;
}

// Type info: SST::Core::Interprocess::SHMChild<ArielTunnel> * tunnelmgr;
ArielFrontendHandler::ArielFrontendHandler(std::string n) : shmemName(n),
  traceRank(0), tunnelmgr(NULL), tunnel(NULL) {
}

ArielFrontendHandler::~ArielFrontendHandler() {
    // tunnel deleted during FinalizeTool/FinalizeTunnel
}

void ArielFrontendHandler::FinalizeTunnel() {
    ArielCommand ac;
    ac.command = ARIEL_PERFORM_EXIT;
    ac.instPtr = (uint64_t) 0;
    tunnel->writeMessage(0, ac);
    delete tunnelmgr;
}

// Create an instance of an Ariel tunnel (can only be done ONCE -- one thread,
// one rank, one image (TODO))
void ArielFrontendHandler::InitializeTunnel() {
    if (tunnel != NULL)
        return;
    tunnelmgr = new SST::Core::Interprocess::SHMChild<ArielTunnel>(shmemName);
    tunnel = tunnelmgr->getTunnel();
}

// Give a handler access to an already-created Ariel tunnel
void ArielFrontendHandler::InitializeTunnel(ArielFrontendHandler& h) {
    if (tunnel != NULL)
        return;
    if (h.tunnel == NULL) {
        ErrorExit("Attempting to initialize a tunnel user but the creator has "
          "not created a tunnel", MetasimError_None);
    }

    tunnel = h.tunnel;
}

void ArielFrontendHandler::Print(ofstream& f){
    f << "ArielFrontendHandler" << ENDL;
}

uint32_t ArielFrontendHandler::Process(void* stats, uint64_t memSeq, 
  bool ldstFlag, uint64_t* addresses, uint64_t length, bool memvecFlag) {
    ArielStats* s = (ArielStats*)stats;
    ArielCommand ac;

    static std::set<uint64_t> reportedMemSeqs;

    if (length <= 0)
        return 0;

    if (GetTaskId() != traceRank)
        return 0;

    if (tunnel == NULL) {
        fprintf(stderr, "ERROR: Rank %d is attempting to process the buffer "
          "but no tunnel has been initialized\n", GetTaskId());
    }

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

    auto myInstClass = ac.inst.instClass;
    ac.inst.simdElemCount = length;
    tunnel->writeMessage(s->GetThread(), ac);

    //if (length > 1) {
    //    if (ldstFlag)
    //        fprintf(stderr, "ACC: Found a gather: %d with element size %d\n", memSeq, s->GetSize(memSeq));
    //    else
    //        fprintf(stderr, "ACC: Found a scatter: %d with element size %d\n", memSeq, s->GetSize(memSeq));
    //}

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

           // if (ldstFlag)
           //     fprintf(stderr, "ARIEL_PERFORM_READ: %#lx, %d, %#lx\n",
           //       ac.instPtr, ac.inst.size, ac.inst.addr);
           // else
           //     fprintf(stderr, "ARIEL_PERFORM_WRITE: %#lx, %d, %#lx\n",
           //       ac.instPtr, ac.inst.size, ac.inst.addr);
            tunnel->writeMessage(s->GetThread(), ac);
        }
    }
    if (reportedMemSeqs.count(memSeq) == 0) {
        fprintf(stderr, "ACC_ARIEL: instPtr=%#lx instClass=%d simdElemCount=%d command=%d size=%d length=%d\n", s->GetInstPtr(memSeq), myInstClass, length, ac.command, ac.inst.size, length);
        reportedMemSeqs.insert(memSeq);

    }

    ac.command = ARIEL_END_INSTRUCTION;
    ac.instPtr = memSeq;
    tunnel->writeMessage(s->GetThread(), ac);
    return 0;

}

void ArielFrontendHandler::ProcessInstructions(void* stats, uint64_t memSeq,
  uint64_t numInsns) {
    ArielStats* s = (ArielStats*)stats;
    ArielCommand ac;

    if (GetTaskId() != traceRank)
        return;

    // Send NOOP instruction for each non memory instruction
    ac.command = ARIEL_NOOP;
    ac.instPtr = memSeq;
    for (auto i = 0; i < numInsns; i++)
        tunnel->writeMessage(s->GetThread(), ac);
}
