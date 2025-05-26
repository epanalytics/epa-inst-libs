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
#include <MemoryLogger.hpp>

#include <iostream>
#include <cstring>
#include <cassert>
#include <sstream>

using namespace std;

void MemoryLoggerTool::AddNewHandlers(AddressStreamStats* stats) {
    string logFileName;
    MemoryLogName(stats, logFileName);
    TryOpen(logFileStream, logFileName.c_str());
    MemoryLoggerHandler* oldHandler = (MemoryLoggerHandler*)(handlers[0]);
    MemoryLoggerHandler* newHandler = new MemoryLoggerHandler(*oldHandler);
    stats->Handlers[indexInStats] = newHandler;
}

void MemoryLoggerTool::AddNewStreamStats(AddressStreamStats* stats) {
    stats->Stats[indexInStats] = new LoggerStats();
    ((LoggerStats*)(stats->Stats[indexInStats]))->SetSize(stats->SizeInBytes);
}

uint32_t MemoryLoggerTool::CreateHandlers(uint32_t index, StringParser* parser) {
    indexInStats = index;
    handlers.push_back(new MemoryLoggerHandler(this));
    return handlers.size();
}

void MemoryLoggerTool::FinalizeTool(DataManager<AddressStreamStats*>* AllData,
  SamplingMethod* Sampler) {
    logFileStream.close();
    return;
}

void MemoryLoggerTool::MemoryLogName(AddressStreamStats* stats,
  std::string& oFile) {
    oFile.clear();
    oFile.append(stats->Application);
    oFile.append(".r");
    AppendRankString(oFile);
    oFile.append(".t");
    AppendTasksString(oFile);
    oFile.append(".log");
}

void MemoryLoggerTool::WriteLog(string s) {
    logFileStream << s;
}

MemoryLoggerHandler::MemoryLoggerHandler(MemoryLoggerTool* tool) :
  memLoggerTool(tool) {
}

MemoryLoggerHandler::MemoryLoggerHandler(MemoryLoggerHandler& h) {
    memLoggerTool = h.memLoggerTool;
}

MemoryLoggerHandler::~MemoryLoggerHandler() {
}

void MemoryLoggerHandler::Print(ofstream& f){
    f << "MemoryLoggerHandler" << ENDL;
}

uint32_t MemoryLoggerHandler::Process(void* stats, uint64_t memSeq, 
  bool ldstFlag, uint64_t* addresses, uint64_t length, bool memvecFlag) {
    LoggerStats* s = (LoggerStats*)stats;

    if (length <= 0)
        return 0;

    for(int i = 0; i < length; i++) {
        uint64_t addr = addresses[i];
        stringstream logStream;
        logStream << "type=";
        // if load
        if (ldstFlag)
            logStream << "R";
        else
            logStream << "W";

        logStream << "\t" << "memSeq=" << std::dec << memSeq;
        logStream << "\t" << "size=" << std::dec << s->GetSize(memSeq);
        logStream << "\t" << "addr=0x" << std::hex << addr;
        logStream << std::endl;
        memLoggerTool->WriteLog(logStream.str());
    }

    return 0;

}
