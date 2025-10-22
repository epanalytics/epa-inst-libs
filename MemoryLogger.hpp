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

#ifndef _MemoryLogger_hpp_
#define _MemoryLogger_hpp_

#include <AddressStreamBase.hpp>
#include <string>


class MemoryLoggerTool : public AddressStreamTool {
  protected:
    std::ofstream logFileStream;
  public:
    MemoryLoggerTool() : AddressStreamTool() {}
    virtual void AddNewHandlers(AddressStreamStats* stats);
    virtual void AddNewStreamStats(AddressStreamStats* stats);
    virtual uint32_t CreateHandlers(uint32_t index, StringParser* parser);
    virtual void FinalizeTool(DataManager<AddressStreamStats*>* AllData,
      SamplingMethod* Sampler);
    virtual void MemoryLogName(AddressStreamStats* stats, std::string& oFile);
    void WriteLog(std::string s);
};

class LoggerStats : public StreamStats {
private:
    uint32_t* sizeInBytes;

public:

    LoggerStats() {}
    virtual ~LoggerStats() {}

    uint64_t GetAccessCount(uint32_t memid) { return 0; }
    uint32_t GetSize(uint32_t memid) { return sizeInBytes[memid]; }
    void SetSize(uint32_t* newSize) { sizeInBytes = newSize; }

    bool Verify() { return 0; }
};

class MemoryLoggerHandler : public MemoryStreamHandler {
protected:
    MemoryLoggerTool* memLoggerTool;
public:
    MemoryLoggerHandler(MemoryLoggerTool* tool);
    MemoryLoggerHandler(MemoryLoggerHandler& h);
    ~MemoryLoggerHandler();

    void Print(std::ofstream& f);
    uint32_t Process(void* stats, uint64_t memSeq, bool ldstFlag,
      uint64_t* addresses, uint64_t length, bool memvecFlag);
    bool Verify() { return true; }
};


#endif /* _MemoryLogger_hpp_ */

