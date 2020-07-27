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

#include <InstrumentationCommon.hpp>
#include <DataManager.hpp>
#include <DynamicInstrumentation.hpp>
#include <Metasim.hpp>
#include <ThreadedCommon.hpp>
//#include <ReuseDistance.hpp>  // external SpatialLocality
#include <SpatialLocality.hpp>
#include <SpatialLocalityPerMemOp.hpp>
#include <AddressStreamStats.hpp>

#include <cstring>

using namespace std;

void SpatialLocalityPerMemOpTool::AddNewHandlers(AddressStreamStats* stats) {
    SpatialLocalityPerMemOpHandler* oldHandler = (SpatialLocalityPerMemOpHandler*)(handlers[0]);
    SpatialLocalityPerMemOpHandler* newHandler = new SpatialLocalityPerMemOpHandler(
      *oldHandler);
    stats->Handlers[indexInStats] = newHandler;
}

void SpatialLocalityPerMemOpTool::AddNewStreamStats(AddressStreamStats* stats) {
    stats->Stats[indexInStats] = new SpatialStreamStats(stats);
}

uint32_t SpatialLocalityPerMemOpTool::CreateHandlers(uint32_t index, StringParser* parser) {
    indexInStats = index;

    uint32_t spatialWindow;
    uint32_t spatialBin;
    uint32_t spatialNMAX;
    if (!(parser->ReadEnvUint32("METASIM_SPATIAL_WINDOW", &spatialWindow))) {
        spatialWindow = 1;
    }
    if (!(parser->ReadEnvUint32("METASIM_SPATIAL_BIN", &spatialBin))) {
        spatialBin = 1;
    }
    if (!(parser->ReadEnvUint32("METASIM_SPATIAL_NMAX", &spatialNMAX))) {
        spatialNMAX = ReuseDistance::Infinity;
    }

    handlers.push_back(new SpatialLocalityPerMemOpHandler(spatialWindow, spatialBin,
      spatialNMAX));

    return handlers.size();
}

void SpatialLocalityPerMemOpTool::FinalizeTool(DataManager<AddressStreamStats*>* 
  AllData, SamplingMethod* sampler) {
    string oFile;
    const char* fileName;

    AddressStreamStats* stats = AllData->GetData(pthread_self());

    ofstream SpatialLocFile;
    SpatialLocalityPerMemOpFileName(stats, oFile);
    fileName = oFile.c_str();
    inform << "Printing spatial locality results to " << fileName << ENDL;
    TryOpen(SpatialLocFile, fileName);

    for (set<image_key_t>::iterator iit = AllData->allimages.begin(); iit != AllData->allimages.end(); iit++){
        for(DataManager<AddressStreamStats*>::iterator it = AllData->begin(*iit); it != AllData->end(*iit); ++it){

            thread_key_t thread = it->first;
            AddressStreamStats* s = it->second;

            SpatialLocFile << "IMAGE" << TAB << hex << (*iit) 
              << TAB << "THREAD" << TAB << dec 
              << AllData->GetThreadSequence(thread) << ENDL;

            SpatialLocalityPerMemOpHandler* sd = (SpatialLocalityPerMemOpHandler*)(s->Handlers[
              indexInStats]);
            assert(sd);

            inform << "Spatial locality bins for " 
              << hex << s->Application << " Thread " 
              << AllData->GetThreadSequence(thread) << ENDL;

            pebil_map_type<uint64_t, vector<uint64_t>*> blockMapper;
            vector<uint64_t> keys;
            sd->PrintHeader(SpatialLocFile);

            for (pebil_map_type<uint64_t, uint64_t>::iterator mem2blockIT =
              sd->memopBlockMapper->begin(); mem2blockIT != sd->memopBlockMapper->end();
              mem2blockIT++) {
                uint64_t memOp = mem2blockIT->first;
                uint64_t block = mem2blockIT->second;
                if (blockMapper.find(block) == blockMapper.end()){
                    vector<uint64_t>* newVec = new vector<uint64_t>();
                    newVec->push_back(memOp);
                    std::pair<uint64_t, vector<uint64_t>*> kvp (block, newVec);
                    blockMapper.insert(kvp);
                    keys.push_back(block);
                } else {
                    vector<uint64_t>* oldVec = blockMapper.find(block)->second;
                    oldVec->push_back(memOp);
                }
            }

            sort(keys.begin(), keys.end());
            //vector<vector<uint64_t>*> inOrder;
            for(vector<uint64_t>::iterator blockIT = keys.begin(); blockIT != keys.end();
              blockIT++){
                uint64_t block = *blockIT;
                vector<uint64_t>* blockVec = blockMapper.find(block)->second;
                sd->PrintBlockInfo(SpatialLocFile, block, blockVec);
                //inOrder.push_back(blockVec);
                
            }

            /*SpatialLocFile << "-------------------------------------" << ENDL;

            for (pebil_map_type<uint64_t, ReuseDistance*>::iterator itbit =
              sd->mapInternalHandler->begin(); itbit != sd->mapInternalHandler->end();
              itbit++){

                itbit->second->Print(SpatialLocFile);
            }*/
            //sd->Print(SpatialLocFile);
        }
    }
    SpatialLocFile.close();

}

void SpatialLocalityPerMemOpHandler::PrintHeader(ostream& f){
    //line 1
    f << "# SPATIALSTATS <block_count> <block_hash> <window_size> <bin_indiv> "
      << "<max_track> <id_count> <tot_access> <tot_miss>" << ENDL;

    //line2
    f << "#" << TAB << "SPATIALID <id> <id_access> <id_miss>" << ENDL;

    //line3
    f << "#" << TAB << TAB << "<bin_lower_bound> <bin_upper_bound> <bin_count>" << ENDL;
}

void SpatialLocalityPerMemOpTool::SpatialLocalityPerMemOpFileName(AddressStreamStats* stats, 
  string& oFile){
    oFile.clear();
    oFile.append(stats->Application);
    oFile.append(".r");
    AppendRankString(oFile);
    oFile.append(".t");
    AppendTasksString(oFile);
    oFile.append(".spatialPerMemOp");
}

SpatialLocalityPerMemOpHandler::SpatialLocalityPerMemOpHandler(uint64_t w, uint64_t b, uint64_t
  n) : ReuseDistanceHandler(0, 0) {
    delete internalHandler;
    internalHandler = nullptr;
    //internalHandler = new SpatialLocality(w, b, n);
    mapInternalHandler = new pebil_map_type<uint64_t, ReuseDistance*>();
    //blockMemopMapper = new pebil_map_type<uint64_t, set<uint64_t>*>();
    memopBlockMapper = new pebil_map_type<uint64_t, uint64_t>();
    window = w;
    bin = b;
    nmax = n;
}

SpatialLocalityPerMemOpHandler::SpatialLocalityPerMemOpHandler(SpatialLocalityPerMemOpHandler& h) :
  ReuseDistanceHandler(0, 0) {
    delete internalHandler;
    internalHandler = nullptr;
    //internalHandler = new SpatialLocality((SpatialLocality*)h.internalHandler);
    mapInternalHandler = new pebil_map_type<uint64_t, ReuseDistance*>();
    //blockMemopMapper = new pebil_map_type<uint64_t, set<uint64_t>*>();
    memopBlockMapper = new pebil_map_type<uint64_t, uint64_t>();
    window = h.window;
    bin = h.bin;
    nmax = h.nmax;
}

//TODO
SpatialLocalityPerMemOpHandler::~SpatialLocalityPerMemOpHandler(){

}

uint32_t SpatialLocalityPerMemOpHandler::Process(void* stats, BufferEntry* access){

    ReuseStreamStats* s = (ReuseStreamStats*)stats;
    uint64_t memOp = access->memseq;
    ReuseDistance* individualMemOpSpatialLocality;

    //not in map, add it 
    if (mapInternalHandler->find(memOp) == mapInternalHandler->end()) {
        individualMemOpSpatialLocality = new SpatialLocality(window,bin,nmax);
        std::pair<uint64_t, ReuseDistance*> kvp (memOp, individualMemOpSpatialLocality);
        mapInternalHandler->insert(kvp);
        uint64_t block = s->GetBlock(memOp);
        std::pair<uint64_t, uint64_t> kvp2 (memOp, block);
        memopBlockMapper->insert(kvp2);
    } else { //in map, grap it
        individualMemOpSpatialLocality = mapInternalHandler->find(memOp)->second;
    }

    /*uint64_t block = s->GetBlock(memOp);

    //keep track of blockid -> memop 
    if (blockMemopMapper->find(block) == blockMemopMapper->end()) {
        set<uint64_t>* newSet = new set<uint64_t>();
        newSet->insert(memOp);
        std::pair<uint64_t, set<uint64_t>*> kvp (block, newSet);
        blockMemopMapper->insert(kvp);
    } else {
        set<uint64_t>* oldSet = blockMemopMapper->find(block)->second;
        oldSet->insert(memOp);
    }*/

    //process it
    ReuseEntry entry = ReuseEntry();
    entry.id = s->GetHash(access->memseq);
    entry.address = access->address;
    individualMemOpSpatialLocality->Process(entry);
}

void SpatialLocalityPerMemOpHandler::SkipAddresses(uint32_t numToSkip){
    //map look up? TODO
    for( pebil_map_type<uint64_t, ReuseDistance*>::iterator it = mapInternalHandler->begin();
      it != mapInternalHandler->end(); it++){
        ReuseDistance* current = it->second;
        current->SkipAddresses(numToSkip);
    }
}

void SpatialLocalityPerMemOpHandler::PrintBlockInfo(std::ostream& f, uint64_t block, 
  std::vector<uint64_t>* vec){

    reuse_map_type<uint64_t,uint64_t> BinTotal;
     
    uint64_t blockHash;

    uint64_t tot = 0, mis = 0, count = 0;
    for (vector<uint64_t>::const_iterator it = vec->begin(); it != vec->end(); it++){
        uint64_t memOp = (*it);
        ReuseDistance* current = mapInternalHandler->find(memOp)->second;
        vector<uint64_t> indices;
        current->GetIndices(indices);
        ReuseStats* r = current->GetStats(indices[0]);
        if (r != nullptr){
            tot += r->GetAccessCount();
            mis += r->GetMissCount();
            blockHash = indices[0]; 
        }
    }

    f << "SPATIALSTATS"
      << TAB << block << dec
      << TAB << hex << blockHash << dec
        //capacity is size
      << TAB << dec << window//capacity
        //binindividual is bin
      << TAB << bin //binindividual
        //maxtracking is nmax
      << TAB << nmax //maxtracking
      << TAB << vec->size()
      << TAB << tot
      << TAB << mis
      << ENDL;

    for (vector<uint64_t>::const_iterator it = vec->begin(); it != vec->end(); it++){
        uint64_t memOp = (*it);
        ReuseDistance* current = mapInternalHandler->find(memOp)->second;
        PrintMemOpInfo(f, memOp, current, BinTotal);
    }
}

void SpatialLocalityPerMemOpHandler::PrintMemOpInfo( ostream& f, 
  uint64_t memop, ReuseDistance* rd, reuse_map_type<uint64_t,uint64_t> BinTotal){

    vector<uint64_t> indices;
    rd->GetIndices(indices);
    ReuseStats* r = rd->GetStats(indices[0]);
    
    if (r != nullptr) {
        f << TAB << "SPATIALID"
          << TAB << memop << dec
          << TAB << r->GetAccessCount()
          << TAB << r->GetMissCount()
          << ENDL;

        r->Print(f, BinTotal);
    }
}
