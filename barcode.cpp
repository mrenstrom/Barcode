//
//  barcode.cpp
//  readBarcode
//
//  Created by mark enstrom on 3/18/17.
//  Copyright © 2017 Mark Enstrom. All rights reserved.
//
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <thread>
#include <ctime>

#include "barcode.hpp"
#include "KMer5.hpp"
#include "IBarcode.hpp"
#include "readLibrary.hpp"


using namespace std;

//int Barcode::q()
//{
//    int totalQ = 0;
//    for (int i = 0; i < strlen(qual); i++)
//    {
//        int qScore = (int)qual[i] - 33;
//        totalQ += qScore;
//    }
//    
//    return (int((float)totalQ/(float)(strlen(qual))));
//}
//
//int Barcode::markQuality()
//{
//    int Q = 3;
//    for (int i = 0; i < (strlen(qual)-1); i++)
//    {
//        int qScore = (int)qual[i] - 33;
//        if (qScore < 30)
//        {
//            //printf("bad Seq = %i   %c\n",qScore,seq[i]);
//            seq[i] = 'N';
//            Q--;
//        }
//    }
//    if (Q <= 0) return 0;
//    return 1;
//}



int hammingDistance(const char* p1,const char*  p2,int l)
{
    int dist = 0;
    for (int i = 0; i < l; i++)
    {
        if (p1[i] != p2[i]) dist += 1;
        if (dist >= 5) break;
    }
    return dist;
}


BarcodeMaster::BarcodeMaster() {
    _pbcExact   = new BCodeMap();
    _pqList     = new std::vector<std::string>();
    _pibList    = new std::vector<IBarcode>();
    _pkMap      = new KMER_OBJ();
}



/*--------------------------------------------------------------------------------------------
 * BarcodeMaster::InitFromFile
 *
 *  read in barcode P1 files which contain list of barcode,count of exact matched from original
 *  or exact matches with one or two Ns
 *
 *--------------------------------------------------------------------------------------------*/
bool BarcodeMaster::initFromFile(std::string file){
    std::string masterFile = file + std::string("MasterBarcodeList.txt");
    std::cout << "Init lib from file " << masterFile << "\n";
    _pibList = readLibrary(masterFile);
    if (_pibList == nullptr) return false;
    
    for (auto &ib : *_pibList) {
        _pbcExact->insert(std::pair<std::string,int>(ib.barcode(),ib._seqCount));
        _pqList->push_back(ib.barcode());
    }
    if (!_pkMap->load(file)) {
        std::cout << "build kmer map failed\n";
        return false;
    }
	std::cout << "loaded " << _pbcExact->size() << " seq for master" << "\n";
	std::cout << "loaded " << _pqList->size()   << " q   for master" << "\n";
	
    return true;
}

/*--------------------------------------------------------------------------------------------
 * BarcodeMaster::findNearest
 *   find the closest hamming dist for each barcode
 *
 *
 *--------------------------------------------------------------------------------------------*/
void BarcodeMaster::cleanIBarcodes()
{
    auto &ibList = *_pibList;
    auto _pibNew  = new std::vector<IBarcode>();
    for (int i = 0; i < ibList.size(); i++) {
        auto &ib1 = ibList[i];
        if (ib1._valid) {
            _pibNew->push_back(ib1);
        }
    }
    delete _pibList;
    _pibList = _pibNew;
}

void BarcodeMaster::findNearest(std::vector<IBarcode> &ibVec)
{
    cleanIBarcodes();
//    auto &ibList = *_pibList;
//    //
//    // use mutli-threads to compute hammings
//    //
//    std::thread t1 {F(ibList,ibVec,4,0)};
//    std::thread t2 {F(ibList,ibVec,4,1)};
//    std::thread t3 {F(ibList,ibVec,4,2)};
//    std::thread t4 {F(ibList,ibVec,4,3)};
//    
//    t1.join();
//    t2.join();
//    t3.join();
//    t4.join();
    //
    // assign results
    //
}


