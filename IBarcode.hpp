//
//  IBarcode.hpp
//  barcode2
//
//  Created by mark enstrom on 12/6/17.
//  Copyright © 2017 mark enstrom. All rights reserved.
//

#ifndef IBarcode_hpp
#define IBarcode_hpp

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <list>
#include <stdlib.h>
#include <assert.h>
#include <ctime>
#include <unistd.h>



class IBarcode {
public:
    unsigned long long  _u0;
    unsigned long long  _u1;
    int                 _N;
    std::string         _bc;
    size_t              _seqCount;
    size_t              _seqCountExact;
    size_t              _seqCountN1;
    size_t              _seqCountMerge;
    bool                _valid;
    int                 _index;
    unsigned long       _minHamming;
    
    IBarcode(std::string seq,size_t num);
    std::string barcode();
    unsigned long hammingDist(IBarcode &seq2);
    unsigned long hammingDistSlow(IBarcode &seq2);
    int countN();
    unsigned long hammingDistN(IBarcode &seq2,int numN);
};

typedef std::vector<IBarcode> IBarVector;


//struct F {
//    int _offset;
//    int _stride;
//    int _hFreq[32] = {};
//    int _minList[32] = {};
//    std::vector<IBarcode> &_ibVec1;
//    std::vector<IBarcode> &_ibVec2;
//    
//    F(std::vector<IBarcode> &v1,std::vector<IBarcode> &v2,int s,int o):_stride{s},_offset{o},_ibVec1{v1},_ibVec2{v2}
//    {
//    }
//    void operator()(){barcodeDist();}
//    void barcodeDist();
//};

#endif /* IBarcode_hpp */

