//
//  IBarcode.cpp
//  barcode2
//
//  Created by mark enstrom on 12/6/17.
//  Copyright © 2017 mark enstrom. All rights reserved.
//

#include "IBarcode.hpp"
#include <thread>
#include <functional>

typedef std::chrono::milliseconds TimeT;

unsigned long IBarcode::hammingDist(IBarcode &seq2)
{
    unsigned long d0 = __builtin_popcountll(_u0 ^ seq2._u0) >> 1;
    unsigned long d1 = __builtin_popcountll(_u1 ^ seq2._u1) >> 1;
    return d0+d1;
}
//
// N = 0xFF, N xor (A,C,G,T) = 3 bits set instead of 2
//
unsigned long IBarcode::hammingDistN(IBarcode &seq2,int numN)
{
    unsigned long d0 = __builtin_popcountll(_u0 ^ seq2._u0);
    unsigned long d1 = __builtin_popcountll(_u1 ^ seq2._u1);
    return (d0+d1-numN)/2;
}


unsigned long IBarcode::hammingDistSlow(IBarcode &seq2)
{
    std::string s2 = seq2.barcode();
    unsigned long ret = 0;
    for (int i = 0; i < (int)_bc.length(); i++) {
        auto c1 = _bc[i];
        auto c2 = s2[i];
        if (c1 != c2) ret++;
    }
    return ret;
}

//                 A b C d e f G h i j k l m n o p q r s T u v w x y z
int charMap[26] = {1,0,2,0,0,0,4,0,0,0,0,0,0,15,0,0,0,0,0,8,0,0,0,0,0,0};

IBarcode::IBarcode(std::string seq1,size_t num)
{
    _u0 = 0;
    _u1 = 0;
    _bc = seq1;
    _seqCount = num;
    _seqCountExact = num;
    _seqCountN1 = 0;
    _valid = true;
    _seqCountMerge  = 0;
    _index = 0;
    _minHamming = 20;
    
    const char *seq = seq1.c_str();
    for (int i = 0; i < 16 ; i++){
        unsigned long long c = charMap[seq[i] - 'A'];
        _u0 |= (c << (4 * i));
    }
    for (int i = 16; i < 20 ; i++){
        unsigned long long c = charMap[seq[i] - 'A'];
        _u1 |= (c << (4 * (i-16)));
    }
    _N = countN();
}

//
// inverse translation
//
std::string IBarcode::barcode()
{
    return (_bc);
}
/*--------------------------------------------------------------------------------------------
 * countN
 *
 *
 *
 *--------------------------------------------------------------------------------------------*/
int IBarcode::countN() {
    int ret = 0;
    for (auto c : _bc) {
        if (c == 'N') ret++;
    }
    return ret;
}
/*--------------------------------------------------------------------------------------------
 * F::barcodeDist
 *
 *  multi-thread "functor" for computing hamming
 *
 *--------------------------------------------------------------------------------------------*/
//void F::barcodeDist() {
//    std::vector<IBarcode> &v1 = _ibVec1;
//    std::vector<IBarcode> &v2 = _ibVec2;
//    
//    int offset = _offset;
//    int loopCount = 0;
//    //
//    // report every 120,000,000 barcodes
//    //
//    int loopDisplay  = 1000;
//    
//    auto start = std::chrono::steady_clock::now();
//    unsigned long long  totalCmp = 0;
//    //
//    // use integers to avoid double dist comp  ie 1 <-> 2  == 2 <-> 1
//    //
//    for (int i = offset; i < v2.size() ; i+=_stride) {
//        auto &seq2 = v2[i];
//        int minH = 20;
//        int index = 0;
//        //
//        // how many "N" in query seq
//        //
//        int numN = seq2._N;
//        for (int k = 0; k < v1.size(); k++) {
//            auto &seq1 = v1[k];
//            unsigned long n;
//            if (numN == 0){
//                n = seq1.hammingDist(seq2);
//            } else {
//                n = seq1.hammingDistN(seq2,numN);
//            }
//            //
//            // save min hamming and index of min seq
//            //
//            if ((n != 0) && (n < minH)) {
//                minH= (int)n;
//                index = k;
//            }
//            totalCmp += 1;
//        }
//        seq2._minHamming = minH;
//        seq2._index = index;
//        //
//        // report
//        //
//        if ((++loopCount%loopDisplay)==0){
//            std::cout << "------" << loopCount << "-----" << v2.size()/_stride << "-----------------------\n";
//            auto duration = std::chrono::duration_cast<TimeT>(std::chrono::steady_clock::now() - start);
//            float bcPerSec = (float)totalCmp/(float)(duration.count());
//            std::cout << "processing " << bcPerSec << " bc/mS " << std::endl;
//            totalCmp = 0;
//            start = std::chrono::steady_clock::now();
//        }
//        //if (loopCount > 10000) break;
//    }
//    return;
//}
//
//

