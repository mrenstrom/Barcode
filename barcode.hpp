//
//  barcode.hpp
//  readBarcode
//
//  Created by mark enstrom on 3/18/17.
//  Copyright © 2017 Mark Enstrom. All rights reserved.
//

#ifndef barcode_hpp
#define barcode_hpp
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
#include <algorithm>
//#include "IBarcode.hpp"
//#include "KMer5.hpp"
//#include "readLibrary.hpp"
//#include "readFile.hpp"
//
// kMer MAP using x base mers
//   every mer has a map of seq_pos : vector<string> for all of the barcodes
//   that match current mer at each sequence position
//
// for example mer 'AAAAA' will have a map of 5 starting positions 0,5,10,15
// so position 0 will have a vector of all barcodes that have 'AAAAA' at position 0
//
//
typedef std::unordered_map<std::string,size_t> bcodeMap;


//typedef std::vector<std::string> TagVector;
//
// std data structure definitions
//


//class Barcode {
//public:
//    char hdr[100];
//    char seq[100];
//    char qual[100];
//    char tag[100];
//    bool achorFound;
//public:
//    Barcode():achorFound{false}{};
//    void shorten(std::vector<std::string> tags);
//    int q();
//    int markQuality();
//};
//
//
//
//
//typedef std::unordered_map<std::string,std::vector<Barcode>> BarcodeMap;
//
//typedef std::list<Barcode> BarcodeList;


//void compareBarcodes2(bcodeMap& bcExact, bcodeMap& bcN);



//
// kMer MAP using x base mers
//   every mer has a map of seq_pos : vector<string> for all of the barcodes
//   that match current mer at each sequence position
//
// for example mer 'AAAAA' will have a map of 5 starting positions 0,5,10,15
// so position 0 will have a vector of all barcodes that have 'AAAAA' at position 0
//
//
class MetaFileData;

typedef std::unordered_map<std::string,size_t> BCodeMap;

class SeqMatch;
class KMER_OBJ;
class IBarcode;

int hammingDistance(const char* p1,const char*  p2,int l);

class BarcodeMaster {
public:
    int                         _dbg = 0;
    BCodeMap*                   _pbcExact;   // unordered map of string:size for all exact barcodes
    std::vector<std::string>*   _pqList;     // string barcodes in list form
    std::vector<IBarcode>*      _pibList;    // IBarcode obects to store more info and allow quick hamming
    KMER_OBJ* _pkMap;
    BarcodeMaster();
    bool initFromFile(const std::string fileName);
    void findNearest(std::vector<IBarcode> &ibVec);
private:
    void cleanIBarcodes();
};



#endif /* barcode_hpp */
