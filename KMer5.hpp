
//
//  KMer5.hpp
//  barcodeP2
//
//  Created by mark enstrom on 12/18/17.
//  Copyright © 2017 mark enstrom. All rights reserved.
//

#ifndef KMer5_hpp
#define KMer5_hpp

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
#include "IBarcode.hpp"


typedef std::map<int,std::vector<int>> KMER_POS;
typedef std::map<std::string,KMER_POS> KMER_MAP;

class KMer5 {
public:
    int     index;
    int     position;
};

class SeqMatch {
public:
    std::string seq1;
    std::string seq2;
    int index1;
    int index2;
    int dist;
};


class KMER_OBJ {
public:
    int dbg4 = 0;
    KMER_MAP *_pMap;
	std::vector<std::string> _qList = std::vector<std::string> ();
    KMER_OBJ();
    ~KMER_OBJ();
    void buildMap(std::vector<std::string>qList);
    void save(std::string filename);
    bool load(std::string base);
	std::vector<int> merCompV(std::string seq,int maxDist);
    int  merComp(std::string seq,int maxDist);
	int  merCompMax(std::string seq,int maxDist,const std::vector<IBarcode> &bar);
	
};

#endif /* KMer5_hpp */


