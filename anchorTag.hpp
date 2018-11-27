//
//  anchorTag.hpp
//  processBarcodeSample
//
//  Created by mark enstrom on 1/6/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//

#ifndef anchorTag_hpp
#define anchorTag_hpp
#include "KMer5.hpp"

class MergeLogEntry {
public:
    std::string   bc1;
    unsigned long count1;
    std::string   bc2;
    unsigned long count2;
    unsigned long hd;
    MergeLogEntry(std::string sb1, unsigned long c1,std::string sb2, unsigned long c2) {bc1=sb1;count1=c1;bc2=sb2;count2=c2;hd=0;}
};

class AnchorTags {
public:
    std::string _startIndex;
    std::string _startAnchor;
    std::string _endAnchor;
    std::string _fileName;
    bcodeMap    _exactMap;
    bcodeMap    _nMap;
    //
    // stats and errors from fastQ, these are sequences that
    // passed  index comparison but then failed start/end anchor or
    // #N in 20bp barcode
    int totalIndexReads;
    std::vector<std::string> _failAnchor  = std::vector<std::string>();
    std::vector<std::string> _failBCQual  = std::vector<std::string>();
    //
    // barcodes matched and error corrected against
    // master list
    //
    std::vector<IBarcode> _vecFound = std::vector<IBarcode>();
    std::vector<IBarcode> _vecNotFound = std::vector<IBarcode>();
    std::vector<IBarcode> _vecUNMATCH = std::vector<IBarcode>();
    std::vector<IBarcode> _vecDiscard = std::vector<IBarcode>();
    std::vector<IBarcode> _vecMerge = std::vector<IBarcode>();
    std::vector<MergeLogEntry> _mergeLog = std::vector<MergeLogEntry>();
	
    int findClose(IBarcode &ib, std::vector<IBarcode> sortList,int start,unsigned long maxHamming);
    
    void correctSmall();
    void compareToMaster(BarcodeMaster &master);
    void exactHamming(BarcodeMaster &master);
    void nHamming(BarcodeMaster &master);
    void handleUnmatched(BarcodeMaster &master);
    void combineDuplicates();
    void prepareMasterSource();
    void prepareSample(int finalHamming);
    void writeLog();
    void writeOutput();
    void writeOutputPre();
	void displayStats();
    
    
    AnchorTags(std::string i1,std::string a1, std::string a2, std::string fname) {
        _startIndex  = i1;
        _startAnchor = a1;
        _endAnchor   = a2;
        _fileName    = fname;
        _exactMap    = bcodeMap();
        _nMap        = bcodeMap();
        totalIndexReads = 0;
    }
};

//typedef std::pair<std::string,std::string> AnchorTags;


typedef std::vector<AnchorTags> AnchorTagVec;
#endif /* anchorTag_hpp */
