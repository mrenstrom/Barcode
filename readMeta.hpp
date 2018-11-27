//
//  readMeta.hpp
//  barcode2
//
//  Created by mark enstrom on 12/4/17.
//  Copyright © 2017 mark enstrom. All rights reserved.
//

#ifndef readMeta_hpp
#define readMeta_hpp



class MetaFileData {
public:
    std::string _masterFile;
    std::string _inputFileName;
    std::vector<std::string> _barcode = std::vector<std::string>();
    AnchorTagVec _tagVec;
	bool _buildMaster;
    bool _sequenceError;
    bool _useMaster;
    bool _samplePrep;
    bool _hammingDiag;
    //
    // error results: quality before index/anchor can be established
    // no valid index or anchor for runs with no index
    //
    std::vector<std::string> _failQual   = std::vector<std::string>();
    std::vector<std::string> _failAnchor = std::vector<std::string>();
    int totalReads = 0;
public:
    bool open(std::string metaFileName);
    MetaFileData() {
	_buildMaster = false; 
	_useMaster = false; 
	_sequenceError = false;
	_samplePrep = false;
        _hammingDiag = false;
    };
};

#endif /* readMeta_hpp */

