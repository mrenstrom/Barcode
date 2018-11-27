//
//  readMeta.cpp
//  barcode2
//
//  Created by mark enstrom on 12/4/17.
//  Copyright © 2017 mark enstrom. All rights reserved.
//

#include "barcode.hpp"
#include "IBarcode.hpp"
#include "anchorTag.hpp"
#include "readMeta.hpp"
//
// metafile very brittle:requires:
//
// master/noMaster
// path_to_process_dir
// fastq_input_filename
// start_tag end_tag output_filename
// ...
// start_tag end_tag output_filename
//
bool MetaFileData::open(std::string metaFileName)
{
    std::ifstream inFile;
    
    inFile.open(metaFileName);
    if (!inFile) {
        std::cout << "Can't open metafile " << metaFileName << "\n";
        return false;
    }
    
    char master[100];
    
    inFile.getline(master, 100);
    if (!inFile) {
        return false;
    }
    if (strcmp(master,"MasterPrep") == 0){
        _hammingDiag   = false;
    } else if (strcmp(master,"MasterPrepDiag") == 0){
        _hammingDiag   = true;
    } else if (strcmp(master,"SamplePrep") == 0){
        _samplePrep = true;
        _hammingDiag   = false;
    } else if (strcmp(master,"SamplePrepDiag") == 0){
        _samplePrep = true;
        _hammingDiag = true;
    }    
    char name[100];
    
    inFile.getline(name, 100);
    if (!inFile) {
        return false;
    }
    _masterFile = std::string(name);
    
    inFile.getline(name, 100);
    if (!inFile) {
        return false;
    }
    _inputFileName = std::string(name);
	//
	// if building master then done
	//
	if (_buildMaster == true) {
		return true;
    }
    //
    // seq error needs target barcode
    //
    if (_sequenceError) {
        do {
            char bc[100];
            inFile.getline(bc, 100);
            if (!inFile) {
                return true;
            }
            std::string bcS = std::string(bc);
            _barcode.push_back(bcS);
        } while (true);
    }
    //
    // index and start and end tags
    //
    // each line contains [startTag] [endTag] [filename]
    // build a structure to hold tags and barcode storage
    //
    do {
        char line[256];
        inFile.getline(line,256);
        if (!inFile) {
            break;
        }
        
        std::string strLine = std::string(line);
        std::size_t pos = strLine.find(" ");
        std::string indexTag = std::string("");
        std::string startAnchor = strLine.substr(0,pos);
        std::string endAnchor = strLine.substr(pos+1);
        //
        // next space
        //
        pos = endAnchor.find(" ");
        std::string outputFile = endAnchor.substr(pos+1);
        endAnchor = endAnchor.substr(0,pos);
        std::cout << "output file " << outputFile << "\n";
        //
        // if startAnchor = 8 or 10 then there is an index tag
        //
        if (startAnchor.length() == 8) {
            indexTag = startAnchor.substr(0,2);
            startAnchor = startAnchor.substr(2);
        } else if (startAnchor.length() == 10) {
            indexTag = startAnchor.substr(0,2);
            startAnchor = startAnchor.substr(2);
        }
        //
        // build s structure to hold start and end anchor tags and
        // storage vectors for read barcodes
        //
        AnchorTags anchorTags = AnchorTags(indexTag,startAnchor,endAnchor,outputFile);
        //
        // save in vector
        //
        _tagVec.push_back(anchorTags);
        
    } while(true);
    std::cout << "Anchor tags are :\n";
    
    for (auto l : _tagVec) {
        std::cout << l._startIndex << " " <<  l._startAnchor << " - " << l._endAnchor << "\n";
    }
    
    if (_tagVec.size() == 0) {
        return false;
    }
    return true;
}
