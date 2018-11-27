//
//  readLibrary.cpp
//  readBacodeTrim
//
//  Created by mark enstrom on 12/20/17.
//  Copyright © 2017 Mark Enstrom. All rights reserved.
//

#include "readLibrary.hpp"
#include "IBarcode.hpp"

std::vector<IBarcode> *readLibrary(std::string libName)
{
    std::ifstream inFile;
    inFile.open(libName);
    if (!inFile) {
        std::cout << "Can't open master library file " << libName << "\n";
        return nullptr;
    }
    //
    // Verify first line Barcode,Count
    //
    char line[100]{};
    inFile.getline(line, 100);
    if (strcmp(line,"Barcode,Count") != 0) {
        std::cout << "barcode master file not verified" << "\n";
        std::cout << "First line = " << line << "\n";
        return nullptr;
    }
               
    
    
    std::vector<IBarcode> *pVec = new std::vector<IBarcode>();

    do {
        
        char bases[21]{};
        char sn1[100]{};

        
        inFile.getline(bases, 100,',');
        if (!inFile) {
            break;
        }
        inFile.getline(sn1, 20);
		
        std::string seq = std::string(bases);
        int seqCount = std::stoi(sn1);
        IBarcode ibc = IBarcode(seq, seqCount);
        pVec->push_back(ibc);
    } while (true);
    return pVec;
}

