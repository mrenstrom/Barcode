//
//  BuildLibrary.hpp
//  processBarcodeSample
//
//  Created by mark enstrom on 1/9/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//

#ifndef BuildLibrary_hpp
#define BuildLibrary_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <unordered_map>
#include "barcode.hpp"
#include "KMer5.hpp"
#include "IBarcode.hpp"

std::vector<std::string> split(const std::string &s, char delim);

class BuildLibrary {
public:
	bcodeMap    _exactMap;
	bcodeMap    _nMap;
	bool init(std::string masterListFile);
	bool save(std::string masterFile);
	bool addFileToLib(std::string filename);
	BuildLibrary(){
		_exactMap    = bcodeMap();
		_nMap        = bcodeMap();
	};
};

#endif /* BuildLibrary_hpp */
