//
//  readLibrary.hpp
//  readBacodeTrim
//
//  Created by mark enstrom on 12/20/17.
//  Copyright © 2017 Mark Enstrom. All rights reserved.
//

#ifndef readLibrary_hpp
#define readLibrary_hpp
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include "vector"
#include "IBarcode.hpp"

std::vector<IBarcode> *readLibrary(std::string libName);


#endif /* readLibrary_hpp */

