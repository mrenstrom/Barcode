
//  main.cpp
//  readBarcode
//
//  Created by mark enstrom on 3/18/17.
//  Copyright © 2017 Mark Enstrom. All rights reserved.
//

#include <iomanip>
#include "barcode.hpp"
#include "IBarcode.hpp"
#include "anchorTag.hpp"
#include "readFile.hpp"
#include "readMeta.hpp"
#include "BuildLibrary.hpp"

/*------------------------------------------------------------------
 *
 * Read in raw barcode file :
 *   remove barcades with low Quality
 *   remove barcodes with bad start index or anchors
 *   de-multiplex barcodes into groups based on index tag
 *   trim start & end anchors to make 20 bit barcode sequence
 *   combine barcodes with no N bases into one dictionary
 *       dictionary combines all exact matches into one entry
 *   combine all barcodes with one or two N bases into dict
 *
 -----------------------------------------------------------------*/
int main(int argc, const char * argv[]) {
    
    if (argc < 2) {
        std::cout << "Useage readBarcodeSample metaFile" << "\n";
        exit(0);
    }
    //
    // get meta parameters
    //
    MetaFileData meta = MetaFileData();
    if (!meta.open(std::string(argv[1]))) {
        std::cout << "metafile not correct\n";
        exit(0);
    }
	//
	// building library?
	//
	if (meta._buildMaster) {
		BuildLibrary bldLib = BuildLibrary();
		bldLib.init(meta._inputFileName);
		bldLib.save(meta._masterFile);
		exit(0);
	}
    //
    // init library
    //
    BarcodeMaster master = BarcodeMaster();
    if (meta._useMaster) {
        if (!master.initFromFile(meta._masterFile)) {
            std::cout << "Failed to init from library " << argv[2] << "\n";
            exit(0);
        }
    }
    //
    // read arrays of matching index/start...end anchor into meta._tagVec
    //
    if (!readFastq(meta)) {
        exit(0);
    }
    //
    // for barocdes in each index tag:
    //
    for (auto &anchors : meta._tagVec)
    {
        anchors.displayStats();
        //
        // prepare sample - means Z09132 and Z08103 sample
        // other option is virus/plasmid prep with many more barcodes so
        // only correct to hamming 3
        //
        if (meta._samplePrep == true) {
            anchors.prepareSample(5);
            anchors.writeOutputPre();
            anchors.writeLog();
        } else {
            anchors.prepareSample(3);
            anchors.writeOutputPre();
            anchors.writeLog();
        }
    }
}
