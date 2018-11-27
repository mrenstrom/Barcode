//
//  BuildLibrary.cpp
//  processBarcodeSample
//
//  read all the masterLibrary prep files and build into one library
//
//
//  Created by mark enstrom on 1/9/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//

#include "BuildLibrary.hpp"
#include "iomanip"
#include "math.h"

bool BuildLibrary::init(std::string masterListFile)
{
	//
	// open master list file
	//
	std::ifstream inFile;
	
	inFile.open(masterListFile);
	if (!inFile) {
		std::cout << "Can't open library master list file " << masterListFile << "\n";
		return false;
	}
	//
	// every entry is a library prep file
	//
	do {
		char prepName[1000];
		inFile.getline(prepName, 1000);
		if (!inFile) {
			break;
		}
		//
		// add file to prep
		//
		addFileToLib(std::string(prepName));
	} while(true);
	
	std::cout << "Lib size = " << _exactMap.size() << "\n";
	
	return true;
}

std::vector<std::string> split(const std::string &s, char delim) {
	std::stringstream ss(s);
	
	std::string item;
	std::vector<std::string> tokens;
	while (getline(ss, item, delim)) {
		tokens.push_back(item);
	}
	return tokens;
}

bool BuildLibrary::addFileToLib(std::string filename)
{
	std::cout << filename << "\n";
	
	std::ifstream inFile;
	
	inFile.open(filename);
	if (!inFile) {
		std::cout << "Can't open library prep file " << filename << "\n";
		return false;
	}
    char line[1000];
    //
    // skip first line
    //
    inFile.getline(line, 1000);
    if (!inFile) {
        return false;
    }
	do {
		inFile.getline(line, 1000);
		if (!inFile) {
			break;
		}
		std::vector<std::string> lnv = split(std::string(line),',');
		//
		// format is:
		// [barcode],[count],[minHamming],[index]
		//
		//std::cout << lnv[0] << " " << lnv[1] << "\n";
		std::string seq = lnv[0];
		std::string fred = lnv[1];
		unsigned long count = (unsigned long)stoi(lnv[1]);
		
		int nbases = 0;
		for (int i = 0; i < seq.size(); i++)
		{
			if (seq[i] == 'N') {
				nbases++;
			}
		}
		//
		// add to appropriate dictionary
		//
		if (nbases == 0) {
			try
			{
				//
				// this will throw out_of_range if seq not in master
				//
				auto oldCount = _exactMap.at(seq);
				_exactMap[seq] = oldCount + count;
				//std::cout << seq << " add " << count << " to existing " << oldCount << "\n";
			}
			catch (std::out_of_range)
			{
				_exactMap[seq] = count;
				//std::cout << seq << " new " << count << "\n";
			}
		} else {
			try
			{
				//
				// this will throw out_of_range if seq not in master
				//
				auto oldCount = _nMap.at(seq);
				_nMap[seq] = oldCount + count;
				//std::cout << seq << " add N " << count << " to existing " << oldCount << "\n";
			}
			catch (std::out_of_range)
			{
				_nMap[seq] = count;
				//std::cout << seq << " new N " << count << "\n";
			}
		}
		
	} while(true);
	return true;
}

struct myMasterSort {
    bool operator() (IBarcode i,IBarcode j) { return (i._seqCount > j._seqCount);}
} masterCompObj;


bool BuildLibrary::save(std::string masterBase) {
	std::string listName = masterBase + std::string("MasterBarcodeList.txt");
	std::string mapName = masterBase + std::string("MasterMap.txt");
	std::cout << "Save master library to " << listName << "\n";
	std::cout << "Save master map    to " << mapName << "\n";
	
	std::ofstream masterList;
	masterList.open (listName);
	//
    // build list from dictionary then sort
    //
    std::vector<IBarcode> vecSort = std::vector<IBarcode>();
    int totalCount = 0;
    for (auto &it : _exactMap) {
        totalCount += it.second;
    }
    //
    // this threshold is very low!!!
    //
    //int threshold = int(round(((float)totalCount / 100000000.0)));
    int threshold = 3;
    std::cout << "total barcode count =  " << std::setw(10) << totalCount << "\n";
    std::cout << "threshold           =  " << std::setw(10) << threshold << "\n";
    //
    // only keep barcodes above threshold
    //
	for (auto &it : _exactMap) {
        IBarcode ib = IBarcode(it.first,it.second);
        if (ib._seqCount >= threshold) {
            vecSort.push_back(ib);
        }
	}
    //
    //
    //
    std::sort(vecSort.begin(),vecSort.end(),masterCompObj);
    //
    // finish file
    //
    masterList << "Barcode" << "," << "Count" << "\n";
    for (auto &ib:vecSort) {
            masterList << ib.barcode() << "," << ib._seqCount << "\n";
    }
	masterList.close();
	//
	// build final KMER and save
	//
	std::cout << "Build Map" << "\n";
	std::vector<std::string> qList = std::vector<std::string>();
	for (auto &ib : vecSort) {
		qList.push_back(ib.barcode());
	}
	//
	// build lookup obj
	//
	KMER_OBJ kmer = KMER_OBJ();
	kmer.buildMap(qList);
	kmer.save(mapName);
	
	
	return true;
}

