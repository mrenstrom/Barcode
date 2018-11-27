//
//  KMer5.cpp
//  barcodeP2
//
//  Created by mark enstrom on 12/18/17.
//  Copyright © 2017 mark enstrom. All rights reserved.
//

#include "KMer5.hpp"
#include <set>
#include <algorithm>

int hammingDistanceNormal(const char* p1,const char*  p2,int l)
{
    int dist = 0;
    for (int i = 0; i < l; i++)
    {
        if (p1[i] != p2[i]) dist += 1;
        if (dist >= 5) break;
    }
    return dist;
}
//typedef std::map<int,std::vector<int>> KMER_POS;
//typedef std::map<std::string,KMER_POS> KMER_MAP;
/*-------------------------------------------------------------------------------------
 * save lib into text file
 *
 *
 *-------------------------------------------------------------------------------------*/
void KMER_OBJ::save(std::string slog) {
    std::ofstream masterF;
    masterF.open (slog);
	for (auto &seq : _qList) {
		masterF << "#" << seq << "\n";
	}
    for (auto &pair : *_pMap) {
        std::string mer = pair.first;
        KMER_POS &km = pair.second;
        // for each mer 4 possible starting positions
        masterF << mer << "\n";
        for (auto &merPair : km) {
            int pos = merPair.first;
            std::vector<int> &vec = merPair.second;
            for (int i : vec) {
                masterF << pos << " " << i << "\n";
            }
        }
    }
    masterF.close();
}
/*-------------------------------------------------------------------------------------
 * restore lib from text file
 *
 * mer-seq
 * pos seq
 *              AAAAA
 *               0 20628
 *               5 32456
 *              10 23245
 *              15 234234
 *              AAAAC
 *
 * typedef std::map<int,std::vector<int>> KMER_POS;
 * typedef std::map<std::string,KMER_POS> KMER_MAP
 *-------------------------------------------------------------------------------------*/
bool KMER_OBJ::load(std::string base) {
	_qList.clear();
    std::string slog = base + std::string("MasterMap.txt");
    std::ifstream masterF;
    masterF.open (slog);
    if (!masterF) {
        std::cout << "Can't open MKER file " << slog << "\n";
        return false;
    }
    //
    // assume starting AAAAA
    // assume starting loc = 0  (of 0,5,10,15)
    //
    int         loc = -1;
    KMER_POS *merMap = nullptr;
    std::vector<int> *merPos = nullptr;
    do {
        char line[100]{};
        masterF.getline(line, 100);
        if (!masterF) {
            break;
        }
		if (line[0] == '#') {
			std::string l = std::string(line);
			std::string seq = l.substr(1);
			_qList.push_back(seq);
			continue;
		}
        //
        // new mer
        //
        if (isalpha(line[0])) {
            merMap = &_pMap->at(line);
            loc = -1;
        } else {
            //
            // sequence line 0 1232134345
            //
            std::string sLine = std::string(line);
            std::string sPos = sLine.substr(0,sLine.find(' '));
            std::string sSeq = sLine.substr(sLine.find(' ')+1);
            int newPos = stoi(sPos);
            int seq = stoi(sSeq);
            if (newPos != loc) {
                loc = newPos;
                merPos = &merMap->at(loc);
            }
            merPos->push_back(seq);
        }
    } while (true);
    masterF.close();
    
    return true;
}
/*-------------------------------------------------------------------------------------
 * init blank 5-mer map for a 20 bp barcode
 *
 *   5-mer alinged to position 0,5,10,15
 *       [AAAAA]
 *               0: <list> - list of all barcodes in library that have AAAAA at pos[0]
 *            1: <list>
 *            2: <list>
 *            3: <list>
 *        [AAAAC]
 *        .
 *        .
 *        [TTTTT]
 *-------------------------------------------------------------------------------------*/
KMER_OBJ::KMER_OBJ () {
    
    _pMap = new KMER_MAP();
    
    char mers[] = {'A','C','G','T'};
    for (auto c1 : mers) {
        for (auto c2 : mers) {
            for (auto c3: mers){
                for (auto c4: mers){
                    for (auto c5: mers){
                        char mer[] = {c1,c2,c3,c4,c5,'\0'};
                        KMER_POS kp = KMER_POS();
                        
                        for (int i = 0; i < 20; i+=5) {
                            kp.insert(std::pair<int,std::vector<int>>(i,std::vector<int>()));
                        }
                        _pMap->insert(std::pair<std::string,KMER_POS>(mer,kp));
                    }
                }
            }
        }
    }
}

KMER_OBJ::~KMER_OBJ() {
    if (_pMap != nullptr) {
        delete _pMap;
        _pMap = nullptr;
    }
}
/*-------------------------------------------------------------------------------------
 * fill in blank 5-mer map with list of barcodes in qList
 *
 *
 *
 *-------------------------------------------------------------------------------------*/
void KMER_OBJ::buildMap(std::vector<std::string>qList)
{
	std::cout << "KMER build map size = " << qList.size() << "\n";
    //
    // build KMER db
    //
    for (int iq = 0; iq < qList.size(); iq++)
    {
		//
		// save seq-index
		//
		_qList.push_back(qList[iq]);
		//
		// build
		//
        const char * seq = qList[iq].c_str();
        for (int i = 0; i < 20;i+=5)
        {
            char mer[6] = {seq[i],seq[i+1],seq[i+2],seq[i+3],seq[i+4],'\0'};
            KMer5 km;
            km.index = iq;
            km.position = i;
            
            try
            {
                //
                // this could fail with an N
                //
                auto& merMap = _pMap->at(mer);
                auto firstVec  = merMap.at(i);
                //
                // store seq ID in list for mer at i
                //
                merMap[i].push_back(iq);
                //
                // check
                //
                auto checkMap = _pMap->at(mer);
                auto tagVec  = checkMap.at(i);
            }
            catch (std::out_of_range)
            {
                std::cout << "An exception occurred in Build. Exception Nr. " << mer << '\n';
                std::cout << "Error " << mer << '\n';
            }
        }
    }
    //
    // sort lists
    //
    for (auto & item : *_pMap){
        KMER_POS &km = item.second;
        for (int i = 0; i < 20;i+=5)
        {
            std::vector<int> &lv = km.at(i);
            std::sort(lv.begin(),lv.end());
        }
    }
}

/*--------------------------------------------------------------------------------------------
 * merComp
 *
 *  compares a barcode to a table of 5-mers. All barcodes in library where at least 2 5-mers
 *  match are compared for full hamming dist. Barcodes with 1 hamming dist are added to list
 *  for later merge
 *
 * match to seq using kmer table
 * then calculate actual hamming dist to qList
 *
 *--------------------------------------------------------------------------------------------*/

std::vector<int> KMER_OBJ::merCompV(std::string seq,int maxDist)
{
    assert(maxDist <= 2);
	

    std::vector<SeqMatch> matchList = std::vector<SeqMatch>();
    
    if (dbg4 == 1) {
        std::cout << "merMap " << seq << "\n";
    }
    std::map<int, int> fred = std::map<int, int>();
    std::vector<int> fredList[4];
    //
    // record all matches of seq (query barcode) to 5-mers
    //
    for (int i = 0; i < 20; i+=5)
    {
        const char mer[6] = {seq[i],seq[i+1],seq[i+2],seq[i+3],seq[i+4],'\0'};
        try
        {
            auto& kit = _pMap->at((char *)(&mer[0]));
            auto& tagList = kit.at(i);
            if (dbg4 == 1) {
                std::cout << "tagList count = " << tagList.size() << "\n";
            }
            // which ones match the current position?
            //
            // list of barcodes with this mer at this position
            //
            fredList[i/5] = tagList;
        }
        catch (std::out_of_range)
        {
            if (dbg4 == 1)
                std::cout << "An exception occurred. Exception Nr. " << mer << '\n';
        }
    }
    //
    // calc interseection of any 2 5-mers
    //
    if (dbg4 == 1) {
        for (auto item : fredList) {
            std::cout << item.size() << "\n";
        }
    }
    
    std::vector<int> v0 = fredList[0];
    std::vector<int> v1 = fredList[1];
    std::vector<int> v2 = fredList[2];
    std::vector<int> v3 = fredList[3];
    
    std::vector<int> v_intersection10;
    std::vector<int> v_intersection20;
    std::vector<int> v_intersection30;
    std::vector<int> v_intersection12;
    std::vector<int> v_intersection13;
    std::vector<int> v_intersection23;
    
    std::vector<int> v_intersection;
    
    std::set_intersection(fredList[0].begin(), fredList[0].end(),
                          fredList[1].begin(), fredList[1].end(),
                          std::back_inserter(v_intersection10));
    
    std::set_intersection(v0.begin(), v0.end(),
                          v2.begin(), v2.end(),
                          std::back_inserter(v_intersection20));
    
    std::set_intersection(v0.begin(), v0.end(),
                          v3.begin(), v3.end(),
                          std::back_inserter(v_intersection30));
    
    std::set_intersection(v1.begin(), v1.end(),
                          v2.begin(), v2.end(),
                          std::back_inserter(v_intersection12));
    
    std::set_intersection(v1.begin(), v1.end(),
                          v3.begin(), v3.end(),
                          std::back_inserter(v_intersection13));
    
    std::set_intersection(v2.begin(), v2.end(),
                          v3.begin(), v3.end(),
                          std::back_inserter(v_intersection23));
    
    std::set<int> retSet = std::set<int>();
    
    for (auto v : v_intersection10) {
        retSet.insert(v);
    }
    for (auto v : v_intersection20) {
        retSet.insert(v);
    }
    for (auto v : v_intersection30) {
        retSet.insert(v);
    }
    for (auto v : v_intersection12) {
        retSet.insert(v);
    }
    for (auto v : v_intersection13) {
        retSet.insert(v);
    }
    for (auto v : v_intersection23) {
        retSet.insert(v);
    }
    
    if (dbg4 == 1) {
        std::cout << "Set size = " << retSet.size() << "\n";
    }
    //
	// !!! could return list of matches (so caller can decide which merge)
	//
	std::vector<int> rValue = std::vector<int>();
	
    for (auto iSeq2 : retSet) {
        int n = hammingDistanceNormal(seq.c_str(), _qList[iSeq2].c_str(),20);
        if (dbg4 == 1) {
            std::cout << seq << "\n";
            std::cout << _qList[iSeq2] << "\n";
            std::cout << "Hamming  = " << n << "\n";
            
        }
        //
        // only match hamming 1
        //

        if (n <= maxDist) {
            
            if (dbg4 == 1) {
                
                std::cout << "---------------------\n";
                std::cout << iSeq2 << " " << n << "\n";
                std::cout << seq << "\n";
                std::cout << _qList[iSeq2] << "\n\n";
            }
			rValue.push_back(iSeq2);
        }
    }
    return rValue;
}


int KMER_OBJ::merComp(std::string seq,int maxDist)
{
	auto v = merCompV(seq, maxDist);
	if (v.size() == 0) {
		return -1;
	}
	return v[0];
	
}

int KMER_OBJ::merCompMax(std::string seq,int maxDist,const std::vector<IBarcode> &bar)
{
	auto v = merCompV(seq, maxDist);
	if (v.size() == 0) {
		return -1;
	}
	//
	// find max
	//
	int max = 0;
	int mi  = -1;
	for (int i = 0; i < v.size(); i++) {
		int index = v[i];
		auto &ib = bar[index];
		if (ib._seqCount > max) {
			max = (int)ib._seqCount;
			mi = index;
		}
	}
	return mi;
}
