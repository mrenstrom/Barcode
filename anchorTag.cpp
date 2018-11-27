//
//  anchorTag.cpp
//      Store all barcodes for a given index tag. Provide routines to match barcodes to
//      master list and also to error-correct barcodes
//
//  Created by mark enstrom on 1/6/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//
#include "barcode.hpp"
#include "IBarcode.hpp"
#include "anchorTag.hpp"
#include "iomanip"
//
// sort IBarcode helper class
//
struct myclassS {
    bool operator() (IBarcode i,IBarcode j) { return (i._seqCount < j._seqCount);}
} smallCompObj;

struct myclassB {
    bool operator() (IBarcode i,IBarcode j) { return (i._seqCount > j._seqCount);}
} bigCompObj;

/*--------------------------------------------------------------------------------------------
 * AnchorTags::findClose(IBarcode &ib, std::vector(<IBarcode>) sortList,int start)
 *
 *  build a mapping object with larger count barcodes
 *
 *  match low count barcodes
 *
 *--------------------------------------------------------------------------------------------*/
int AnchorTags::findClose(IBarcode &ib, std::vector<IBarcode> sortList,int start,unsigned long maxHamming)
{
    //
    // start at largest barcode and search down
    //
    for (int i = (int)(sortList.size() - 1); i > start; i--) {
        IBarcode &ib2 = sortList[i];
        unsigned long h = ib2.hammingDist(ib);
        if (h <= maxHamming) {
            return i;
        }
    }
    return -1;
}
/*--------------------------------------------------------------------------------------------
 * AnchorTags::prepareSample
 *
 *  1: add all barcode counts > threshold (1 or 2) to initial "good" list
 *  2: build mapping obj from good list
 *  3: try to match all under threshold barcodes to "good" list (1-2 hamming)
 *  4: try to match all "N" barcode to "good" list (1-2 hamming)
 *  5:
 *
 *--------------------------------------------------------------------------------------------*/
void AnchorTags::prepareSample(int finalHamming)
{
    std::cout << "-----------------------------------------------------\n";
    std::cout << _fileName << "\n";
    std::cout << "Prepare Sample Final Count              = " << std::setw(8) << _exactMap.size() << "\n";
    std::vector<std::string> qList= std::vector<std::string> ();
    std::vector<IBarcode> vecSmall = std::vector<IBarcode>();
    std::vector<IBarcode> vecGood = std::vector<IBarcode>();
    bcodeMap _resultMap = bcodeMap();
    unsigned long threshold1 = 1;  // keep count over this value
    unsigned long threshold2 = 2;  // merge barcodes within
    //
    // get barcodes with count > threshold
    //
    for (auto &it : _exactMap) {
        const std::string &seq = it.first;
        unsigned long count = it.second;
        IBarcode ib = IBarcode(seq,count);
        
        if (count > threshold1) {
            vecGood.push_back(ib);
        } else {
            vecSmall.push_back(ib);
        }
    }
    //
    // build fast-mapping obj   !!! vecGood must not be re-ordered
    // whem using the mapping object!!!
    //
    std::sort(vecGood.begin(),vecGood.end(),smallCompObj);
    for (auto &ib:vecGood) {
        qList.push_back(ib.barcode());
    }
    KMER_OBJ kmer = KMER_OBJ();
    kmer.buildMap(qList);
    //
    // try to match low barcode count barcodes to higher count ones
    //
    std::cout << "try to correct small count barcodes     = " << std::setw(8) << vecSmall.size() << "\n";
    int found = 0;
    int notFound = 0;
    for (auto &smallb: vecSmall) {
        int index = kmer.merComp(smallb.barcode(), (int)threshold2);
        if (index != -1) {
            //
            // merge with big
            //
            auto &ibBig = vecGood[index];
            MergeLogEntry me = MergeLogEntry(smallb.barcode(),smallb._seqCount,ibBig.barcode(),ibBig._seqCount);
            me.hd = 2;
            _mergeLog.push_back(me);
            ibBig._seqCount += smallb._seqCount;
            found++;
        } else {
            //
            // no match, discard
            //
            //
            // maybe need to add these back!!!
            //
            _vecDiscard.push_back(smallb);
            notFound ++;
        }
    }
    std::cout << "Found                                     "  << std::setw(8) << found << " error correcting small barcodes\n";
    std::cout << "Not Found                                 " << std::setw(8) << notFound << " small barcode sequences\n";
    //
    // try to error correct "N" barcodes to good barcodes
    //
    std::cout << "look for N barcodes, count =              " << std::setw(8) << _nMap.size() << "\n";
    found = 0;
    notFound = 0;
    for (auto &it : _nMap) {
        const std::string &seq = it.first;
        unsigned long count = it.second;
        int index = kmer.merComp(seq, (int)threshold2);
        if (index != -1) {
            //
            // match, add count to known seq
            //
            auto &ibBig = vecGood[index];
            MergeLogEntry me = MergeLogEntry(seq,count,ibBig.barcode(),ibBig._seqCount);
            me.hd = 2;
            _mergeLog.push_back(me);
            ibBig._seqCount += count;
            found++;
        } else {
            // no match, discard
            //
            //
            IBarcode ib = IBarcode(seq,count);
            // really discard..._vecDiscard.push_back(ib);
            notFound ++;
        }
    }
    std::cout << "Found                                     " << std::setw(8) << found << " error correcting barcodes\n";
    std::cout << "Not Found                                 " << std::setw(8) << notFound << " N barcode sequences\n";
    //
    // now starting from smallest count barcode and working up: look
    // for close barcodes and merge any within (1-2) to the close
    // barcode with the largest count
    //
    // since order of index changes due to sort, a new KMER object must be build
    //
    qList.clear();
    std::sort(vecGood.begin(),vecGood.end(),smallCompObj);
    for (auto &ib:vecGood) {
        qList.push_back(ib.barcode());
    }
    KMER_OBJ kmer2 = KMER_OBJ();
    kmer2.buildMap(qList);
    
    std::cout << "Now search all barcodes within 2         " << std::setw(8) << vecGood.size() << "\n";
    int merged = 0;
    for (int iGood = 0; iGood < vecGood.size();iGood++) {
        auto &ib = vecGood[iGood];
        if (ib._seqCount > 100000) {
            break;
        }
        int index = kmer2.merCompMax(ib.barcode(), (int)threshold2,vecGood);
        if ((index != -1) && (index != iGood)) {
            auto &ibBig = vecGood[index];
            MergeLogEntry me = MergeLogEntry(ib.barcode(),ib._seqCount,ibBig.barcode(),ibBig._seqCount);
            //
            // debug entry
            //
            me.hd = ibBig.hammingDist(ib);
            _mergeLog.push_back(me);
            ibBig._seqCount = ibBig._seqCount + ib._seqCount;
            ib._valid = false;
            merged++;
        }
    }
    std::cout << "merged                                  = " << std::setw(8) << merged << "\n";
    std::cout << "Search for hamming " << finalHamming << "\n";

    //
    // add unmatched singles in  
    //
    for (auto ib:_vecDiscard) {
      vecGood.push_back(ib);
      //std::cout << "add back in 1 bc " << ib._bc << "\n";
    }
    std::sort(vecGood.begin(),vecGood.end(),smallCompObj);
    //
    // one more search: from largest down...merge from smallest up
    //
    int loopCount = 0;
    for (int iBig = (int)vecGood.size()-1; iBig > 0;iBig--) {
        auto &ibBig = vecGood[iBig];
        //std::cout << "BIG " << iBig << " " << ibBig._seqCount << " " << ibBig._valid << "\n";
        if (ibBig._valid) {
            for (int iSmall = 0; iSmall < iBig;iSmall++) {
                if (iSmall == iBig) continue;
                auto &ibSmall = vecGood[iSmall];
                if (ibSmall._valid) {
                    unsigned long h = ibBig.hammingDist(ibSmall);
                    if (h <= finalHamming) {
                        //std::cout << "Merge d = " << h << " " << ibSmall._seqCount << " > " << ibBig._seqCount << "\n";
                        //std::cout << ibSmall._bc << "\n";
                        //std::cout << ibBig._bc << "\n";
                        //
                        // merge small to big
                        //
                        MergeLogEntry me = MergeLogEntry(ibSmall.barcode(),ibSmall._seqCount,ibBig.barcode(),ibBig._seqCount);
                        //
                        // debug entry
                        //
                        me.hd = ibBig.hammingDist(ibSmall);
                        _mergeLog.push_back(me);
                        ibBig._seqCount = ibBig._seqCount + ibSmall._seqCount;
                        ibSmall._valid = false;
                        merged++;
                    }
                }
            }
        }
        loopCount++;
        if ((loopCount % 1000) == 0) {
          std::cout << "Loop " << loopCount << " of " << vecGood.size() << " m = " << merged << "\n";
        }
    }
    std::cout << "merged                                  = " << std::setw(8) << merged << "\n";
    //
    // build _vecFound after removing invalid (merged) barcodes.
    //
    for (auto &ib : vecGood) {
        if (ib._valid) {
            _vecFound.push_back(ib);
        }
    }
    std::sort(_vecFound.begin(),_vecFound.end(),bigCompObj);
    
    //
    // vecFound are results that exacly match master lib
    // vecNotFound don't match exactly but may still match with a diff of 1 or 2
    //
    std::cout << "Unique Barcodes                           " << std::setw(8) << _vecFound.size()    << "\n";
    std::cout << "-----------------------------------------------------\n";
}

/*--------------------------------------------------------------------------------------------
 * AnchorTags::writeOutput
 *
 *
 *
 *--------------------------------------------------------------------------------------------*/
void AnchorTags::writeOutput()
{
    std::cout << "write output to " << _fileName << "\n";
    std::ofstream myOut;
    std::string sOut = _fileName;
    myOut.open (sOut);
    std::sort (_vecFound.begin(), _vecFound.end(), bigCompObj);
    myOut << "barcode,count,minHamming,index\n";
    for (auto &ib : _vecFound) {
        myOut << ib.barcode() << "," << ib._seqCount << std::endl;
    }
}
/*--------------------------------------------------------------------------------------------
 * AnchorTags::writeOutput
 *
 *
 *
 *--------------------------------------------------------------------------------------------*/
void AnchorTags::writeOutputPre()
{
    std::ofstream myOut;
    std::string sOut = _fileName;
    myOut.open (sOut);
    std::sort (_vecFound.begin(), _vecFound.end(), bigCompObj);
    myOut << "barcode,count,minHamming,index\n";
    for (auto &ib : _vecFound) {
        myOut << ib.barcode() << "," << ib._seqCount << std::endl;
    }
    myOut.close();
    //
    // write out raw summarized data
    //
    std::ofstream ofRaw;
    std::string rawName = _fileName + ".raw";
    ofRaw.open(rawName);
    ofRaw << "barcode,count\n";
    for (auto &it:_exactMap) {
        ofRaw << it.second << "," << it.first << "\n";
    }
    for (auto &it:_nMap) {
        ofRaw << it.first << "," << it.second << "\n";
    }
    ofRaw.close();
}
/*--------------------------------------------------------------------------------------------
 * AnchorTags::writeLog
 *
 *
 *
 *--------------------------------------------------------------------------------------------*/
void AnchorTags::writeLog()
{
    //
    // log
    //
    std::ofstream myLog;
    std::string slog = _fileName + ".log";
    std::cout << "Writing output to " << slog << "\n";
    myLog.open (slog);
    myLog << "#Lib 1 initial exact size       = " << _exactMap.size() << "\n";
    myLog << "#Lib 1 initial N size           = " << _nMap.size() << "\n";
    myLog << "#------------------------------------------------\n";
    myLog << "# good barcodes                 = " << _vecFound.size() << "\n";
    myLog << "#discarded barcodes             = " << _vecDiscard.size() << "\n";
    myLog << "#merged barcodes                = " << _vecMerge.size() << "\n";
    myLog << "#---------------------------------\n";
    for (auto &ml : _mergeLog) {
        myLog << ml.bc1 << " " << ml.count1 << " >> " << ml.bc2 << " " << std::setw(7) << ml.count2 << " hd " << ml.hd << "\n";
    }
    for (auto &ib : _vecDiscard) {
        myLog << "discard " << ib.barcode() << "\n";
    }
    myLog.close();
}

/*--------------------------------------------------------------------------------------------
 * AnchorTags::displayStats
 *      master
 *
 *
 *--------------------------------------------------------------------------------------------*/
void AnchorTags::displayStats()
{
    std::cout << "-----------------------------------------------------\n";
    std::cout << "Stats for                  " << _exactMap.size() <<  " barcodes \n";
    unsigned long countArray[15] = {0};
    std::string labelArray[15] = {"    0","    1","    2","    3","    4","    5",
        "    6","    7","    8","    9","   10","  100"," 1000","10000"," more"};
    
    for (auto &it : _exactMap) {
        unsigned long count = it.second;
        if (count <= 10) {
            countArray[count]++;
        }
        else if (count < 100) {
            countArray[11]++;
        }
        else if (count < 1000) {
            countArray[12]++;
        }
        else if (count < 10000) {
            countArray[13]++;
        } else {
            countArray[14]++;
        }
    }
    //
    // vecFound are results that exacly match master lib
    // vecNotFound don't match exactly but may still match with a diff of 1 or 2
    //
    for (int i = 1;i <= 14; i++) {
        std::cout << labelArray[i]  << " count = " << std::setw(8) << countArray[i] << "\n";
    }
    std::cout << "-----------------------------------------------------\n";
}
