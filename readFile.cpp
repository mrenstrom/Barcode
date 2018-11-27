//
//  readFile.cpp
//  Created by mark enstrom on 12/3/17.
//  Copyright © 2017 mark enstrom. All rights reserved.
//
#include "barcode.hpp"
#include "IBarcode.hpp"
#include "anchorTag.hpp"
#include "readMeta.hpp"
#include "readFile.hpp"
#include "iomanip"

struct myclassSX {
    bool operator() (IBarcode i,IBarcode j) { return (i._seqCount < j._seqCount);}
} smallCompObjX;

/*--------------------------------------------------------------------------------------------
 *
 * like strncmp
 * compare a dna(char) sequence and return # differences
 *
 *
 *--------------------------------------------------------------------------------------------*/
int dnancmp(char *p1, char *p2,size_t n) {
    int error = 0;
    for (int i = 0; i < n; i++) {
        if (p1[i] != p2[i]) {
            error += 1;
        }
    }
    return error;
}
/*--------------------------------------------------------------------------------------------
 * FASTQ Read Format
 * [index tag][start anchor][20 bases][end anchor]
 *
 * read fastq seq and qual lines
 *  - mark bases with qual score < 30 as N
 *  - any line with 3 or more Ns is discarded
 *  - find start anchor (one substutution allowed)
 *  - find end anchor at 20bp beyond start anchor (total one error in start + end)
 *  - cut out 20bp barcode
 *  - if barcode has no "N" store in exact unordered_map
 *    - otherwise store in N unordered map
 *
 *--------------------------------------------------------------------------------------------*/

bool readFastq(MetaFileData &meta)
{
    std::ifstream inFile;
    inFile.open(meta._inputFileName);
    if (!inFile) {
        std::cout << "Can't open input file " << meta._inputFileName << "\n";
        return false;
    }
    std::cout << "Read from source file " << meta._inputFileName << "\n";
    AnchorTagVec &tagVec                      = meta._tagVec;
    std::vector<std::string> &failQual        = meta._failQual;
    
    int reads = 0;
    clock_t cl = clock();
    clock_t cl1 = clock();
    //
    // read in all data.
    // keep track of metrics
    //
    int goodQuality   = 0;
    int badQuality    = 0;
    int goodAnchor    = 0;
    int badAnchor     = 0;
    int goodIndex     = 0;
    int Nbarcode      = 0;
    int exactBarcode  = 0;
    
    do {
        reads ++;
        bool track = false;
        char hdr[100];
        char oriSeq[100];
        char seq[100];
        char mid[2];
        char qual[100];
        
        inFile.getline(hdr, 100);
        if (!inFile) {
            break;
        }
        inFile.getline(seq, 100);
        if (!inFile) {
            break;
        }
        inFile.getline(mid, 2);
        if (!inFile) {
            break;
        }
        inFile.getline(qual, 100);
        if (!inFile) {
            break;
        }
        strcpy(oriSeq,seq);
        std::string sub1 = std::string(seq);
        sub1 = sub1.substr(0,10);
        //
        // debug tracking
        //
        //        if (sub1 == std::string("CATTTCTAGA")) {
        //            track = true;
        //        }
        //
        // mark bases with qual < (N Threshold) as "N"
        //
        int badBases = 0;
        for (int i = 0; i < (strlen(qual)-1); i++)
        {
            int qScore = (int)qual[i] - 33;
            if (qScore < 20)
            {
                seq[i] = 'N';
            }
            //
            // some bases are originally 'N' from fastq so count
            // separately
            //
            if (seq[i] == 'N') {
                badBases++;
            }
        }
        //
        // allow 3 N bases, a good sequence could have 2 N in 20bp barcode and
        // one N in anchor seq.
        //
        if (badBases > 3) {
            badQuality += 1;
            if (track) {
                std::cout << "fail quality " << seq << "\n";
                std::cout << "Q          = " << qual << "\n";
            }
            failQual.push_back(oriSeq);
            failQual.push_back(qual);
            //
            // done with this sequence
            //
            continue;
        }
        //
        // good qaulity
        //
        goodQuality += 1;
        //
        // look for exact match on index (if present)
        // allow 1 mismatch in start+end anchor
        //
        bool bHasMatchingIndex = false;
        for (auto &anchors : tagVec)
        {
            std::string indexTag    = anchors._startIndex;
            std::string startAnchor = anchors._startAnchor;
            std::string endAnchor   = anchors._endAnchor;
            //
            // if there is an index tag it must match exactly
            //      if an error correction is made it could change index
            //
            // there may not be a indexTag if multiple samples were not combined
            //
            if (indexTag.length() > 0) {
                int startIndexError = dnancmp((char *)indexTag.c_str(),(char *) seq,indexTag.length());
                // index fail, may match other index
                if (startIndexError != 0) {
                    //skip to next anchor
                    continue;
                }
                // index match
                bHasMatchingIndex = true;
                goodIndex++;
                anchors.totalIndexReads++;
            }
            //
            // can have one/zero error in startTag + endTag
            //
            size_t startTagError = dnancmp((char *)startAnchor.c_str(),
                                           &seq[indexTag.length()],
                                           startAnchor.length());
            
            size_t endTagError   = dnancmp(&seq[indexTag.length() + startAnchor.length() + 20],
                                         (char *)endAnchor.c_str(),
                                         endAnchor.length());
            //
            // start/end anchor error
            //
            if ((startTagError+endTagError) > 0) {
                //
                // was index matched?
                //
                if (bHasMatchingIndex) {
                    // index specific
                    anchors._failAnchor.push_back(seq);
                    badAnchor++;
                    if (track) {
                        std::cout << "fail start+end anchor " << startAnchor.c_str() << "\n";
                    }
                    // sequence has bad anchors, done with sequence
                    break;
                }

                // was a non-index anchor but either bad anchor or index archor line
                // skip to next anchor
                continue;
            }
            if (track) {
                std::cout << "match start + end " << startAnchor.c_str() << "\n";
            }
            goodAnchor++;
            //
            // inc total reads for non-index anchor match
            //
            if (indexTag.length() == 0) {
                anchors.totalIndexReads++;
            }
            //
            // pull out 20bp barcode and check for N bases
            //
            std::string bc20 = std::string(seq,indexTag.length() + startAnchor.length(),20);
            //
            // were there any Ns in the 20 bp barcode itself (as opposed to anchors)
            //
            int badBarcodeBases = 0;
            for (int i = 0; i < bc20.size(); i++)
            {
                if (bc20[i] == 'N') {
                    badBarcodeBases++;
                }
            }
            //
            // add to appropriate dictionary
            //
            if (badBarcodeBases == 0) {
                ++anchors._exactMap[bc20];
                exactBarcode += 1;
                if (track) {
                    std::cout << "exact count = " << anchors._exactMap.size() << bc20 << "\n";
                }
            } else if (badBarcodeBases <= 2){
                ++anchors._nMap[bc20];
                Nbarcode += 1;
                if (track) {
                    std::cout << "N count = " << anchors._nMap.size() << bc20 << "\n";
                }
            } else {
                anchors._failBCQual.push_back(seq);
                if (track) {
                    std::cout << "3 or more N in seq " << bc20 << "\n";
                }
            }
            //
            // done with this sequence
            //
            break;
        }

        if ((reads%5000000) == 0) {
            clock_t cl2 = clock();
            double lineTime = (cl2-cl1)/(double)CLOCKS_PER_SEC;
            cl1= cl2;
            std::cout << "--------------------------------------------\n";
            std::cout << "reads             " << std::setw(10) << reads << " " << lineTime << "s \n";
        }
        
    } while (true);
    
    clock_t cl2 = clock();
    double lineTime = (cl2-cl)/(double)CLOCKS_PER_SEC;
    std::cout << "-------FINAL Counts---------------------------\n";
    std::cout << "reads             " << std::setw(8) << reads << " " << lineTime << "\n";
    std::cout << "goodQuality       " << std::setw(8) << goodQuality << "\n";
    std::cout << "badQuality        " << std::setw(8) << badQuality << "\n";
    std::cout << "goodAnchor        " << std::setw(8) << goodAnchor << "\n";
    std::cout << "badAnchor         " << std::setw(8) << badAnchor << "\n";
    std::cout << "Nbarcode          " << std::setw(8) << Nbarcode << "\n";
    std::cout << "exactBarcode      " << std::setw(8) << exactBarcode << "\n";
    //
    // summarize each tag
    //
    for (auto &anchors : tagVec)
    {
        std::cout << anchors._startIndex  << "\n";
        std::cout << anchors._startAnchor << "\n";
        std::cout << "total reads for index " << std::setw(8) << anchors.totalIndexReads << "\n";
        std::cout << "index anchor fail     " << std::setw(8) << anchors._failAnchor.size() << "\n";
        std::cout << "index bc qual fail    " << std::setw(8) << anchors._failBCQual.size() << "\n";
        std::cout << "unique Exact BCs      " << std::setw(8) << anchors._exactMap.size()  << "\n";
        std::cout << "unique N BCs          " << std::setw(8) << anchors._nMap.size()  << "\n";
        unsigned long count = 0;
        for (auto &it:anchors._exactMap) {
            count = count + it.second;
        }
        for (auto &it:anchors._nMap) {
            count = count + it.second;
        }
        std::cout << "total barcodes        " << std::setw(8) << count << "\n";
    }
    
    
    //
    // diagnostic : compare barcode dist
    //
    if (meta._hammingDiag) {

      std::cout << "Run Hamming distance diagnostic\n";
      //
      // global distances
      //
      std::map<unsigned long,int> distMap = std::map<unsigned long,int>();
      for (auto &anchors : tagVec)
      {
          std::cout << anchors._fileName << "\n";
          //
          // convert map to vector
          //
          std::vector<IBarcode> vecGood = std::vector<IBarcode>();
          //
          // get barcodes w
          //
          for (auto &it : anchors._exactMap) {
              const std::string &seq = it.first;
              unsigned long count = it.second;
              IBarcode ib = IBarcode(seq,count);
              vecGood.push_back(ib);
          }
          //
          // build fast-mapping obj   !!! vecGood must not be re-ordered
          // whem using the mapping object!!!
          //
          std::sort(vecGood.begin(),vecGood.end(),smallCompObjX);
          std::cout << "Calc hamming on vector of size " << vecGood.size() << "\n";          
          std::string path = anchors._fileName;
          size_t l = path.find_last_of("/");

          std::string base = anchors._fileName.substr(l+1);
          std::cout << "Base = " << base << std::endl;
          l = base.find_last_of(".");
          std::string file1 = base.substr(0,l) + "_ham.txt";
          std::string file2 = base.substr(0,l) + "_global_ham.txt";

          std::ofstream myLog;
          std::string slog = meta._masterFile + file1;
          myLog.open (slog);
          myLog << "#" << anchors._startIndex << " " << anchors._fileName << "\n";
          
          for (int i = 0; (i < 1000) && (i < vecGood.size()) ; i++) {
              std::map<unsigned long,int> tagDistMap = std::map<unsigned long,int>();
              // init map
              for (int x = 0; x < 19; x++) {
                tagDistMap[x] = 1;
              }
              IBarcode &ib = vecGood[i];
              for (int j = 0; j < vecGood.size() ; j++) {
                  IBarcode &ib2 = vecGood[j];
                  unsigned long h = ib2.hammingDist(ib);
                  // add to local and global dist map
                  tagDistMap[h] += 1;
                  distMap[h] += 1;
              }
              for (auto pr:tagDistMap) {
                 myLog << i << "," << pr.first << "," << pr.second << "\n";
              }
          }

          myLog.close();
          //
          // global
          //
          slog = meta._masterFile + file2;
          myLog.open(slog);
          for (auto pr:distMap) {
             myLog << pr.first << "," << pr.second << "\n";
          }
          myLog.close();
      }
    }
    return true;
}

