#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "ryggrad/src/base/Logger.h"
#include "RestSiteCoreUnit.h"


int RestSiteMapCore:: CreateRSitesPerString(const string& origString, const string& origName, RSiteReads& reads, bool addRC) const {
  if (origString == "" && origName == "") {
    return 0;
  }
  svec<int> mm;
  RSiteRead rr;
  rr.Name() = origName;
  int motifLen = m_motif.length();
  int origLen  = origString.length();
  mm.reserve(origLen/motifLen);
  bool wrotePrefix = false;
  int n = -1;
  for (int i=0; i<origLen-motifLen ; i++) {
    int j = 0;
    for (j=0; j<motifLen; j++) {
      if (m_motif[j] != origString[i+j])
        break;
    }
    if (j == motifLen) {
      if (n >= 0) {
        // Obtain the pre/post & dmer values
        if (!wrotePrefix) {
          rr.PreDist() = n; // prefix (number of trailing bits before the first motif location)
          wrotePrefix = true; 
        }
        mm.push_back(i-n);
      } 
      n = i;
    }
    if (origString != "") {
      if(wrotePrefix) { 
        rr.PostDist() = origLen - n - 1; // postfix (number of leading bits after last motif location
      } 
    }
  } 
  rr.Dist() = mm;
  int readIdx = reads.AddRead(rr);
  FILE_LOG(logDEBUG3) << "Adding Read: " << readIdx << "  " << rr.Name();
  if(addRC) {  
    rr.Flip();
    rr.Name() += "_RC";
    readIdx = reads.AddRead(rr);
    FILE_LOG(logDEBUG3) << "Adding Read: " << readIdx << "  " << rr.Name();
    return 2*mm.isize(); // Return the total number of sites that have been added
  }
  return mm.isize(); // Return the total number of sites that have been added
}

void RestSiteMapCore::BuildDmers() { 
  int dimCount = pow(TotalSiteCount()*3, 1.0/m_modelParams.DmerLength());   // Number of bins per dimension
  if(dimCount<2) { dimCount = 2; } // implementation ease
  if(pow(dimCount, m_modelParams.DmerLength()) > 1900000000) {              //TODO parameterise
    FILE_LOG(logWARNING) << "Input data size is too large";
    dimCount = pow(1900000000, 1.0/m_modelParams.DmerLength()); 
  }
  FILE_LOG(logINFO) << "Estimated number of Dmers and dimension size for dmer storage: " << TotalSiteCount() << "  " << dimCount; 
  m_dmers.BuildDmers(m_rReads , m_modelParams.DmerLength(), m_modelParams.MotifLength(), dimCount); 
}

int RestSiteMapCore::FindMapInstances(float indelVariance, map<int, map<int,int>>& checkedSeqs) const {
  int counter       = 0;
  double matchCount = 0;
  int loopLim       = m_dmers.NumCells();
  for (int iterIndex=0; iterIndex<loopLim; iterIndex++) {
    counter++;
    if (counter % 100000 == 0) {
//      cout << "\rLOG Progress: " << 100*(double)iterIndex/(double)loopLim << "%" << flush;
    }
    if(!m_dmers[iterIndex].empty()) {
      FILE_LOG(logDEBUG2) << "Number of dmers in cell " << iterIndex << " " << m_dmers[iterIndex].isize(); 
      matchCount += HandleMappingInstance(m_dmers[iterIndex], indelVariance, checkedSeqs, false);
    }
  }
  return matchCount;
}

int RestSiteMapCore::FindSingleReadMapInstances(const RSiteRead& read, int rIdx, float indelVariance, map<int, map<int,int>>& checkedSeqs) const {
  svec<Dmer> dmers;
  m_dmers.GenerateDmers(read, rIdx, dmers);
  return HandleMappingInstance(dmers, indelVariance, checkedSeqs, true);
}

int RestSiteMapCore::HandleMappingInstance(const svec<Dmer>& dmers, float indelVariance, map<int, map<int,int>>& checkedSeqs, bool acceptSelf) const {
  int matchCount = 0;
  svec<int> neighbourCells;
  neighbourCells.reserve(pow(2, m_modelParams.DmerLength()));
  svec<int> deviations;
  deviations.resize(m_modelParams.DmerLength());
  for(Dmer dm1:dmers) {
    dm1.CalcDeviations(deviations, indelVariance, m_modelParams.CNDFCoef()); //TODO this does not need to be redone every time!
    int merLoc = m_dmers.MapNToOneDim(dm1.Data());
    neighbourCells.clear();
    m_dmers.FindNeighbourCells(merLoc, dm1, deviations, neighbourCells); 
    for (int nCell:neighbourCells) {
      for (auto dm2:m_dmers[nCell]) {
        int offset = abs(dm1.Pos() - dm2.Pos());
        map<int, map<int, int>>::iterator it1 = checkedSeqs.find(dm1.Seq());
        map<int, int>::iterator it2;
        if(it1 != checkedSeqs.end()) { 
          it2 = it1->second.find(dm2.Seq()); 
          if(it2 != it1->second.end() && it2->second<=offset) { continue; } //Check if current sequence and offset have already been checked 
        }
        FILE_LOG(logDEBUG3) << "Checking dmer match: dmer1 - " << dm1.ToString() << " dmer2 - " << dm2.ToString() << endl;
        checkedSeqs[dm1.Seq()][dm2.Seq()] = offset; 
        if(dm1.IsMatch(dm2, deviations, acceptSelf)) {
          // Refinement check
          FILE_LOG(logDEBUG3) << "verifying match" << endl;
          MatchInfo matchInfo;
          ValidateMatch(dm1, dm2, matchInfo); 
          HandleMatch(dm1, dm2, matchInfo);
          matchCount++;
          FILE_LOG(logDEBUG3) << "Matched: " << RSToString(dm1.Seq(), 0) << endl << RSToString(dm2.Seq(), 0);
        }
      }
    } 
  }
  return matchCount;
}

void RestSiteMapCore::ValidateMatch(const Dmer& dmer1, const Dmer& dmer2, MatchInfo& matchInfo) const {
  DPMatcher validator;
  float matchScore = validator.FindMatch(dmer1, dmer2, Reads(), matchInfo);
}

void RestSiteMapCore::HandleMatch(const Dmer& dm1, const Dmer& dm2, const MatchInfo& matchInfo) const {
  int offsetBase1 = GetBasePos(dm1.Seq(), dm1.Pos(), false);
  int offsetBase2 = GetBasePos(dm2.Seq(), dm2.Pos(), false);
  int offset = offsetBase1 - offsetBase2;
  bool dir   = (offset>=0)?true:false;
  char dirSign = dir?'+':'-';
  int startBase1 = GetBasePos(dm1, matchInfo.GetFirstMatchPos1(), false); 
  int endBase1   = GetBasePos(dm1, matchInfo.GetLastMatchPos1(), true); 
  int startBase2 = GetBasePos(dm2, matchInfo.GetFirstMatchPos2(), false); 
  int endBase2   = GetBasePos(dm2, matchInfo.GetLastMatchPos2(), true); 
  cout << dm1.Seq() << " " << dm2.Seq() << " " << startBase1 << " " << endBase1 << " "
       << startBase2 << " "  << endBase2 << " "<< dirSign << endl;
}

int RestSiteMapCore::GetBasePos(int seqIdx, int rsPos, bool inclusive) const {
  const RSiteRead& rSites = GetRead(seqIdx);
  int cmPos = rSites.PreDist(); //cumulative position
  int upto = (inclusive? rsPos: rsPos-1);
  for(int i=0; i<=upto; i++) {
    cmPos += rSites[i];
  }
  return cmPos;
}

int RestSiteMapCore::GetBasePos(const Dmer& dm, int rsPos, bool inclusive) const {
  return GetBasePos(dm.Seq(), rsPos, inclusive);
}

string RestSiteMapCore::RSToString(int rIdx, int offset) const {
  return m_rReads[rIdx].ToString(offset);
} 

string RestSiteMapCore::RSToString(const Dmer& dmer) const {
  return RSToString(dmer.Seq(), dmer.Pos());
} 
