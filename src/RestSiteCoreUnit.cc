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
        mm.push_back(i+motifLen/2-n);
      } 
      n = i + motifLen/2; // Set point to middle of motif so that sequences are reversible
    }
    if (origString != "") {
      if(wrotePrefix) { 
        rr.PostDist() = origLen - n - 1; // postfix (number of leading bits after last motif location
      } 
    }
  } 
  rr.Dist() = mm;
  int readIdx = reads.AddRead(rr);
  FILE_LOG(logDEBUG3) << "Adding Read: " << readIdx << "  " << rr.Name() << " " << rr.Ori();
  if(addRC) {  
    rr.Flip();
    readIdx = reads.AddRead(rr);
    FILE_LOG(logDEBUG3) << "Adding Read: " << readIdx << "  " << rr.Name() << " " << rr.Ori();
    return 2*mm.isize(); // Return the total number of sites that have been added
  }
  return mm.isize(); // Return the total number of sites that have been added
}

//TODO this is a very rough way of estimating memory and should be improved 
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

int RestSiteMapCore::FindMapInstances(float indelVariance, map<int, map<int,bool>>& checkedSeqs) const {
  int counter       = 0;
  double matchCount = 0;
  int loopLim       = m_dmers.NumCells();
  svec<int> neighbourCells;
  neighbourCells.reserve(pow(2, m_modelParams.DmerLength()));
  svec<int> deviations;
  deviations.resize(m_modelParams.DmerLength());

  for (int iterIndex=0; iterIndex<loopLim; iterIndex++) {
    counter++;
    if (counter % 100000 == 0) {
//      cout << "\rLOG Progress: " << 100*(double)iterIndex/(double)loopLim << "%" << flush;
    }
    if(!m_dmers[iterIndex].empty()) {
      FILE_LOG(logDEBUG2) << "Number of dmers in cell " << iterIndex << " " << m_dmers[iterIndex].isize(); 
      matchCount += HandleMappingInstance(m_dmers[iterIndex], indelVariance, checkedSeqs, neighbourCells, deviations, false);
    }
  }
  return matchCount;
}

int RestSiteMapCore::HandleMappingInstance(const svec<Dmer>& dmers, float indelVariance, map<int, map<int,bool>>& checkedSeqs,
                                           svec<int>& neighbourCells, svec<int>& deviations, bool acceptSameIdx) const {
  int matchCount = 0;
   for(Dmer dm1:dmers) {
    neighbourCells.clear();
    deviations.clear();
    dm1.CalcDeviations(deviations, indelVariance, m_modelParams.CNDFCoef1()); //TODO this does not need to be redone every time!
    int merLoc = m_dmers.MapNToOneDim(dm1.Data());
    m_dmers.FindNeighbourCells(merLoc, dm1, deviations, neighbourCells); 
    for (int nCell:neighbourCells) {
      for (auto dm2:m_dmers[nCell]) {
        if(checkedSeqs[dm1.Seq()][dm2.Seq()]) {// || checkedSeqs[dm2.Seq()][dm1.Seq()]) {
          continue;  //Check if current pair has not been matched already 
        }
        int offset = abs(dm1.Pos() - dm2.Pos());
        FILE_LOG(logDEBUG3) << "Checking dmer match: dmer1 - " << dm1.ToString() << " dmer2 - " << dm2.ToString() << " offset: " << offset << endl;
        if(dm1.IsMatch(dm2, deviations, acceptSameIdx)) {
          // Refinement check
          FILE_LOG(logDEBUG3) << "verifying match" << endl;
          MatchInfo matchInfo;
          float side1Score, side2Score = 0;
          ValidateMatch(dm1, dm2, indelVariance, matchInfo, side1Score, side2Score); 
          bool passed = WriteMatchPAF(dm1, dm2, matchInfo, side1Score, side2Score);
          if(passed) {
            checkedSeqs[dm1.Seq()][dm2.Seq()]=true;
            matchCount++;
            FILE_LOG(logDEBUG3) << "Matched: " << RSToString(dm1.Seq(), 0) << endl << RSToString(dm2.Seq(), 0);
          }
        }
      }
    } 
  }
  return matchCount;
}

void RestSiteMapCore::ValidateMatch(const Dmer& dmer1, const Dmer& dmer2, float indelVariance, MatchInfo& matchInfo,
                                    float& side1Score, float& side2Score) const {
  DPMatcher validator;
  float matchScore = validator.FindMatch(dmer1, dmer2, Reads(), indelVariance, m_modelParams.CNDFCoef2(), matchInfo, side1Score, side2Score);
}

bool RestSiteMapCore::WriteMatchPAF(const Dmer& dm1, const Dmer& dm2, const MatchInfo& matchInfo, 
                                    float& side1Score, float& side2Score) const {
  string name_query    = GetRead(dm2.Seq()).Name();
  int length_query     = GetBasePos(dm2, GetRead(dm2.Seq()).Size(), true); //This function will find the total length of the sequence in bases
  int startBase_query  = GetBasePos(dm2, matchInfo.GetFirstMatchPos2(), false); 
  int endBase_query    = GetBasePos(dm2, matchInfo.GetLastMatchPos2(), true); 
  char strand_query    = (GetRead(dm2.Seq()).Ori()>0? '+': '-');
  // Items useful for assembly
  int preDist_query    = GetRead(dm2.Seq()).PreDist();
  int postDist_query   = length_query - GetRead(dm2.Seq()).PostDist();
 
  string name_target   = GetRead(dm1.Seq()).Name();
  int length_target    = GetBasePos(dm1, GetRead(dm1.Seq()).Size(), true); //This function will find the total length of the sequence in bases
  int startBase_target = GetBasePos(dm1, matchInfo.GetFirstMatchPos1(), false); 
  int endBase_target   = GetBasePos(dm1, matchInfo.GetLastMatchPos1(), true); 
  // Items useful for assembly
  int preDist_target   = GetRead(dm1.Seq()).PreDist();
  int postDist_target  = length_target - GetRead(dm1.Seq()).PostDist();
  
  float matchScore     = matchInfo.GetIdentScore();
  //int  numMatches      = matchInfo.GetNumMatches();
  int  alignBlockLen   = max(endBase_query-startBase_query, endBase_target-startBase_target);
  int  mappingQual     = 255;

  char delim = '\t';

  if(matchScore>GetThresholdScore()) {
    cout << name_query << delim << length_query << delim << startBase_query 
        << delim << endBase_query << delim << strand_query << delim << name_target 
        << delim << length_target << delim << startBase_target << delim << endBase_target
        << delim << matchScore << delim << alignBlockLen << delim << mappingQual << delim;
    //Auxillary info:
    cout << "queryPreDist:" << preDist_query << delim << "queryPostDist:" << postDist_query << delim 
         << "targetPreDist:" << preDist_target << delim << "targetPostDist:" << postDist_target;
    cout << endl;
    return true;
  }
  return false; //Did not pass acceptance threshold
}

float RestSiteMapCore::GetThresholdScore() const { 
  float thresh = m_modelParams.ScoreThreshold(); 
  if(thresh<0) { //Automatic computation
    thresh = 0.2 + 0.1*(1-exp(-2*(m_modelParams.CNDFCoef2()-1)));  //TODO parameterise and save value not to recalculate everytime
  }
  return thresh;
}

float RestSiteMapCore::GetRandomMatchProb() const {
  float rho       = 1.0/pow(m_modelParams.AlphabetSize(), m_modelParams.MotifLength());
  float theta     = 1 - rho;
  int maxSep      = (1.0/rho)*10.0;
  float alphaPrem = 0;
  for(int sep=0; sep<maxSep; sep++) {
    float thetaSum = 0;
    int tsLim    = m_modelParams.CNDFCoef2()*sqrt(sep*m_dataParams.IndelErr());
    for(int i=0; i<=tsLim; i++) {
      thetaSum += pow(theta, i);
    }
    float toAdd = rho*pow(theta, sep) * rho*pow(theta, sep) * thetaSum;
    alphaPrem  += toAdd; 
  }
  float alpha = 2*alphaPrem - rho/2;
  return alpha;
}

int RestSiteMapCore::GetBasePos(int seqIdx, int rsPos, bool inclusive) const {
  const RSiteRead& rSites = GetRead(seqIdx);
  int cmPos = rSites.PreDist(); //cumulative position
  int upto = (inclusive? rsPos: rsPos-1);
  bool includePostDist = false;
  if(rsPos==rSites.Size()) { 
    upto--;
    includePostDist = true;
  } 
  for(int i=0; i<=upto; i++) {
    cmPos += rSites[i];
  }
  if(includePostDist) { cmPos += rSites.PostDist(); }
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
