#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "ryggrad/src/base/Logger.h"
#include "DPMatcher.h"
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
        rr.PostDist() = origLen - n - 1; // postfix (number of leading bits after last motif location, last item in dmer sequence)
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
  if(pow(dimCount, m_modelParams.DmerLength()) > 1900000000) {              //TODO parameterise
    FILE_LOG(logWARNING) << "Input data size is too large";
    dimCount = pow(1900000000, 1.0/m_modelParams.DmerLength()); 
  }
  FILE_LOG(logINFO) << "Estimated number of Dmers and dimension size for dmer storage: " << TotalSiteCount() << "  " << dimCount; 
  m_dmers.BuildDmers(m_rReads , m_modelParams.DmerLength(), m_modelParams.MotifLength(), dimCount); 
}

void RestSiteMapCore::FindMapInstances(float indelVariance) const {
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
      cout << "\rLOG Progress: " << 100*(double)iterIndex/(double)loopLim << "%" << flush;
    }
    if(!m_dmers[iterIndex].empty()) {
      FILE_LOG(logDEBUG2) << "Number of dmers in cell " << iterIndex << " " << m_dmers[iterIndex].isize(); 
      for (Dmer dm1:m_dmers[iterIndex]) {
        neighbourCells.clear();
        dm1.CalcDeviations(deviations, indelVariance, m_modelParams.CNDFCoef()); //TODO this does not need to be redone every time!
        m_dmers.FindNeighbourCells(iterIndex, dm1, deviations, neighbourCells); 
        for (int nCell:neighbourCells) {
          FILE_LOG(logDEBUG3) << "Checking Neighbour cell: " << nCell; 
          for (Dmer dm2:m_dmers[nCell]) {
            FILE_LOG(logDEBUG4) << endl << "Checking dmer match: dmer1 - " << dm1.ToString() << " dmer2 - " << dm2.ToString();
            if(dm1.IsMatch(dm2, deviations, false)) {
              // Refinement check
              if(ValidateMatch(dm1, dm2)) { 
                cout << "Dmer 1: " << dm1.ToString() << endl;
                cout << "Dmer 2: " << dm2.ToString() << endl;
                FILE_LOG(logDEBUG4) << "Found Dmer Match";
                matchCount++;
              }
            }
          }
        }
      }
    }
  }
  cout << "Total number of matches recorded: " << matchCount << endl;
}

void RestSiteMapCore::FindSingleReadMapInstances(const RSiteRead& read, int rIdx, float indelVariance) const {
  svec<int> neighbourCells;
  neighbourCells.reserve(pow(2, m_modelParams.DmerLength()));
  svec<int> deviations;
  deviations.resize(m_modelParams.DmerLength());
  svec<Dmer> dmers;
  m_dmers.GenerateDmers(read, rIdx, dmers);
  for(Dmer dm1:dmers) {
    dm1.CalcDeviations(deviations, indelVariance, m_modelParams.CNDFCoef()); //TODO this does not need to be redone every time!
    int merLoc = m_dmers.MapNToOneDim(dm1.Data());
    neighbourCells.clear();
    m_dmers.FindNeighbourCells(merLoc, dm1, deviations, neighbourCells); 
    for (int nCell:neighbourCells) {
      for (auto dm2:m_dmers[nCell]) {
        FILE_LOG(logDEBUG4) << endl << "Checking dmer match: dmer1 - " << dm1.ToString() << " dmer2 - " << dm2.ToString();
        if(dm1.IsMatch(dm2, deviations, true)) {
          // Refinement check
          if(ValidateMatch(dm1, dm2)) { 
            cout << "Dmer 1: " << dm1.ToString() << endl;
            cout << "Dmer 2: " << dm2.ToString() << endl;
          }
        }
      }
    } 
  }
}
bool RestSiteMapCore::ValidateMatch(const Dmer& dmer1, const Dmer& dmer2) const {
  DPMatcher validator;
  float matchScore = validator.FindMatchScore(MappedInstance(dmer1, dmer2), Reads());
  cout << matchScore << endl;
  if(matchScore>0.9) { return true; }
  else { return false; }
}

string RestSiteMapCore::RSToString(int rIdx, int offset) {
  return m_rReads[rIdx].ToString(offset);
} 

string RestSiteMapCore::RSToString(const Dmer& dmer) {
  return RSToString(dmer.Seq(), dmer.Pos());
} 