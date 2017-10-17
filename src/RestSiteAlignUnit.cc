#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <cmath>
#include <algorithm>
#include <sstream>
#include "RestSiteAlignUnit.h"

string RSiteRead::ToString() const {
  stringstream ss;
  ss << PreDist() << ", ";
  for (int i=0; i<m_dist.isize(); i++) {
    ss << m_dist[i] << ", ";
  }
  ss << PostDist();
  return ss.str();
}

void RSiteRead::Flip() {
  m_ori = -m_ori;
  svec<int> tmp;
  tmp.resize(m_dist.isize());
  int k = m_dist.isize()-1;
  for (int i=0; i<m_dist.isize(); i++) {
    tmp[k] = m_dist[i];
    k--;
  }
  m_dist = tmp;
  int tmp_pp; //Swap pre/postfix
  tmp_pp     = m_preDist;
  m_preDist  = m_postDist;
  m_postDist = tmp_pp;
}

int RSiteReads::AddRead(const RSiteRead& rr) {
  m_rReads.push_back(rr);
  m_readCount++;
  return m_readCount-1;
}

string RSiteReads::ToString() const {
  string strOut;
  for(int i=0; i<m_readCount; i++) {
    strOut += m_rReads[i].ToString();
    strOut += "\n";
  }
  return strOut;
}

bool Dmer::operator < (const Dmer & m) const {
  for (int i=0; i<m_data.isize(); i++) {
    if (m_data[i] != m.m_data[i])
      return m_data[i] < m.m_data[i];
  }
  return false;
}
  
bool Dmer::operator != (const Dmer & m) const {
  for (int i=0; i<m_data.isize(); i++) {
    if (m_data[i] != m.m_data[i])
      return true;
  }
  return false;
}

bool Dmer::operator == (const Dmer & m) const {
  for (int i=0; i<m_data.isize(); i++) {
    if (m_data[i] != m.m_data[i])
      return false;
  }
  return true;
}

void Dmer::CalcDeviations(svec<int>& deviations, float indelVariance, float deviationCoeff) const {
  if(deviations.isize() < m_data.isize()) { deviations.resize(m_data.isize()); }
  for(int i=0; i<m_data.isize(); i++) {
    int deviation = sqrt(m_data[i]*indelVariance)*deviationCoeff;
    deviations[i] = deviation;
  }
}
 

string Dmer::ToString() const {
  stringstream ss;
  for (int i=0; i<m_data.isize(); i++)
    ss << " " << m_data[i];
  ss << " seq: " << m_seq << " pos: " << m_pos;
  return ss.str();
}

void Dmers::BuildDmers(const RSiteReads& rReads , int dmerLength, int motifLength, int countPerDimension) { 
  m_dmerLength = dmerLength;
  m_dimCount   = countPerDimension;
  m_mers.resize(pow(m_dimCount, m_dmerLength)); // TODO Check to be within memory limit
  SetRangeBounds(motifLength);
  cout << "Building dmers ..." << endl;
  FILE_LOG(logINFO) << "LOG Build mer list...";
  for (int rIdx=0; rIdx<rReads.NumReads(); rIdx++) {
    AddSingleReadDmers(rReads[rIdx], rIdx);
  }
  cout << "Total number of dmers: " << NumMers() << endl;
  FILE_LOG(logINFO) << "Total number of dmers: " << NumMers();
}

void Dmers::SetRangeBounds(int motifSize) {
  double p1     = 1.0 / pow(4, motifSize); // for example for a motif size of 4 this will be 1/256
  double p2     = 1.0 - p1;                // for example for a motif size of 4 this will be 255/256
  double pC     = 0;                       // Cumulative probability
  int rangeLim  = 0;

  for(int dim=0; dim<m_dimCount-1; dim++) {
    while(pC < (double)(dim+1)/m_dimCount) {
      double pi = pow(p2, rangeLim+1) * p1;
      pC += pi;
      rangeLim++;
      m_dmerCellMap[rangeLim] = dim;
    }
    m_dimRangeBounds.push_back(rangeLim); 
    FILE_LOG(logDEBUG1) << "Dimension Range: " << m_dimRangeBounds.isize()-1 << "  " << rangeLim; 
  }
}

void Dmers::AddSingleReadDmers(const RSiteRead& rRead, int rIdx) {
  svec<Dmer> dmers;
  GetDmers(rRead, rIdx, dmers);
  for (Dmer dm:dmers) {
    int merLoc = MapNToOneDim(dm.Data());
    m_mers[merLoc].push_back(dm);
    m_dmerCount++;
  }
  FILE_LOG(logDEBUG3) << "Read Index: " << rIdx << " total dmers so far: " << m_dmerCount << endl;
}

void Dmers::GetDmers(const RSiteRead& rRead, int rIdx, svec<Dmer>& dmers) const {
  Dmer mm;
  mm.Seq() = rIdx;
  mm.Data().resize(m_dmerLength);
  int loopLim = rRead.Dist().isize() - m_dmerLength;
  for (int i=0; i<=loopLim; i++) {
    mm.Pos() = i;
    for (int j=0; j<m_dmerLength; j++) {
      mm.Data()[j] = rRead.Dist()[i+j];
    }
    dmers.push_back(mm);
  }
}

int Dmers::MapNToOneDim(const svec<int>& nDims) const {
  // This function does not do bound checking and assumes that nDims size is m_dmerLength and values are between 0 and m_dimCount
  int mapVal = 0;
  int coeff  = pow(m_dimCount, m_dmerLength-1);
  for(int i=0; i<m_dmerLength; i++) {
    int qVal = m_dimCount - 1; // First set it to the highest possible digit and then check if it belongs in another cell 
    if(nDims[i] < m_dimRangeBounds[m_dimCount-2])  { qVal = m_dmerCellMap.at(nDims[i]); } // the last digit range bound is in m_dimCount-2
    mapVal += qVal * coeff;
    coeff  /= m_dimCount;
  }
  return mapVal;
}

svec<int> Dmers::MapOneToNDim(int oneDMappedVal) const {
  svec<int> nDims;
  nDims.resize(m_dmerLength);
  int coeff  = m_dimCount;
  for(int i=m_dmerLength-1; i>=0; i--) {
    nDims[i] = oneDMappedVal % coeff;
    oneDMappedVal /= coeff;
  }
  return nDims;
}

void Dmers::FindNeighbourCells(int initVal, const Dmer& dmer, const svec<int>& deviations, svec<int>& result) const {
  FindNeighbourCells(initVal, dmer, deviations, m_dmerLength-1, result);
}

void Dmers::FindNeighbourCells(int initVal, const Dmer& dmer, const svec<int>& deviations, int depth, svec<int>& result) const {
  if(depth == -1) { 
    result.push_back(initVal);
    return;
  }
  FindNeighbourCells(initVal, dmer, deviations, depth-1, result);
  svec<int> tempResult = result;
  int currDigit = m_dimCount - 1; // First set it to the highest possible digit and then check if it belongs in another cell 
  if(dmer[depth] < m_dimRangeBounds[m_dimCount-2])  { currDigit = m_dmerCellMap.at(dmer[depth]); } // the last digit range bound is in m_dimCount-2
  if((currDigit < m_dimCount-2)  //only add one to the current digit if it has room to be increased 
    && (dmer[depth]+deviations[depth] > m_dimRangeBounds[currDigit])) { // only try one cell up if the deviation limits don't fall within the same cell
    for(int elem:tempResult) {
      int newElem = elem + pow(m_dimCount, (m_dmerLength-1-depth));
      result.push_back(newElem);
    }
  }
}

string MatchCandid::ToString() const {
  stringstream ss;
  ss << GetFirstReadIndex() << "\t" <<  GetFirstMatchPos()  << "\t" 
     << GetSecondReadIndex() << "\t" << GetSecondMatchPos(); 
  return ss.str();
}

void MatchCandids::AddCandidSort(int rIdx1, int rIdx2, int rPos1, int rPos2) {
/*
  // Make sure that rIdx1, rIdx2 are in increasing order (so that sorting will bring all relevant pairs together)
  if(rIdx1>rIdx2) {
    int temp = rIdx1;
    rIdx1 = rIdx2;
    rIdx2 = temp; // Swap rIdx1 & rIdx2
    temp  = rPos1;
    rPos1 = rPos2;
    rPos2 = temp;
  }
  m_candids.push_back(MatchCandid(rIdx1, rIdx2, rPos1, rPos2));
*/
  m_candids[rIdx1][rIdx2]++;
  FILE_LOG(logDEBUG3) << "Adding overlap candidate: " << MatchCandid(rIdx1, rIdx2, rPos1, rPos2).ToString();
}
 
string MatchCandids::ToString() const {
  stringstream ss;
  for(auto const &ent1 : m_candids) {
    for(auto const &ent2 : ent1.second) {
      ss << ent1.first << " " << ent2.first << " " <<  ent2.second << endl;
    }
  }
  return ss.str();
}

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

void RestSiteMapCore::FindLapCandids(float indelVariance, MatchCandids& matchCandids) const {
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
              FILE_LOG(logDEBUG4) << "Found Dmer Match";
              matchCount++;
              matchCandids.AddCandidSort(dm1.Seq(), dm2.Seq(), dm1.Pos(), dm2.Pos());
            }
          }
        }
      }
    }
  }
  cout << "Total number of matches recorded: " << matchCount << endl;
}

void RestSiteMapCore::FindSingleReadMatchCandids(const RSiteRead& read, int rIdx, float indelVariance, MatchCandids& matchCandids) const {
  svec<int> neighbourCells;
  neighbourCells.reserve(pow(2, m_modelParams.DmerLength()));
  svec<int> deviations;
  deviations.resize(m_modelParams.DmerLength());
  svec<Dmer> dmers;
  m_dmers.GetDmers(read, rIdx, dmers);
  for(Dmer dm1:dmers) {
    dm1.CalcDeviations(deviations, indelVariance, m_modelParams.CNDFCoef()); //TODO this does not need to be redone every time!
    int merLoc = m_dmers.MapNToOneDim(dm1.Data());
    neighbourCells.clear();
    m_dmers.FindNeighbourCells(merLoc, dm1, deviations, neighbourCells); 
    for (int nCell:neighbourCells) {
      for (auto dm2:m_dmers[nCell]) {
        FILE_LOG(logDEBUG4) << endl << "Checking dmer match: dmer1 - " << dm1.ToString() << " dmer2 - " << dm2.ToString();
        if(dm1.IsMatch(dm2, deviations, true)) {
          matchCandids.AddCandidSort(dm1.Seq(), dm2.Seq(), dm1.Pos(), dm2.Pos());
        }
      }
    } 
  }
}

/*
void RestSiteMapCore::FinalOverlaps(const MatchCandids& matchCandids, int tolerance,  MatchCandids& finalMatches) {
  cout << "LOG Refine overlap candidates... " << matchCandids.NumCandids() << endl;
  finalMatches.ReserveInit(matchCandids.NumCandids()/2); // Rough estimate 
  MatchCandid currCandid;
  int rejectCnt = 0;
  for(int i=0; i<matchCandids.NumCandids(); i++) {
    if(matchCandids[i]==currCandid) { continue; }
    currCandid = matchCandids[i];  
    const RSiteRead& read1    = m_rReads[currCandid.GetFirstReadIndex()]; 
    const RSiteRead& read2    = m_rReads[currCandid.GetSecondReadIndex()]; 
    int delta     = currCandid.GetOffsetDelta();
    int readPos1  = delta>=0 ? 0 : -delta;
    int readPos2  = delta>=0 ? delta : 0;
    bool overlaps = true;
    while(readPos1 < read1.Size() && readPos2 < read2.Size()) {
      if(abs(read1[readPos1] - read2[readPos2]) <= tolerance) {
        readPos1++;
        readPos2++;
      } else {
        overlaps = false;
        rejectCnt++;
        break;
      }
    }
    if(overlaps) { finalMatches.AddCandid(currCandid); }
  }
  cout << "LOG Total number of overlap candidates rejected: " << rejectCnt << endl;
  cout << "LOG Final number of overlaps: " << finalMatches.NumCandids() << endl;
  WriteLapCandids(finalMatches);
}
*/

void RestSiteGeneral::GenerateMotifs() {
  m_motifs.reserve(m_modelParams.NumOfMotifs());
  vector<char> alphabet = {'A', 'C', 'G', 'T'}; //Should be in lexographic order
  map<char, char> RCs = {{'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'}}; 
  vector<vector<char>> tempMotifs; 
  CartesianPower(alphabet, m_modelParams.MotifLength(), tempMotifs);
  for(vector<char> sElem:tempMotifs){
    string motif = "";
    for(char cElem:sElem) {
      motif += cElem;
    }
    if(motif.length() == m_modelParams.MotifLength()) { 
      if(ValidateMotif(motif, alphabet, RCs)) {
        FILE_LOG(logDEBUG1) << "Motif: "  << m_motifs.isize() << " " << motif;
        m_motifs.push_back(motif); 
        if(m_motifs.isize() == m_modelParams.NumOfMotifs()) {
          break;
        }
      }
    } 
  }
  if(m_motifs.isize() < m_modelParams.NumOfMotifs()) {
    FILE_LOG(logWARNING) << "Could not generate the number of requested motifs - maximum " << m_motifs.isize() << " being used";
    m_modelParams.ChangeNumOfMotifs(m_motifs.isize());
  }
}
  
bool RestSiteGeneral::ValidateMotif(const string& motif, const vector<char>& alphabet, const map<char, char>& RCs) const {
  int motifLen = motif.length();

  // 1. Simplicity Filter
  map<char, int> alphabetCnt;
  bool simple = false;
  for(int i=0; i<motifLen-1; i++) {
    alphabetCnt[motif[i]]++;
    if(alphabetCnt[motif[i]]>=motifLen/2 ||
      (i<motifLen-1 && motif[i]==motif[i+1])) {
      simple = true; 
      break;
    }
  }
  if(simple) { return false; }

  // 2. RC Filter (Palindrome)
  if(motifLen%2 != 0) { return false; }
  for(int ii=0; ii<motifLen/2; ii++) {
    if(RCs.at(motif[ii]) != motif[motifLen-ii-1]) { return false; } 
  }

  return true;
}

void RestSiteGeneral::CartesianPower(const vector<char>& input, unsigned k, vector<vector<char>>& result) const {
  if (k == 1) {
    for (char value: input) {
      result.push_back( {value} );
    }
    return;
  } else {
    CartesianPower(input, k - 1, result);
    vector<vector<char>> smallerPower = result;
    for (int elem: input) {
      for (vector<char> sublist: smallerPower) {
        sublist.push_back(elem);
        result.push_back(sublist);
      }
    }
    return;
  }
} 

void RestSiteGeneral::WriteMatchCandids(const map<int, map<int, int> >& candids) const {
  for(auto const &ent1 : candids) {
    for(auto const &ent2 : ent1.second) {
      cout << ent1.first << " " << ent2.first << " " <<  ent2.second << endl;
    }
  }
}

string RestSiteGeneral::GetTargetName(int readIdx) const {
  if(m_rsaCores.empty())  { return ""; }
  else { return m_rsaCores.begin()->second.Reads()[readIdx].Name(); }
}

void RestSiteGeneral::SetTargetSites(const string& fileName, bool addRC) {
  for(int motifIdx=0; motifIdx<m_modelParams.NumOfMotifs(); motifIdx++) {
    string motif = m_motifs[motifIdx];
    m_rsaCores[motif] = RestSiteMapCore(motif, m_modelParams, m_dataParams);
  }
 
  FlatFileParser parser;
  parser.Open(fileName);
  string l;
  string name;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.Line()[0] == '>') {
      transform(l.begin(), l.end(), l.begin(), ::toupper);
      for(int motifIdx=0; motifIdx<m_modelParams.NumOfMotifs(); motifIdx++) {
        string motif = m_motifs[motifIdx];
        int totSiteCnt = m_rsaCores[motif].CreateRSitesPerString(l, name, m_rsaCores[motif].Reads(), addRC);
        m_rsaCores[motif].IncTotalSiteCount(totSiteCnt);
      }
      l.clear();
      name = parser.Line();
    }
    l += parser.Line();
  }
  if( l != "") {
    transform(l.begin(), l.end(), l.begin(), ::toupper);
    for(int motifIdx=0; motifIdx<m_modelParams.NumOfMotifs(); motifIdx++) {
      string motif = m_motifs[motifIdx];
      int totSiteCnt = m_rsaCores[motif].CreateRSitesPerString(l, name, m_rsaCores[motif].Reads(), addRC);
      m_rsaCores[motif].IncTotalSiteCount(totSiteCnt);
    }
  }
  for(int motifIdx=0; motifIdx<m_modelParams.NumOfMotifs(); motifIdx++) {
    string motif = m_motifs[motifIdx];
    m_rsaCores[m_motifs[motifIdx]].BuildDmers();
  }
}

void RestSiteMapper::FindMatches(const string& fileNameQuery, const string& fileNameTarget, MatchCandids& finalMatches) {
  GenerateMotifs(); 
  SetTargetSites(fileNameTarget, true);
  MatchCandids matchCandids;
  for(int motifIdx=0; motifIdx<m_modelParams.NumOfMotifs(); motifIdx++) {
    string motif = m_motifs[motifIdx];
    FILE_LOG(logDEBUG1) << "Finding matches based on motif: " << motif;
    m_rsaCores[motif].FindLapCandids(0.08, matchCandids); //TODO parameterise data params
  }
}

void RestSiteDBMapper::SetQuerySites(const string& fileName) {
  int totSiteCnt = 0;
  FlatFileParser parser;
  parser.Open(fileName);
  string l;
  string name;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.Line()[0] == '>') {
      transform(l.begin(), l.end(), l.begin(), ::toupper);
      for(int motifIdx=0; motifIdx<m_modelParams.NumOfMotifs(); motifIdx++) {
        string motif = m_motifs[motifIdx];
        totSiteCnt += m_rsaCores[motif].CreateRSitesPerString(l, name, m_queryReads[motif], false);
      }
      l.clear();
      name = parser.Line();
    }
    l += parser.Line();
  }
  if( l != "") {
    transform(l.begin(), l.end(), l.begin(), ::toupper);
    for(int motifIdx=0; motifIdx<m_modelParams.NumOfMotifs(); motifIdx++) {
      string motif = m_motifs[motifIdx];
      totSiteCnt += m_rsaCores[motif].CreateRSitesPerString(l, name, m_queryReads[motif], false);
    }
  }
}

void RestSiteDBMapper::FindMatches(const string& fileNameQuery, const string& fileNameTarget, MatchCandids& finalMatches) {
  GenerateMotifs(); 
  SetTargetSites(fileNameTarget, true);
  SetQuerySites(fileNameQuery);  
  MatchCandids matchCandids;
  for(int motifIdx=0; motifIdx<m_modelParams.NumOfMotifs(); motifIdx++) {
    string motif = m_motifs[motifIdx];
    int readCnt   = m_queryReads[motif].NumReads();
    int reportCnt = readCnt/1000;
    if(reportCnt == 0) { reportCnt = 1; }
    for(int rIdx=0; rIdx<readCnt; rIdx++) {
      FILE_LOG(logDEBUG3) << "Finding dmer match candidates for read: " << rIdx; 
      m_rsaCores[motif].FindSingleReadMatchCandids(m_queryReads[motif][rIdx], rIdx, 0.08, matchCandids);
      if (rIdx % reportCnt== 0) {
        cout << "\rLOG Progress: " << 100*(double)rIdx/(double)readCnt << "%" << flush;
      }
    }
  }
  cout << endl;
  WriteMatchCandids(matchCandids.GetAllCandids());
}

string RestSiteDBMapper::GetQueryName(int readIdx) const {
  if(m_queryReads.empty())  { return ""; }
  else { return m_queryReads.begin()->second[readIdx].Name(); }
}


void RestSiteDBMapper::WriteMatchCandids(const map<int, map<int, int> >& candids) const {
  for(auto const &ent1 : candids) {
    int savedMaxScore = 0;
    int savedIdx      = 0;
    for(auto const &ent2 : ent1.second) {
      if(ent2.second > savedMaxScore) { 
        savedMaxScore = ent2.second;
        savedIdx      = ent2.first;
      }
    }
    cout << GetQueryName(ent1.first) << " " << GetTargetName(savedIdx) << " " <<  savedMaxScore << endl;
  }
}

