#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <cmath>
#include <algorithm>
#include "RestSiteAlignUnit.h"

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
    cout<< "Motif: " << motif << endl;
    m_rsaCores[m_motifs[motifIdx]].BuildDmers();
  }
}

void RestSiteMapper::FindMatches(const string& fileNameQuery, const string& fileNameTarget) {
  GenerateMotifs(); 
  map<int, map<int,int>> checkedSeqs;  // Flagset for sequences that have been searched for a given sequence index and from a specific offset
  int matchCount = 0;
  SetTargetSites(fileNameTarget, !m_modelParams.IsSingleStrand());
  for(int motifIdx=0; motifIdx<m_modelParams.NumOfMotifs(); motifIdx++) {
    string motif = m_motifs[motifIdx];
    FILE_LOG(logDEBUG1) << "Finding matches based on motif: " << motif;
    matchCount += m_rsaCores[motif].FindMapInstances(0.08, checkedSeqs); //TODO parameterise data params
  }
  cout << "Total number of matches recorded: " << matchCount << endl;
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

void RestSiteDBMapper::FindMatches(const string& fileNameQuery, const string& fileNameTarget) {
  GenerateMotifs(); 
  SetTargetSites(fileNameTarget, !m_modelParams.IsSingleStrand());
  SetQuerySites(fileNameQuery);  
  map<int, map<int,int>> checkedSeqs;  // Flagset for sequences that have been searched for a given sequence index and from a specific offset
  int matchCount = 0;
  for(int motifIdx=0; motifIdx<m_modelParams.NumOfMotifs(); motifIdx++) {
    string motif = m_motifs[motifIdx];
    int readCnt   = m_queryReads[motif].NumReads();
    int reportCnt = readCnt/1000;
    if(reportCnt == 0) { reportCnt = 1; }
    for(int rIdx=0; rIdx<readCnt; rIdx++) {
      FILE_LOG(logDEBUG3) << "Finding dmer match candidates for read: " << rIdx; 
      matchCount += m_rsaCores[motif].FindSingleReadMapInstances(m_queryReads[motif][rIdx], rIdx, 0.08, checkedSeqs);
      if (rIdx % reportCnt== 0) {
//        cout << "\rLOG Progress: " << 100*(double)rIdx/(double)readCnt << "%" << flush;
      }
    }
  }
  cout << endl;
  cout << "Total number of matches recorded: " << matchCount << endl;
  //WriteMatchCandids(matchCandids.GetAllCandids());
}

string RestSiteDBMapper::GetQueryName(int readIdx) const {
  if(m_queryReads.empty())  { return ""; }
  else { return m_queryReads.begin()->second[readIdx].Name(); }
}

/*
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
*/

