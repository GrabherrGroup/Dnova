#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "ryggrad/src/base/Logger.h"
#include "Dmers.h"

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
  double p1     = 1.0 / pow(4, motifSize); // for example for a motif size of 4 this will be 1/256 //TODO parameterise 4
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
  GenerateDmers(rRead, rIdx, dmers);
  for (Dmer dm:dmers) {
    int merLoc = MapNToOneDim(dm.Data());
    m_mers[merLoc].push_back(dm);
    m_dmerCount++;
  }
  FILE_LOG(logDEBUG3) << "Read Index: " << rIdx << " total dmers so far: " << m_dmerCount << endl;
}

void Dmers::GenerateDmers(const RSiteRead& rRead, int rIdx, svec<Dmer>& dmers) const {
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

