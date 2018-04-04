#ifndef DMER_H
#define DMER_H

#include <map>
#include <string>
#include "RSiteReads.h"

class Dmer {
public:
  Dmer() {
    m_seq = -1;
    m_pos = -1;
  }

  bool operator <  (const Dmer & m) const;
  bool operator != (const Dmer & m) const; 
  bool operator == (const Dmer & m) const; 
  int operator[](int idx) const { return m_data[idx];     }
  int& operator[](int idx)      { return m_data[idx];     }

  int & Seq() {return m_seq;}
  int & Pos() {return m_pos;}
 
  const int & Seq() const {return m_seq;}
  const int & Pos() const {return m_pos;}
 
  svec<int> & Data() {return m_data;}
  const svec<int> & Data() const {return m_data;}

  inline bool IsMatch(const Dmer& otherDmer, const svec<int>& deviations, bool allowSame) const {
    if(!allowSame && m_seq == otherDmer.m_seq) { return false; } // Same sequence is not a real match
    for(int i=0; i<m_data.isize(); i++) {
      if (m_data[i] < otherDmer.m_data[i]-deviations[i] || m_data[i] > otherDmer.m_data[i]+deviations[i])
        return false;
    }
    return true;
  }
  void CalcDeviations(svec<int>& deviations, float indelVariance, float deviationCoeff) const; 
  string ToString() const; 

private:
  svec<int> m_data;
  int m_seq;
  int m_pos;
};

class Dmers {
public:
  Dmers(): m_mers(), m_dimCount(0), m_dmerLength(0), m_dimRangeBounds(), m_dmerCellMap(), m_dmerCount(0) {}

  int NumMers() const                      { return m_dmerCount;     }
  int NumCells() const                     { return m_mers.isize();  }
  svec<Dmer> operator[](int index) const   { return m_mers[index];   }
  svec<Dmer>& operator[](int index)        { return m_mers[index];   }

  void BuildDmers(const RSiteReads& rReads, int dmerLength, int motifLength, int countPerDimension); 
  void FindNeighbourCells(int initVal, const Dmer& dmer, const svec<int>& deviations, svec<int>& result) const; 
  void GenerateDmers(const RSiteRead& rRead, int rIdx, svec<Dmer>& dmers) const;
  int MapNToOneDim(const svec<int>& nDims) const;
  svec<int> MapOneToNDim(int oneDMappedVal) const;

protected:
  void SetRangeBounds(int motifLength);
  void AddSingleReadDmers(const RSiteRead& rRead, int rIdx);
  void FindNeighbourCells(int initVal, const Dmer& dmer, const svec<int>& deviations, int depth, svec<int>& result) const; 

private:
  svec<svec<Dmer> > m_mers;    /// Multi-dimensional matrix representation of dmers projected on to dimensions
  int m_dimCount;              /// Number of cells in each dimension (this is dependent on the site values and the reduction coefficient)
  int m_dmerLength;            /// Number of dimensions in the matrix (i.e. dmer length)
  svec<int> m_dimRangeBounds;  /// The range limits for dmer values to be placed in each dimennsion
  map<int, int> m_dmerCellMap; /// Mapping every dmer value to the relevant cell placement
  int m_dmerCount;             /// Total number of dmers
};

#endif //DMER_H

