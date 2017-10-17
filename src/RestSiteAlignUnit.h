#ifndef OPTIMAPALIGNUNIT_H
#define OPTIMAPALIGNUNIT_H

#include <string>
#include <parallel/algorithm>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"

class RSiteRead
{
public:
  RSiteRead() {
    m_ori = 1;
  }
  
  int operator[](int idx) const  { return m_dist[idx];    }
  svec<int> & Dist()             { return m_dist;         }
  string & Name()                { return m_name;         }
  int& PreDist()                 { return m_preDist;      }
  int& PostDist()                { return m_postDist;     }
  int & Ori()                    { return m_ori;          }
  const svec<int> & Dist() const { return m_dist;         }
  const int& PreDist() const     { return m_preDist;      }
  const int& PostDist() const    { return m_postDist;     }
  const string & Name() const    { return m_name;         }
  const int & Ori() const        { return m_ori;          }
  int Size() const               { return m_dist.isize(); }

  string ToString() const;
  void Flip();

private:
  svec<int> m_dist;       /// Distmer values
  int m_preDist;          /// Number of bits prior to the start of the first distmer value
  int m_postDist;         /// Number of bits left over after the last distmer value
  string m_name;          /// Name of optiRead
  int m_ori;              /// Orientation
};

class RSiteReads 
{
public:
  // Default Ctor
  RSiteReads(): m_readCount(0), m_rReads() {}

  void Reserve(int size)                     { m_rReads.reserve(size);   }

  const RSiteRead& operator[](int idx) const { return m_rReads[idx];     }
  RSiteRead& operator[](int idx)             { return m_rReads[idx];     }
  int NumReads() const                       { return m_readCount;       }

  int AddRead(const RSiteRead& rr); 
  string ToString() const;

private:
  int m_readCount;
  svec<RSiteRead> m_rReads;  /// List of site reads
};

class Dmer
{
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
      if (m_data[i] > otherDmer.m_data[i]+deviations[i] || m_data[i] < otherDmer.m_data[i]-deviations[i])
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

class Dmers 
{
public:
  Dmers(): m_mers(), m_dimCount(0), m_dmerLength(0), m_dimRangeBounds(), m_dmerCellMap(), m_dmerCount(0) {}

  int NumMers() const                      { return m_dmerCount;     }
  int NumCells() const                     { return m_mers.isize();  }
  svec<Dmer> operator[](int index) const   { return m_mers[index];   }
  svec<Dmer>& operator[](int index)        { return m_mers[index];   }

  void BuildDmers(const RSiteReads& rReads, int dmerLength, int motifLength, int countPerDimension); 
  inline void FindNeighbourCells(int initVal, const Dmer& dmer, const svec<int>& deviations, svec<int>& result) const; 

  void SetRangeBounds(int motifLength);
  void AddSingleReadDmers(const RSiteRead& rRead, int rIdx);
  void GetDmers(const RSiteRead& rRead, int rIdx, svec<Dmer>& dmers) const;
  inline int MapNToOneDim(const svec<int>& nDims) const;
  inline svec<int> MapOneToNDim(int oneDMappedVal) const;
  inline void FindNeighbourCells(int initVal, const Dmer& dmer, const svec<int>& deviations, int depth, svec<int>& result) const; 

private:
  svec<svec<Dmer> > m_mers;    /// Multi-dimensional matrix representation of dmers projected on to dimensions
  int m_dimCount;              /// Number of cells in each dimension (this is dependent on the site values and the reduction coefficient)
  int m_dmerLength;            /// Number of dimensions in the matrix (i.e. dmer length)
  svec<int> m_dimRangeBounds;  /// The range limits for dmer values to be placed in each dimennsion
  map<int, int> m_dmerCellMap; /// Mapping every dmer value to the relevant cell placement
  int m_dmerCount;             /// Total number of dmers
};

class MatchCandid
{
public:
  MatchCandid(): m_rIdx1(-1), m_rIdx2(-1), m_rPos1(-1), m_rPos2(-1) {}
  MatchCandid(int idx1, int idx2, int rp1, int rp2): m_rIdx1(idx1), m_rIdx2(idx2), m_rPos1(rp1), m_rPos2(rp2) {}

  int GetFirstReadIndex() const  { return m_rIdx1;       }
  int GetSecondReadIndex() const { return m_rIdx2;       }
  int GetFirstMatchPos() const   { return m_rPos1;       }
  int GetSecondMatchPos() const  { return m_rPos2;       }

  inline bool operator < (const MatchCandid& rhs) const {
    return(tie(m_rIdx1, m_rIdx2, m_rPos1)
       < tie(rhs.m_rIdx1, rhs.m_rIdx2, rhs.m_rPos1)); //keep the same order
  }

  inline bool operator == (const MatchCandid& rhs) const {
    return(tie(m_rIdx1, m_rIdx2, m_rPos1)
       == tie(rhs.m_rIdx1, rhs.m_rIdx2, rhs.m_rPos2)); //keep the same order
  }

  string ToString() const;

private:
  int m_rIdx1;        /// The index of the first read
  int m_rIdx2;        /// The index of the second read 
  int m_rPos1;        /// Position of match in first read 
  int m_rPos2;        /// Position of match in second read 
};

class MatchCandids
{
public:
  MatchCandids(): m_candids() {}
  
  //int NumCandids() const { return m_candids.isize(); }

  void AddCandidSort(int rIdx1, int rIdx2, int rPos1, int rPos2);
  
  string ToString() const; 

  const map<int, map<int, int > >& GetAllCandids() const { return m_candids; }
private:
  map<int, map<int, int> > m_candids;
};

class RestSiteDataParams 
{
public:
  RestSiteDataParams( int totalNumReads=10000000, int meanReadLength=10000, int minMapLength=4000, 
                      float deletionErr=0.03, float insertionErr=0.03, float substitutionErr=0.03)
                     :m_totalNumReads(totalNumReads), m_meanReadLength(meanReadLength), m_minMapLength(minMapLength),
                      m_deletionErr(deletionErr), m_insertionErr(insertionErr), m_substitutionErr(substitutionErr) { }

    bool   TotalNumReads() const     { return m_totalNumReads;   }
    int    MeanReadLength() const    { return m_meanReadLength;  }  
    int    MeanMapLength() const     { return m_minMapLength;    } 
    float  DeletionErr() const       { return m_deletionErr;     }
    float  InsertionErr() const      { return m_insertionErr;    }
    float  SubstitutionErr() const   { return m_substitutionErr; }

private: 
  bool    m_totalNumReads;    /// Flag specifying whether the reads are single or double strand
  int     m_meanReadLength;   /// Length of each motif
  int     m_minMapLength;   /// Number of motifs to generate/use
  float   m_deletionErr;      /// The length of distmers to use for seed finding
  float   m_insertionErr;     /// The length of distmers to use for seed finding
  float   m_substitutionErr;  /// The length of distmers to use for seed finding
};

class RestSiteModelParams 
{
public:
  RestSiteModelParams(bool singleStrand=false, int motifLength=4, int numOfMotifs=1, 
                      int dmerLength=6, float cndfCoef=1.0, const vector<char>& alphabet= {'A', 'C', 'G', 'T' }) 
                     :m_singleStrand(singleStrand), m_motifLength(motifLength), m_numOfMotifs(numOfMotifs),
                      m_dmerLength(dmerLength), m_cndfCoef(cndfCoef), m_alphabet(alphabet) { }

  bool   IsSingleStrand() const        { return m_singleStrand;   }
  int    MotifLength() const           { return m_motifLength;    }  
  int    NumOfMotifs() const           { return m_numOfMotifs;    } 
  int    DmerLength() const            { return m_dmerLength;     }
  float  CNDFCoef() const              { return m_cndfCoef;       }
  const vector<char>& Alphabet() const { return m_alphabet;       }

  void ChangeNumOfMotifs(int motifCnt) { m_numOfMotifs = motifCnt; }
private: 
  bool    m_singleStrand;   /// Flag specifying whether the reads are single or double strand
  int     m_motifLength;    /// Length of each motif
  int     m_numOfMotifs;    /// Number of motifs to generate/use
  int     m_dmerLength;     /// The length of distmers to use for seed finding
  float   m_cndfCoef;       /// Cumulative Normal Distribution Function coefficeint used for estimating similarity 
  vector<char>  m_alphabet; /// Alphabet containing base letters used in the reads/motifs in lexographic order 
};

//Forward Declaration
class RestSiteGeneral;

class RestSiteMapCore 
{
  friend RestSiteGeneral;

public:
  //Default Ctor
  RestSiteMapCore(): m_motif(), m_totalSiteCnt(0), m_rReads(), m_dmers() {}

  //Ctor 1
  RestSiteMapCore(string motif, const RestSiteModelParams& mp, const RestSiteDataParams& dp)
                   : m_motif(motif), m_modelParams(mp), m_dataParams(dp), m_rReads(), m_dmers() {}

  int  TotalSiteCount() const              { return m_totalSiteCnt; }
  void IncTotalSiteCount(int cnt)          { m_totalSiteCnt += cnt; }
  const RSiteRead& GetRead(int rIdx) const { return m_rReads[rIdx]; }

  int  CreateRSitesPerString(const string& origString, const string& origName, RSiteReads& reads, bool addRC) const; 

  void BuildDmers(); 
  void FindLapCandids(float indelVariance, MatchCandids& matchCandids) const;
  void FindSingleReadMatchCandids(const RSiteRead& reads, int rIdx, float indelVariance, MatchCandids& matchCandids) const; 
  void FinalOverlaps(const MatchCandids& matchCandids, int tolerance, MatchCandids& finalMatches) const;

protected:
  RSiteReads& Reads()             { return m_rReads; }
  const RSiteReads& Reads() const { return m_rReads; }

private:
  string m_motif;                    /// Vector of all motifs for which restriction site reads have been generated
  RestSiteModelParams m_modelParams; /// Model Parameters
  RestSiteDataParams m_dataParams;   /// Model Parameters
  double  m_totalSiteCnt;            /// The total of restriction site count over all reads
  RSiteReads m_rReads;               /// Restriction Site reads per motif
  Dmers  m_dmers;                    /// To build dmers from restriction site reads
};

class RestSiteGeneral 
{
public:
  RestSiteGeneral(): m_rsaCores(), m_motifs(), m_modelParams(), m_dataParams() {}
  RestSiteGeneral(const RestSiteModelParams& mParams): m_motifs(), m_modelParams(mParams), m_dataParams() {}

  /* Generate Permutation of the given alphabet to reach number of motifs required */
  void GenerateMotifs();  
  bool ValidateMotif(const string& motif, const vector<char>& alphabet, const map<char, char>& RCs) const; 
  void SetTargetSites(const string& fileName, bool addRC); 
  string GetTargetName(int readIdx) const;


  virtual void WriteMatchCandids(const map<int, map<int, int> >& candids) const; 
  virtual void FindMatches(const string& fileNameQuery, const string& fileNameTarget, MatchCandids& finalMatches) = 0; 

protected:
  void CartesianPower(const vector<char>& input, unsigned k, vector<vector<char>>& result) const; 
  map<string, RestSiteMapCore> m_rsaCores;   /// Mapping engine (core data and functionality) per motif
  svec<string> m_motifs;                     /// Vector of all motifs for which restriction site reads have been generated
  RestSiteModelParams m_modelParams;         /// Model Parameters
  RestSiteDataParams m_dataParams;           /// Model Parameters
};

class RestSiteMapper : public RestSiteGeneral 
{
public:
  RestSiteMapper() {} 
  RestSiteMapper(const RestSiteModelParams& mParams): RestSiteGeneral(mParams) {}

  virtual void FindMatches(const string& fileNameQuery, const string& fileNameTarget, MatchCandids& finalMatches); 

private:
};

class RestSiteDBMapper : public RestSiteGeneral 
{
public:
  RestSiteDBMapper() {} 
  RestSiteDBMapper(const RestSiteModelParams& mParams): RestSiteGeneral(mParams) {}

  virtual void FindMatches(const string& fileNameQuery, const string& fileNameTarget, MatchCandids& finalMatches); 
  void SetQuerySites(const string& fileName); 
  string GetQueryName(int readIdx) const; 
  void WriteMatchCandids(const map<int, map<int, int> >& candids) const; 

private:
  map<string, RSiteReads> m_queryReads;          /// Query Restriction Site reads per motif 
};

#endif //OPTIMAPALIGNUNIT_H
