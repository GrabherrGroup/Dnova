#ifndef RESTSITECOREUNIT_H
#define RESTSITECOREUNIT_H

#include <string>
#include "RSiteReads.h"
#include "Dmers.h"
#include "DPMatcher.h"
#include "MappedInstance.h"

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
  float  IndelErr() const          { return m_insertionErr + m_substitutionErr; }

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
                      int dmerLength=6, float cndfCoef1=2.0, float cndfCoef2=1.0, 
                      float sThresh =0.2, const vector<char>& alphabet= {'A', 'C', 'G', 'T' }) 
                     :m_singleStrand(singleStrand), m_motifLength(motifLength), m_numOfMotifs(numOfMotifs),
                      m_dmerLength(dmerLength), m_cndfCoef1(cndfCoef1), m_cndfCoef2(cndfCoef2), 
                      m_scoreThresh(sThresh), m_alphabet(alphabet) { }

  bool   IsSingleStrand() const        { return m_singleStrand;    }
  int    MotifLength() const           { return m_motifLength;     }  
  int    NumOfMotifs() const           { return m_numOfMotifs;     } 
  int    DmerLength() const            { return m_dmerLength;      }
  float  CNDFCoef1() const             { return m_cndfCoef1;       }
  float  CNDFCoef2() const             { return m_cndfCoef2;       }
  float  ScoreThreshold() const        { return m_scoreThresh;     }
  int    AlphabetSize() const          { return m_alphabet.size(); }
  const vector<char>& Alphabet() const { return m_alphabet;        }

  void ChangeNumOfMotifs(int motifCnt) { m_numOfMotifs = motifCnt; }
private: 
  bool    m_singleStrand;   /// Flag specifying whether the reads are single or double strand
  int     m_motifLength;    /// Length of each motif
  int     m_numOfMotifs;    /// Number of motifs to generate/use
  int     m_dmerLength;     /// The length of distmers to use for seed finding
  float   m_cndfCoef1;      /// Cumulative Normal Distribution Function coefficeint used for estimating similarity at filtering stage 
  float   m_cndfCoef2;      /// Cumulative Normal Distribution Function coefficeint used for estimating similarity at  refinement stage
  float   m_scoreThresh;    /// Score threshold for accepting alignment refinement 
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
  const RSiteReads& Reads() const          { return m_rReads;       }

  string RSToString(int rIdx, int offset) const; //Convert RestSite read to string from given offset 
  string RSToString(const Dmer& dmer) const;     // read index and offset provided as dmer object

  int  CreateRSitesPerString(const string& origString, const string& origName, RSiteReads& reads, bool addRC) const; 

  void BuildDmers(); 
  int FindMapInstances(float indelVariance, map<int, map<int,vector<int>>>& checkedSeqs) const; 
  int FindSingleReadMapInstances(const RSiteRead& read, int rIdx, float indelVariance, map<int, map<int,vector<int>>>& checkedSeqs) const; 
  int HandleMappingInstance(const svec<Dmer>& dmers, float indelVariance, map<int, map<int,vector<int>>>& checkedSeqs, bool acceptSameIdx) const; 
  void ValidateMatch(const Dmer& dmer1, const Dmer& dmer2, float indelVariance, MatchInfo& matchInfo, float& side1Score, float& side2Score) const;
  void WriteMatchBasic(const Dmer& dm1, const Dmer& dm2, const MatchInfo& matchInfo, float& side1Score, float& side2Score) const; 
  void WriteMatchPAF(const Dmer& dm1, const Dmer& dm2, const MatchInfo& matchInfo, float& side1Score, float& side2Score) const;
  int GetBasePos(int seqIdx, int rsPos, bool inclusive) const; 
  int GetBasePos(const Dmer& dm, int rsPos, bool inclusive) const;
  float GetThresholdScore() const; 
  float GetRandomMatchProb() const;

protected:
  RSiteReads& Reads()             { return m_rReads; }

private:
  string m_motif;                    /// Vector of all motifs for which restriction site reads have been generated
  RestSiteModelParams m_modelParams; /// Model Parameters
  RestSiteDataParams m_dataParams;   /// Data Parameters
  double  m_totalSiteCnt;            /// The total of restriction site count over all reads
  RSiteReads m_rReads;               /// Restriction Site reads per motif
  Dmers  m_dmers;                    /// To build dmers from restriction site reads
};

#endif //RESTSITECOREUNIT_H

