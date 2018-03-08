#ifndef DPMATCHER_H
#define DPMATCHER_H

#include "RSiteReads.h"
#include "MappedInstance.h"

class MatchInfo {
public:
  MatchInfo(): m_numMatches(0), m_firstMatchPos1(-1), m_lastMatchPos1(-1),
               m_firstMatchPos2(-1), m_lastMatchPos2(-1), m_seqLen1(0), m_seqLen2(0) {}
  MatchInfo(int numMatch, int fmPos1, int lmPos1, int fmPos2, int lmPos2, int sLen1, int sLen2): 
    m_numMatches(numMatch), m_firstMatchPos1(fmPos1), m_lastMatchPos1(lmPos1),
    m_firstMatchPos2(fmPos2), m_lastMatchPos2(lmPos2), m_seqLen1(sLen1), m_seqLen2(sLen2) {}

  void SetNumMatches(int nm)      { m_numMatches = nm;      }
  void SetFirstMatchPos1(int fmp) { m_lastMatchPos1 = fmp;  }
  void SetLastMatchPos1(int lmp)  { m_lastMatchPos1 = lmp;  }
  void SetFirstMatchPos2(int fmp) { m_lastMatchPos2 = fmp;  }
  void SetLastMatchPos2(int lmp)  { m_lastMatchPos2 = lmp;  }
  void SetSeqLen1(int sl1)        { m_seqLen1 = sl1;        }
  void SetSeqLen2(int sl2)        { m_seqLen1 = sl2;        }
 
  int GetNumMatches() const       { return m_numMatches;     }
  int GetFirstMatchPos1() const   { return m_firstMatchPos1; }
  int GetLastMatchPos1() const    { return m_lastMatchPos1;  }
  int GetFirstMatchPos2() const   { return m_firstMatchPos2; }
  int GetLastMatchPos2() const    { return m_lastMatchPos2;  }
  int GetSeqLen1() const          { return m_seqLen1;        }
  int GetSeqLen2() const          { return m_seqLen2;        }

  float GetOverlapScore() const;
  float GetLocalIdenityScore() const { return ((float)m_numMatches/(m_lastMatchPos1-m_firstMatchPos1)); }

private:
  int m_numMatches;      /// Total number of bases matching between the two sequences
  int m_firstMatchPos1;  /// The position for last rsite value that matches in seq1
  int m_lastMatchPos1;   /// The position for last rsite value that matches in seq1
  int m_firstMatchPos2;  /// The position for last rsite value that matches in seq2
  int m_lastMatchPos2;   /// The position for last rsite value that matches in seq2
  int m_seqLen1;         /// Length of region to be matched from the first sequence
  int m_seqLen2;         /// Length of region to be matched from the second sequence
};

class DPMatcher{
public:
  DPMatcher() {} 

  float FindMatch(const Dmer& dm1, const Dmer& dm2, const RSiteReads& reads, MatchInfo& mInfo) const; 
  float FindMatch(int readIdx1, int readIdx2, int offset1, int offset2, bool matchDir, const RSiteReads& reads, MatchInfo& mInfo) const; 

private:
  int GetScore(int hCoord, int vCoord, const vector<vector<int>>& editGrid) const; 
  void SetScore(int hCoord, int vCoord, int setVal, vector<vector<int>>& editGrid) const; 
};


#endif //DPMATACHER_H
