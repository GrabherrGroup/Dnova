#ifndef DPMATCHER_H
#define DPMATCHER_H

#include "RSiteReads.h"
#include "MappedInstance.h"

class MatchInfo {
public:
  MatchInfo(): m_numMatches(0), m_lastMatchPos1(-1), m_lastMatchPos2(-1), m_seqLen1(0), m_seqLen2(0) {}
  MatchInfo(int numMatch, int lmPos1, int lmPos2, int sLen1, int sLen2): 
    m_numMatches(numMatch), m_lastMatchPos1(lmPos1), m_lastMatchPos2(lmPos2),
    m_seqLen1(sLen1), m_seqLen2(sLen2) {}

  void SetNumMatches(int nm)     { m_numMatches = nm;      }
  void SetLastMatchPos1(int lmp) { m_lastMatchPos1 = lmp;  }
  void SetLastMatchPos2(int lmp) { m_lastMatchPos1 = lmp;  }
  void SetSeqLen1(int sl1)       { m_seqLen1 = sl1;        }
  void SetSeqLen2(int sl2)       { m_seqLen1 = sl2;        }
 
  int GetNumMatches()            { return m_numMatches;    }
  int GetLastMatchPos1()         { return m_lastMatchPos1; }
  int GetLastMatchPos2()         { return m_lastMatchPos2; }
  int GetSeqLen1()               { return m_seqLen1;       }
  int GetSeqLen2()               { return m_seqLen1;       }

  float GetOverlapScore() { return ((float)m_numMatches/min(m_seqLen1, m_seqLen2)); }

private:
  int m_numMatches;     /// Total number of bases matching between the two sequences
  int m_lastMatchPos1;  /// The position for last rsite value that matches
  int m_lastMatchPos2;  /// The position for last rsite value that matches
  int m_seqLen1;        /// Length of region to be matched from the first sequence
  int m_seqLen2;        /// Length of region to be matched from the second sequence
};

class DPMatcher{
public:
  DPMatcher() {} 

  float FindMatchScore(const Dmer& dm1, const Dmer& dm2, const RSiteReads& reads, MatchInfo& mInfo); 

private:
  int GetScore(int hCoord, int vCoord, const vector<vector<int>>& editGrid); 
  void SetScore(int hCoord, int vCoord, int setVal, vector<vector<int>>& editGrid); 

};


#endif //DPMATACHER_H
