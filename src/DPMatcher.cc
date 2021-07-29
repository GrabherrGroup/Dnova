#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <sstream>
#include "ryggrad/src/base/Logger.h"
#include "DPMatcher.h"
#include <math.h>

float DPMatcher::FindMatch(const Dmer& dm1, const Dmer& dm2,
                           const RSiteReads& reads, float indelVariance, float cndfCoef, MatchInfo& mInfo,
                           float& side1Score, float& side2Score) const { 
  MatchInfo mInfo1, mInfo2;
  FindMatch(dm1.Seq(), dm2.Seq(), dm1.Pos(), dm2.Pos(), true, reads, indelVariance, cndfCoef, mInfo1);
  FindMatch(dm1.Seq(), dm2.Seq(), dm1.Pos(), dm2.Pos(), false, reads, indelVariance, cndfCoef, mInfo2);
  int totNumMatches = mInfo1.GetNumMatches() + mInfo2.GetNumMatches() - 1; //Subtracting 1 to cater for double-counting
  int seqLen1       = mInfo1.GetSeqLen1() + mInfo2.GetSeqLen1();
  int seqLen2       = mInfo1.GetSeqLen2() + mInfo2.GetSeqLen2();
  mInfo = MatchInfo(totNumMatches, mInfo2.GetFirstMatchPos1(), mInfo1.GetLastMatchPos1(),
                    mInfo2.GetFirstMatchPos2(), mInfo1.GetLastMatchPos2(), seqLen1, seqLen2);
  side1Score = mInfo1.GetIdentScore();
  side2Score = mInfo2.GetIdentScore();
  FILE_LOG(logDEBUG4) << "Overall MatchInfo: " << mInfo.ToString();
  return (mInfo.GetIdentScore());
}

//Cumulative DPMatcher
float DPMatcher::FindMatch(int readIdx1, int readIdx2, int offset1, int offset2, bool matchDir,
    const RSiteReads& reads, float indelVariance, float cndfCoef, MatchInfo& mInfo) const {

  const RSiteRead& read1  = reads[readIdx1];
  const RSiteRead& read2  = reads[readIdx2];
  
  int baseLength1 = LengthOfBases(readIdx1, offset1, matchDir, reads);
  int baseLength2 = LengthOfBases(readIdx2, offset2, matchDir, reads);
  int length1 = read1.Size() - offset1; 
  int length2 = read2.Size() - offset2; 
  if(!matchDir) { //moving backwards to find match
    length1 = offset1 + 1; 
    length2 = offset2 + 1; 
  }
  if(baseLength1 <= baseLength2) { // Change length2
    length2 = GetRSiteLenForBaseLength(readIdx2, offset2, matchDir, baseLength1, reads);
  } else { // change length1
    length1 = GetRSiteLenForBaseLength(readIdx1, offset1, matchDir, baseLength2, reads);
  }

  int maxCell_score  = 0;
  int maxCell_coord1 = offset1;
  int maxCell_coord2 = offset2;

  int coord1=0, coord2=0;
  bool coord1Changed=true, coord2Changed=true;
  int readCmt1=0, readCmt2=0; //Read value cumulations
  while(coord1 < length1 && coord2 < length2) {
    int readVal1 = read1[offset1+coord1];
    int readVal2 = read2[offset2+coord2];
    if(!matchDir) {
      readVal1 = read1[offset1-coord1];
      readVal2 = read2[offset2-coord2];
    }
    if(coord1Changed) { readCmt1 += readVal1; }
    if(coord2Changed) { readCmt2 += readVal2; }
    float average = (readCmt1 + readCmt2) / 2.0;
    int deviation = sqrt(average*indelVariance)*cndfCoef;
    FILE_LOG(logDEBUG4) << "Using deviation value:" << deviation << " for sites: " << readCmt1 << " " << readCmt2; 
    int matched = (readCmt1>=average-deviation && readCmt1<=average+deviation 
                    && readCmt2>=average-deviation && readCmt2<=average+deviation)? 1 : 0;
    if(matched) { 
      FILE_LOG(logDEBUG4) << "Matched " << coord1 << " " << coord2; 
      maxCell_score++;
      if(matchDir) {
        maxCell_coord1 = offset1 + coord1; 
        maxCell_coord2 = offset2 + coord2;
      } else {
        maxCell_coord1 = offset1 - coord1; 
        maxCell_coord2 = offset2 - coord2;
      }
      coord1++;
      coord2++;
      readCmt1 = 0;
      readCmt2 = 0;
      coord1Changed = true;
      coord2Changed = true;
      continue;
    }
    if(readCmt1 > readCmt2) {
      coord2++; //move other coordinate up only
      coord1Changed = false;
      coord2Changed = true;
    } else {
      coord1++; //move other coordinate up only
      coord1Changed = true;
      coord2Changed = false;
    }
  }
  if(matchDir) {
    mInfo = MatchInfo(maxCell_score, offset1, maxCell_coord1, offset2, maxCell_coord2, length1, length2);
  } else {
    mInfo = MatchInfo(maxCell_score, maxCell_coord1, offset1, maxCell_coord2, offset2, length1, length2);
  }
  FILE_LOG(logDEBUG4) << "MatchInfo: " << mInfo.ToString();

  return mInfo.GetIdentScore();
}

int DPMatcher::LengthOfBases(int readIdx, int offset, bool dir, const RSiteReads& reads) const {
  const RSiteRead& read  = reads[readIdx];
  int totBaseLen         = 0;
  int siteLen            = read.Size() - offset; 
  if(!dir) { //moving backwards 
    siteLen = offset + 1; 
  }
  for(int idx=0; idx<siteLen; idx++) {
    int toAdd = read[offset + idx];
    if(!dir) { toAdd = read[offset - idx]; }
    totBaseLen += toAdd;
  }
  return totBaseLen;
}

int DPMatcher::GetRSiteLenForBaseLength(int readIdx, int offset, bool dir, int totLength, const RSiteReads& reads) const {
  const RSiteRead& read  = reads[readIdx];
  int totBaseLen         = 0;
  int siteLen            = read.Size() - offset; 
  if(!dir) { //moving backwards 
    siteLen = offset + 1; 
  }
  for(int idx=0; idx<siteLen; idx++) {
    int toAdd = read[offset + idx];
    if(!dir) { toAdd = read[offset - idx]; }
    totBaseLen += toAdd;
    if(totBaseLen >= totLength) {
      return idx;
    }
  }
}
 
float MatchInfo::GetOverlapScore() const { 
  float containmentScore = (float)m_numMatches/min(m_seqLen1, m_seqLen2);
  return containmentScore;
//  return max(containmentScore, 0.0);
}

float MatchInfo::GetIdentScore() const { 
  float identScore = (float)m_numMatches/max(m_seqLen1, m_seqLen2);
  return identScore;
}

string MatchInfo::ToString() const {
  stringstream ss;
  ss << m_numMatches << " " << m_firstMatchPos1 << " " << m_lastMatchPos1
     << " " << m_firstMatchPos2 << " " << m_lastMatchPos2 << " " << m_seqLen1 
     << " " << m_seqLen2 << " "  << GetOverlapScore();  
  return ss.str();
}

