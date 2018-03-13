#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <sstream>
#include "ryggrad/src/base/Logger.h"
#include "DPMatcher.h"

float DPMatcher::FindMatch(const Dmer& dm1, const Dmer& dm2,
  const RSiteReads& reads, MatchInfo& mInfo) const {
  MatchInfo mInfo1, mInfo2;
  FindMatch(dm1.Seq(), dm2.Seq(), dm1.Pos(), dm2.Pos(), true, reads, mInfo1);
  FindMatch(dm1.Seq(), dm2.Seq(), dm1.Pos(), dm2.Pos(), false, reads, mInfo2);
  int totNumMatches = mInfo1.GetNumMatches() + mInfo2.GetNumMatches();
  int seqLen1       = mInfo1.GetSeqLen1() + mInfo2.GetSeqLen1();
  int seqLen2       = mInfo1.GetSeqLen2() + mInfo2.GetSeqLen2();
  mInfo = MatchInfo(totNumMatches, mInfo2.GetFirstMatchPos1(), mInfo1.GetLastMatchPos1(),
                    mInfo2.GetFirstMatchPos2(), mInfo1.GetLastMatchPos2(), seqLen1, seqLen2);
  FILE_LOG(logDEBUG4) << "Overall MatchInfo: " << mInfo.ToString();
  return (mInfo.GetOverlapScore());
}

float DPMatcher::FindMatch(int readIdx1, int readIdx2, int offset1, int offset2, bool matchDir,
    const RSiteReads& reads, MatchInfo& mInfo) const {

  const RSiteRead& read1  = reads[readIdx1];
  const RSiteRead& read2  = reads[readIdx2];

  int length1 = read1.Size() - offset1; 
  int length2 = read2.Size() - offset2; 
  if(!matchDir) { //moving backwards to find match
    length1 = offset1 + 1; 
    length2 = offset2 + 1; 
  }
  // lookup table for storing results of subproblems initialized with zeros
  vector<vector<int> > editGrid(length1, vector<int>(length2));
  int maxCell_score  = 0;
  int maxCell_hCoord = offset1;
  int maxCell_vCoord = offset2;
  if(!matchDir) {
    maxCell_hCoord = offset1;
    maxCell_vCoord = offset2;
  }
  int windowThresh = 5; //Parameterize
  for (int hCoord = 0; hCoord < length1; hCoord++) {
    for (int vCoord = 0; vCoord < length2; vCoord++) {
      if(abs(hCoord-vCoord)>windowThresh) { continue; }
      int hMoveTotal = GetScore(hCoord, vCoord-1, editGrid) - 1; 
      int vMoveTotal = GetScore(hCoord-1, vCoord, editGrid) - 1; 
      int dMoveAdd   = (read1[hCoord+offset1]==read2[vCoord+offset2])?1:0;
      if(!matchDir) {
        dMoveAdd = (read1[offset1-hCoord]==read2[offset2-vCoord])?1:0;
      }
      int dMoveTotal = GetScore(hCoord-1, vCoord-1, editGrid) + dMoveAdd; 
      int currScore  = max(max(hMoveTotal, vMoveTotal), dMoveTotal);
      if(currScore > maxCell_score) { 
        maxCell_score  = currScore; 
        if(matchDir) {
          maxCell_hCoord = offset1 + hCoord; 
          maxCell_vCoord = offset2 + vCoord;
        } else {
          maxCell_hCoord = offset1 - hCoord; 
          maxCell_vCoord = offset2 - vCoord;
        }
      }
      SetScore(hCoord, vCoord, currScore, editGrid);
    }
  }
  if(matchDir) {
    mInfo = MatchInfo(maxCell_score, offset1, maxCell_hCoord, offset2, maxCell_vCoord, length1, length2);
  } else {
    mInfo = MatchInfo(maxCell_score, maxCell_hCoord, offset1, maxCell_vCoord, offset2, length1, length2);
  }
  FILE_LOG(logDEBUG4) << "MatchInfo: " << mInfo.ToString();
  //TODO: Overlap specifc score, this needs to be extended to more score types
  return mInfo.GetOverlapScore();
}

int DPMatcher::GetScore(int hCoord, int vCoord, const vector<vector<int>>& editGrid) const {
  if(hCoord < 0 || vCoord < 0) { return 0; } //Beyond limits
  if(hCoord < editGrid.size() && vCoord < editGrid[hCoord].size()) { 
    return editGrid[hCoord][vCoord];
  } else {
    return 0;  //Beyond limits
  }
}

void DPMatcher::SetScore(int hCoord, int vCoord, int setVal, vector<vector<int>>& editGrid) const {
  if(hCoord>=0 && vCoord>=0 && hCoord < editGrid.size() && vCoord < editGrid[hCoord].size()) { 
    editGrid[hCoord][vCoord] = setVal;
  }
}

float MatchInfo::GetOverlapScore() const { 
  float containmentScore = (float)m_numMatches/min(m_seqLen1, m_seqLen2);
  return containmentScore;
//  return max(containmentScore, 0.0);
}

string MatchInfo::ToString() const {
  stringstream ss;
  ss << m_numMatches << " " << m_firstMatchPos1 << " " << m_lastMatchPos1
     << " " << m_firstMatchPos2 << " " << m_lastMatchPos2 << " " << m_seqLen1 
     << " " << m_seqLen2 << " "  << GetOverlapScore();  
  return ss.str();
}

