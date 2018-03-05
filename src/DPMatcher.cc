#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "DPMatcher.h"

float DPMatcher::FindMatchScore(const Dmer& dm1, const Dmer& dm2,
    const RSiteReads& reads, MatchInfo& mInfo) {
  const RSiteRead& read1  = reads[dm1.Seq()];
  const RSiteRead& read2  = reads[dm2.Seq()];
  int readPos1            = dm1.Pos();
  int readPos2            = dm2.Pos();

  int minPos = min(readPos1, readPos2);
  readPos1  -= minPos;
  readPos2  -= minPos;

  bool isMatch = true;
  
  int length1 = read1.Size() - readPos1; 
  int length2 = read2.Size() - readPos2; 
  
  // lookup table for storing results of subproblems initialized with zeros
  vector<vector<int> > editGrid(length1, vector<int>(length2));
  int maxCell_score  = 0;
  int maxCell_hCoord = readPos1;
  int maxCell_vCoord = readPos2;
  int windowThresh = 2; //Parameterize
  for (int hCoord = 0; hCoord <= length1; hCoord++) {
    for (int vCoord = 0; vCoord <= length2; vCoord++) {
      if(abs(hCoord-vCoord)>windowThresh) { continue; }
      int hMoveTotal = GetScore(hCoord, vCoord-1, editGrid); 
      int vMoveTotal = GetScore(hCoord-1, vCoord, editGrid); 
      int dMoveAdd   = (read1[hCoord+readPos1]==read2[vCoord+readPos2])?1:0;
      int dMoveTotal = GetScore(hCoord-1, vCoord-1, editGrid) + dMoveAdd; 
      int currScore  = max(max(hMoveTotal, vMoveTotal), dMoveTotal);
      if(currScore > maxCell_score) { 
        maxCell_score  = currScore; 
        maxCell_hCoord = readPos1 + hCoord; 
        maxCell_vCoord = readPos2 + vCoord;
      }
      SetScore(hCoord, vCoord, currScore, editGrid);
    }
  }
  mInfo = MatchInfo(maxCell_score, maxCell_hCoord, maxCell_vCoord, length1, length2);
  //TODO: Overlap specifc score, this needs to be extended to more score types
  float score = (float)maxCell_score/min(length1, length2);
  return score;
}

int DPMatcher::GetScore(int hCoord, int vCoord, const vector<vector<int>>& editGrid) {
  if(hCoord < 0 || vCoord < 0) { return 0; } //Beyond limits
  if(hCoord < editGrid.size() && vCoord < editGrid[hCoord].size()) { 
    return editGrid[hCoord][vCoord];
  } else {
    return 0;  //Beyond limits
  }
}

void DPMatcher::SetScore(int hCoord, int vCoord, int setVal, vector<vector<int>>& editGrid) {
  if(hCoord>=0 && vCoord>=0 && hCoord < editGrid.size() && vCoord < editGrid[hCoord].size()) { 
    editGrid[hCoord][vCoord] = setVal;
  }
}
