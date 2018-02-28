#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "DPMatcher.h"

float DPMatcher::FindMatchScore(const MappedInstance& mapCandid, const RSiteReads& reads) {
  const RSiteRead& read1  = reads[mapCandid.GetFirstReadIndex()]; 
  const RSiteRead& read2  = reads[mapCandid.GetSecondReadIndex()]; 
  int readPos1            = mapCandid.GetFirstMatchPos();
  int readPos2            = mapCandid.GetSecondMatchPos();

  bool isMatch = true;
  
  int length1 = read1.Size() - readPos1; 
  int length2 = read2.Size() - readPos2; 
  
  // lookup table for storing results of subproblems initialized with zeros
  vector<vector<int> > editGrid(length1, vector<int>(length2));
  int maxCell = 0;
  for (int hCoord = 0; hCoord <= length1; hCoord++) {
    for (int vCoord = 0; vCoord <= length2; vCoord++) {
      int hMoveTotal = GetScore(hCoord, vCoord-1, editGrid); 
      int vMoveTotal = GetScore(hCoord-1, vCoord, editGrid); 
      int dMoveAdd   = (read1[hCoord]==read2[vCoord])?1:0;
      int dMoveTotal = GetScore(hCoord-1, vCoord-1, editGrid) + dMoveAdd; 
      int currScore  = max(max(hMoveTotal, vMoveTotal), dMoveTotal);
      if(currScore > maxCell) { maxCell = currScore; }
      SetScore(hCoord, vCoord, currScore, editGrid);
    }
  }
  float score = (float)maxCell/max(length1, length2);
  return maxCell;
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
