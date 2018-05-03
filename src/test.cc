#include <ctime>
#include "RestSiteAlignUnit.h"

int GetScore(int hCoord, int vCoord, const vector<vector<int>>& editGrid) {
  if(hCoord < 0 || vCoord < 0) { return 0; } //Beyond limits
  if(hCoord < editGrid.size() && vCoord < editGrid[hCoord].size()) { 
    return editGrid[hCoord][vCoord];
  } else {
    return 0;  //Beyond limits
  }
}

void SetScore(int hCoord, int vCoord, int setVal, vector<vector<int>>& editGrid) {
  if(hCoord>=0 && vCoord>=0 && hCoord < editGrid.size() && vCoord < editGrid[hCoord].size()) { 
    editGrid[hCoord][vCoord] = setVal;
  }
}

int main( int argc, char** argv )
{
  vector<int> read1  = {18,22,31,154,195,214,247,257,277,341,448,543,562,574,634,751,1749,1773,1783,1883};
  vector<int> read2  = {18,22,31,150,193,211,243,252,270,327,452,460,532,632,640,1707,1729,1759,1770,2040}; 
  
  int length1 = read1.size();
  int length2 = read2.size();
  // lookup table for storing results of subproblems initialized with zeros
  vector<vector<int> > editGrid(length1, vector<int>(length2));
  int maxCell_score  = 0;
  int maxCell_hCoord = 0;
  int maxCell_vCoord = 0;
  int windowThresh = 1000;   // Parameterize
  for (int hCoord = 0; hCoord < length1; hCoord++) {
    for (int vCoord = 0; vCoord < length2; vCoord++) {
      if(abs(hCoord-vCoord)>windowThresh) { continue; }
      int hMoveTotal = GetScore(hCoord, vCoord-1, editGrid); 
      int vMoveTotal = GetScore(hCoord-1, vCoord, editGrid); 
      int readVal1 = read1[hCoord];
      int readVal2 = read2[vCoord];
      float average = (readVal1 + readVal2) / 2.0;
      int deviation = sqrt(average*0.08)*3;
      int dMoveAdd = (readVal1>=average-deviation && readVal1<=average+deviation 
                      && readVal2>=average-deviation && readVal2<=average+deviation)? 1 : 0;
      int dMoveTotal = GetScore(hCoord-1, vCoord-1, editGrid) + dMoveAdd; 
      int currScore  = max(max(hMoveTotal, vMoveTotal), dMoveTotal);
      if(currScore > maxCell_score) { 
        maxCell_score  = currScore; 
        maxCell_hCoord = hCoord; 
        maxCell_vCoord = vCoord;
      }
      SetScore(hCoord, vCoord, currScore, editGrid);
    }
  }
  cout<<maxCell_score<<endl;
}

