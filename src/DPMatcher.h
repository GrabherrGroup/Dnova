#ifndef DPMATCHER_H
#define DPMATCHER_H

#include "RSiteReads.h"
#include "MappedInstance.h"

class DPMatcher{
public:
  DPMatcher() {} 

  float FindMatchScore(const MappedInstance& mapCandid, const RSiteReads& reads); 

private:
  int GetScore(int hCoord, int vCoord, const vector<vector<int>>& editGrid); 
  void SetScore(int hCoord, int vCoord, int setVal, vector<vector<int>>& editGrid); 

};

#endif //DPMATACHER_H
