#ifndef DPMATCHER_H
#define DPMATCHER_H

#include "RSiteReads.h"
#include "MappedInstance.h"

class DPMatcher{
public:
  DPMatcher() {} 

  float FindMatchScore(const Dmer& dm1, const Dmer& dm2, const RSiteReads& reads); 

private:
  int GetScore(int hCoord, int vCoord, const vector<vector<int>>& editGrid); 
  void SetScore(int hCoord, int vCoord, int setVal, vector<vector<int>>& editGrid); 

};

#endif //DPMATACHER_H
