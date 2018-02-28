#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "MappedInstance.h"


string MappedInstance::ToString() const {
  stringstream ss;
  ss << GetFirstReadIndex() << "\t" <<  GetFirstMatchPos()  << "\t" 
     << GetSecondReadIndex() << "\t" << GetSecondMatchPos(); 
  return ss.str();
}

void MappedInstances::AddCandidSort(int rIdx1, int rIdx2, int rPos1, int rPos2) {
/*
  // Make sure that rIdx1, rIdx2 are in increasing order (so that sorting will bring all relevant pairs together)
  if(rIdx1>rIdx2) {
    int temp = rIdx1;
    rIdx1 = rIdx2;
    rIdx2 = temp; // Swap rIdx1 & rIdx2
    temp  = rPos1;
    rPos1 = rPos2;
    rPos2 = temp;
  }
  m_candids.push_back(MappedInstance(rIdx1, rIdx2, rPos1, rPos2));
*/
  m_candids[rIdx1][rIdx2]++;
  FILE_LOG(logDEBUG3) << "Adding overlap candidate: " << MappedInstance(rIdx1, rIdx2, rPos1, rPos2).ToString();
}
 
string MappedInstances::ToString() const {
  stringstream ss;
  for(auto const &ent1 : m_candids) {
    for(auto const &ent2 : ent1.second) {
      ss << ent1.first << " " << ent2.first << " " <<  ent2.second << endl;
    }
  }
  return ss.str();
}


