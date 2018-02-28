#ifndef MAPINSTANCE_H
#define MAPINSTANCE_H


#include <map>
#include "Dmers.h"

class MappedInstance
{
public:
  MappedInstance(): m_rIdx1(-1), m_rIdx2(-1), m_rPos1(-1), m_rPos2(-1) {}
  MappedInstance(int idx1, int idx2, int rp1, int rp2): m_rIdx1(idx1), m_rIdx2(idx2), m_rPos1(rp1), m_rPos2(rp2) {}
  //Constructor to initialize from dmers
  MappedInstance(const Dmer& dm1, const Dmer& dm2): 
    m_rIdx1(dm1.Seq()), m_rIdx2(dm2.Seq()), m_rPos1(dm1.Pos()), m_rPos2(dm2.Pos()) {}

  int GetFirstReadIndex() const  { return m_rIdx1;       }
  int GetSecondReadIndex() const { return m_rIdx2;       }
  int GetFirstMatchPos() const   { return m_rPos1;       }
  int GetSecondMatchPos() const  { return m_rPos2;       }

  inline bool operator < (const MappedInstance& rhs) const {
    return(tie(m_rIdx1, m_rIdx2, m_rPos1)
       < tie(rhs.m_rIdx1, rhs.m_rIdx2, rhs.m_rPos1)); //keep the same order
  }

  inline bool operator == (const MappedInstance& rhs) const {
    return(tie(m_rIdx1, m_rIdx2, m_rPos1)
       == tie(rhs.m_rIdx1, rhs.m_rIdx2, rhs.m_rPos2)); //keep the same order
  }

  string ToString() const;

private:
  int m_rIdx1;        /// The index of the first read
  int m_rIdx2;        /// The index of the second read 
  int m_rPos1;        /// Position of match in first read 
  int m_rPos2;        /// Position of match in second read 
};

class MappedInstances
{
public:
  MappedInstances(): m_candids() {}
  
  //int NumCandids() const { return m_candids.isize(); }

  void AddCandidSort(int rIdx1, int rIdx2, int rPos1, int rPos2);
  
  string ToString() const; 

  const map<int, map<int, int > >& GetAllCandids() const { return m_candids; }
private:
  map<int, map<int, int> > m_candids;
};


#endif //MAPINSTANCE_H
