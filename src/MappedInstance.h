#ifndef MAPINSTANCE_H
#define MAPINSTANCE_H


#include <map>
#include "Dmers.h"

class MappedInstance
{
public:
  MappedInstance(): m_mappedIndex(-1), m_contactPos(-1), m_leftUncertain(-1), 
                    m_rightUncertain(-1), m_direction(false), m_orient(false)  {}
  MappedInstance(int mIdx, int cPos, int leftUc, int rightUc, bool dir, bool ori): m_mappedIndex(mIdx), m_contactPos(cPos), 
                 m_leftUncertain(leftUc), m_rightUncertain(rightUc), m_direction(dir), m_orient(ori) {}

  int getMappedIndex() const     { return m_mappedIndex;      }  
  int getContactPos() const      { return m_contactPos;       }  
  int getLeftUncertain() const   { return m_leftUncertain;    }  
  int getRightUncertain() const  { return m_rightUncertain;   }  
  int getDirection() const       { return (m_direction)?1:-1; }
  int getOrient() const          { return (m_orient)?1:-1;    }
  bool getOrientBool() const     { return m_orient;           }
  bool getDirectionBool() const  { return m_direction;        } 

  inline bool operator < (const MappedInstance& rhs) const {
    return(tie(m_mappedIndex, m_contactPos)
       < tie(rhs.m_mappedIndex, rhs.m_contactPos)); //keep the same order
  }

  inline bool operator == (const MappedInstance& rhs) const {
    return(tie(m_mappedIndex, m_contactPos)
       == tie(rhs.m_mappedIndex, rhs.m_contactPos)); //keep the same order
  }

  string ToString() const;

private:
    int     m_mappedIndex;    /// The index of the sequence to which this item refers to
    int     m_contactPos;     /// The index in the read where this mapping occurs from 
    int     m_leftUncertain;  /// Number of bases uncertain from the left of the mapping
    int     m_rightUncertain; /// Number of bases uncertain from the right of the mapping 
    bool    m_direction;      /// The direction of the overlap (true (+1): right  false (-1): left) 
    bool    m_orient;         /// The overlap orientat (Whether the reads are of the same strand) - same logic as direction
};

class MappedInstances
{
public:
  MappedInstances(): m_candids() {}
  
  //int NumCandids() const { return m_candids.isize(); }

  void AddCandidSort(int mIdx, int cPos, int leftUc, int rightUc, bool dir, bool ori);  
  
  string ToString() const; 

  const map<int, map<int, MappedInstance>>& GetAllCandids() const { return m_candids; }
private:
  map<int, map<int, MappedInstance>> m_candids;
};


#endif //MAPINSTANCE_H
