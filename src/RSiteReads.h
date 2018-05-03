#ifndef RSITEREADS_H
#define RSITEREADS_H

#include "ryggrad/src/base/SVector.h"

class RSiteRead
{
public:
  RSiteRead() {
    m_ori = 1;
  }
  
  int operator[](int idx) const    { return m_dist[idx];    }
  svec<int> & Dist()               { return m_dist;         }
  string & Name()                  { return m_name;         }
  int& PreDist()                   { return m_preDist;      }
  int& PostDist()                  { return m_postDist;     }
  int & Ori()                      { return m_ori;          }
  const svec<int> & Dist() const   { return m_dist;         }
  const int& PreDist() const       { return m_preDist;      }
  const int& PostDist() const      { return m_postDist;     }
  const string & Name() const      { return m_name;         }
  const int & Ori() const          { return m_ori;          }
  int Size() const                 { return m_dist.isize(); }

  string ToString() const;
  string ToString(int offset) const;
  void Flip();
  void GetCumulative(RSiteRead& cRead, int offset, bool dir) const;

private:
  svec<int> m_dist;       /// Distmer values
  int m_preDist;          /// Number of bits prior to the start of the first distmer value
  int m_postDist;         /// Number of bits left over after the last distmer value
  string m_name;          /// Name of optiRead
  int m_ori;              /// Orientation
};

class RSiteReads 
{
public:
  // Default Ctor
  RSiteReads(): m_readCount(0), m_rReads() {}

  void Reserve(int size)                     { m_rReads.reserve(size);   }

  const RSiteRead& operator[](int idx) const { return m_rReads[idx];     }
  RSiteRead& operator[](int idx)             { return m_rReads[idx];     }
  int NumReads() const                       { return m_readCount;       }

  int AddRead(const RSiteRead& rr); 
  string ToString() const;

private:
  int m_readCount;
  svec<RSiteRead> m_rReads;  /// List of site reads
};

#endif //RSITERREADS_H
