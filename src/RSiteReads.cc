#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <sstream>
#include "RSiteReads.h"

string RSiteRead::ToString() const {
  return ToString(0);
}
string RSiteRead::ToString(int offset) const {
  stringstream ss;
  ss << PreDist() << ", ";
  for (int i=offset; i<m_dist.isize(); i++) {
    ss << m_dist[i] << ", ";
  }
  ss << PostDist();
  return ss.str();
}

void RSiteRead::Flip() {
  m_ori = -m_ori;
  svec<int> tmp;
  tmp.resize(m_dist.isize());
  int k = m_dist.isize()-1;
  for (int i=0; i<m_dist.isize(); i++) {
    tmp[k] = m_dist[i];
    k--;
  }
  m_dist = tmp;
  int tmp_pp; //Swap pre/postfix
  tmp_pp     = m_preDist;
  m_preDist  = m_postDist;
  m_postDist = tmp_pp;
}

void RSiteRead::GetCumulative(RSiteRead& cRead, int offset, bool dir) const {
  cRead.m_preDist  = m_preDist;
  cRead.m_postDist = m_postDist;
  cRead.m_name     = m_name;
  cRead.m_ori      = m_ori;
  
  svec<int> tmp;
  int readLen  = m_dist.isize();
  tmp.resize(readLen);
  int subreadLen = readLen - offset;
  if(!dir) { subreadLen = offset; }
  int acm      = 0; //Accumulation so far
  for (int cnt=0; cnt<subreadLen; cnt++) {
    if(dir) {
      int addValue = m_dist[offset+cnt];
      acm += addValue;
      tmp[offset+cnt] = acm;
    } else { 
      int addValue = m_dist[offset-cnt]; 
      acm += addValue;
      tmp[offset-cnt] = acm;
    }
  }
  cRead.m_dist = tmp;
}

int RSiteReads::AddRead(const RSiteRead& rr) {
  m_rReads.push_back(rr);
  m_readCount++;
  return m_readCount-1;
}

string RSiteReads::ToString() const {
  string strOut;
  for(int i=0; i<m_readCount; i++) {
    strOut += m_rReads[i].ToString();
    strOut += "\n";
  }
  return strOut;
}
