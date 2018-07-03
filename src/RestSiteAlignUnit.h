#ifndef OPTIMAPALIGNUNIT_H
#define OPTIMAPALIGNUNIT_H

#include <string>
#include <parallel/algorithm>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"
#include "RSiteReads.h"
#include "MappedInstance.h"
#include "RestSiteCoreUnit.h"

class RestSiteGeneral 
{
public:
  RestSiteGeneral(): m_rsaCores(), m_motifs(), m_modelParams(), m_dataParams() {}
  RestSiteGeneral(const RestSiteModelParams& mParams): m_motifs(), m_modelParams(mParams), m_dataParams() {}

  /* Generate Permutation of the given alphabet to reach number of motifs required */
  void GenerateMotifs();  
  bool ValidateMotif(const string& motif, const vector<char>& alphabet, const map<char, char>& RCs) const; 
  void SetTargetSites(const string& fileName, bool addRC); 
  string GetTargetName(int readIdx) const;

  virtual void WriteMatchCandids(const map<int, map<int, int> >& candids) const; 
  virtual void FindMatches(const string& fileNameQuery, const string& fileNameTarget) = 0; 

protected:
  void CartesianPower(const vector<char>& input, unsigned k, vector<vector<char>>& result) const; 
  map<string, RestSiteMapCore> m_rsaCores;   /// Mapping engine (core data and functionality) per motif
  svec<string> m_motifs;                     /// Vector of all motifs for which restriction site reads have been generated
  RestSiteModelParams m_modelParams;         /// Model Parameters
  RestSiteDataParams m_dataParams;           /// Model Parameters
};

class RestSiteMapper : public RestSiteGeneral 
{
public:
  RestSiteMapper() {} 
  RestSiteMapper(const RestSiteModelParams& mParams): RestSiteGeneral(mParams) {}

  virtual void FindMatches(const string& fileNameQuery, const string& fileNameTarget); 

private:
};

#endif //OPTIMAPALIGNUNIT_H
