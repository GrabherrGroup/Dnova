#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <ctime>
#include "RestSiteAlignUnit.h"


int main( int argc, char** argv )
{
  commandArg<string> tFileCmmd("-t","Target input fasta file");
  commandArg<string> qFileCmmd("-q","Query input fasta file");
  commandArg<int> dmerCmmd("-d","dmer length", 6);
  commandArg<int> motifLenCmmd("-ml","Motif Length", 4);
  commandArg<int> motifCntCmmd("-mc","Number of motifs to use", 1);
  commandArg<bool> singleStrCmmd("-s", "1: if single strand or 0: if reverse complements should also be included", 0);
  commandArg<double> ndfcCmmd("-nc", "Coefficient to determine how much room to allow for differences in dmers", 1.0);
  commandArg<int>  coreCmmd("-n","Number of Cores to run with", 2);
  commandArg<string> appLogCmmd("-L","Application logging file","application.log");
  commandLineParser P(argc,argv);
  P.SetDescription("Find overlaps in restriction maps.");
  P.registerArg(tFileCmmd);
  P.registerArg(qFileCmmd);
  P.registerArg(dmerCmmd);
  P.registerArg(motifLenCmmd);
  P.registerArg(motifCntCmmd);
  P.registerArg(singleStrCmmd);
  P.registerArg(ndfcCmmd);
  P.registerArg(coreCmmd);
 
  P.parse();
  
  string tFileName  = P.GetStringValueFor(tFileCmmd);
  string qFileName  = P.GetStringValueFor(qFileCmmd);
  int dmerLen       = P.GetIntValueFor(dmerCmmd);
  int motifLen      = P.GetIntValueFor(motifLenCmmd);
  int motifCnt      = P.GetIntValueFor(motifCntCmmd);
  bool singleStrand = P.GetBoolValueFor(singleStrCmmd);
  double ndfCoef    = P.GetDoubleValueFor(ndfcCmmd);
  int numOfCores    = P.GetIntValueFor(coreCmmd);
    string logFile  = P.GetStringValueFor(appLogCmmd);

  FILE* pFile               = fopen(logFile.c_str(), "w");
  Output2FILE::Stream()     = pFile;
  FILELog::ReportingLevel() = logDEBUG4;
#if defined(FORCE_DEBUG)
//    FILELog::ReportingLevel() = logDEBUG4; 
#endif

  omp_set_num_threads(numOfCores); //The sort functions use OpenMP

  RestSiteModelParams mParams(singleStrand, motifLen, motifCnt, dmerLen, ndfCoef); 
  RestSiteDBMapper rsDBMapper(mParams);
  clock_t clock1_optiLoad, clock2_mapInstances;
  clock1_optiLoad = clock();
  rsDBMapper.FindMatches(qFileName, tFileName); 
  clock2_mapInstances = clock();

/*
  cout << "Report runtime duration: " << endl;
  cout << " Loading Optical Reads: "  
       << ((double) (clock2_mapInstances-clock1_optiLoad) / CLOCKS_PER_SEC)
       << endl;
  cout << " Finding Overlap Candidates: " 
       << ((double) (clock3_finalOverlaps-clock2_overlapCand) / CLOCKS_PER_SEC)
       << endl;
  cout << " Refining and finalizing overlaps: " 
       << ((double) (clock4_done - clock3_finalOverlaps) / CLOCKS_PER_SEC)
       << endl;
*/
  return 0;
}
