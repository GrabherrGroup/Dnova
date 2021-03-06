#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <ctime>
#include "RestSiteAlignUnit.h"


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input fasta file");
  commandArg<int> dmerCmmd("-d","dmer length", 4);
  commandArg<int> motifLenCmmd("-ml","Motif Length", 4);
  commandArg<int> motifCntCmmd("-mc","Number of motifs to use", 1);
  commandArg<bool> singleStrCmmd("-s", "1: if single strand or 0: if reverse complements should also be included", 0);
  commandArg<double> ndfcCmmd1("-nc1", "Coefficient to determine how much room to allow for differences in dmers in filtering stage", 2.5);
  commandArg<double> ndfcCmmd2("-nc2", "Coefficient to determine how much room to allow for differences in dmers in refinement stage", 1.0);
  commandArg<double> sThreshCmmd("-t", "Threshold of score for accepting a mapping at refinement stage. Default will be internally computed", -1.0);
  commandArg<int>  coreCmmd("-n","Number of Cores to run with", 2);
  commandArg<string> appLogCmmd("-L","Application logging file","application.log");
  commandLineParser P(argc,argv);
  P.SetDescription("Find overlaps in restriction maps.");
  P.registerArg(fileCmmd);
  P.registerArg(dmerCmmd);
  P.registerArg(motifLenCmmd);
  P.registerArg(motifCntCmmd);
  P.registerArg(singleStrCmmd);
  P.registerArg(ndfcCmmd1);
  P.registerArg(ndfcCmmd2);
  P.registerArg(sThreshCmmd);
  P.registerArg(coreCmmd);
 
  P.parse();
  
  string fileName   = P.GetStringValueFor(fileCmmd);
  int dmerLen       = P.GetIntValueFor(dmerCmmd);
  int motifLen      = P.GetIntValueFor(motifLenCmmd);
  int motifCnt      = P.GetIntValueFor(motifCntCmmd);
  bool singleStrand = P.GetBoolValueFor(singleStrCmmd);
  double ndfCoef1   = P.GetDoubleValueFor(ndfcCmmd1);
  double ndfCoef2   = P.GetDoubleValueFor(ndfcCmmd2);
  double scoreThresh= P.GetDoubleValueFor(sThreshCmmd);
  int numOfCores    = P.GetIntValueFor(coreCmmd);
    string logFile  = P.GetStringValueFor(appLogCmmd);

  FILE* pFile               = fopen(logFile.c_str(), "w");
  Output2FILE::Stream()     = pFile;
  FILELog::ReportingLevel() = logINFO;
  //FILELog::ReportingLevel() = logDEBUG4;
#if defined(FORCE_DEBUG)
//    FILELog::ReportingLevel() = logDEBUG4; 
#endif

  omp_set_num_threads(numOfCores); //The sort functions use OpenMP

  RestSiteModelParams mParams(singleStrand, motifLen, motifCnt, dmerLen, ndfCoef1, ndfCoef2, scoreThresh); 
  RestSiteMapper rsMapper(mParams);

  clock_t clock1_optiLoad, clock2_overlapCand, clock3_finalOverlaps, clock4_done;
  // 1a. Populate the motifs 
  // 1b. Construct Restriction-site Reads
  clock1_optiLoad = clock();

  // 2. Build Optimers and find those that share a seed as cadidates for overlap detection 
  rsMapper.FindMatches("", fileName); 
  clock2_overlapCand = clock();
 
  // 3. Take the overlap candidates and refine to remove false positives
  clock3_finalOverlaps = clock();
  clock4_done = clock();


  cout << "Report runtime duration: " << endl;
  cout << " Loading Optical Reads: "  
       << ((double) (clock2_overlapCand-clock1_optiLoad) / CLOCKS_PER_SEC)
       << endl;
  cout << " Finding Overlap Candidates: " 
       << ((double) (clock2_overlapCand-clock1_optiLoad) / CLOCKS_PER_SEC)
       << endl;
  cout << " Refining and finalizing overlaps: " 
       << ((double) (clock3_finalOverlaps-clock2_overlapCand) / CLOCKS_PER_SEC)
       << endl;

  return 0;
}
