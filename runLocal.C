int runLocal(TString listOfFiles = "local_filelist.txt", Bool_t bIsPbPb = 1, TString sYear = "2010", Bool_t bMCData = 0, UInt_t numLocalFiles = 50)
{
  TString sFuncAnalysis = "AddAnalysisStrangeJetsEmcal";

  TChain* chain = reinterpret_cast<TChain*>(gInterpreter->ProcessLine(Form(".x %s(\"%s\", %d)", gSystem->ExpandPathName("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C"), listOfFiles.Data(), numLocalFiles)));
  if(!chain)
  {
    printf("runLocal: Cannot create the chain\n");
    return 1;
  }

  // Make the analysis manager and connect event handlers
  AliAnalysisManager* mgr = new AliAnalysisManager("LocalTrain", "train");
  mgr->SetUseProgressBar(1, 25);

  // All analysis tasks are loaded and configured here
  long iLoadAnalysis = reinterpret_cast<long>(gInterpreter->ProcessLine(Form(".x %s(%d, \"%s\", %d)", gSystem->ExpandPathName(Form("%s.C", sFuncAnalysis.Data())), bIsPbPb, sYear.Data(), bMCData)));
  if(iLoadAnalysis)
  {
    printf("runLocal: Error in %s, exit code: %ld\n", sFuncAnalysis.Data(), iLoadAnalysis);
    return 1;
  }
  //return 0;

  // counting duration of analysis
  TStopwatch timer;
  timer.Start();

  if(!mgr->InitAnalysis())
  {
    printf("runLocal: Error in InitAnalysis");
    return 1;
  }
  mgr->PrintStatus();
  //return 0;
  mgr->StartAnalysis("local", chain);

  // stop counting and print out
  timer.Stop();
  timer.Print();

  return 0;
}
