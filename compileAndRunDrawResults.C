void compileAndRunDrawResults(TString sNameFileIn, Int_t V0Type = 0, TString sNameFileEff = "", TString sNameFileFD = "", Int_t iSwitchMacro = 0, TString sFlag = "", Int_t iMode = 0)
{
  Bool_t bLoopSys = 0; // will call DrawResultsLoopSys instead of DrawResults in order to loop over cut variations
  TString sPathMacros = ".";
  TString sPathXi = "./official/spectraXi";
  TString sFunctionName = "";

  if(iSwitchMacro == 0) {
    sFunctionName = bLoopSys ? "DrawResultsLoopSys" : "DrawResults";
    reinterpret_cast<long>(gInterpreter->ProcessLine(Form(".x %s(\"%s\", %d, \"%s\", \"%s\",  \"%s\",  %d)", gSystem->ExpandPathName(Form("%s.C", sFunctionName.Data())), sNameFileIn.Data(), V0Type, sNameFileEff.Data(), sNameFileFD.Data(), sFlag.Data(), iMode)));
  } else if(iSwitchMacro == 1) {
    sFunctionName = "DrawRatios";
    reinterpret_cast<long>(gInterpreter->ProcessLine(Form(".x %s(\"%s\", %d)", gSystem->ExpandPathName(Form("%s.C", sFunctionName.Data())), sFlag.Data(), iMode)));
  } else if(iSwitchMacro == 2) {
    sFunctionName = "V0TreeAnalysis";
    reinterpret_cast<long>(gInterpreter->ProcessLine(Form(".x %s(\"%s\")", gSystem->ExpandPathName(Form("%s.C", sFunctionName.Data())), sNameFileIn.Data())));
  } else if(iSwitchMacro == 3) {
    sFunctionName = "CalculateFeedDown";
    reinterpret_cast<long>(gInterpreter->ProcessLine(Form(".x %s(\"%s\", %d, \"%s\", \"%s\", \"%s\")", gSystem->ExpandPathName(Form("%s.C", sFunctionName.Data())), sNameFileIn.Data(), V0Type, sNameFileEff.Data(), sPathXi.Data(), sFlag.Data())));
  } else if(iSwitchMacro == 4) {
    sFunctionName = "ScaleEfficiency";
    reinterpret_cast<long>(gInterpreter->ProcessLine(Form(".x %s(\"%s\", %d, \"%s\")", gSystem->ExpandPathName(Form("%s.C", sFunctionName.Data())), sNameFileIn.Data(), V0Type, sNameFileEff.Data())));
  } else if(iSwitchMacro == 5) {
    sFunctionName = "DrawCutVariations";
    reinterpret_cast<long>(gInterpreter->ProcessLine(Form(".x %s(%d)", gSystem->ExpandPathName(Form("%s.C", sFunctionName.Data())), V0Type)));
  } else if(iSwitchMacro == 6) {
    sFunctionName = "WeightMultiplicity";
    reinterpret_cast<long>(gInterpreter->ProcessLine(Form(".x %s(\"%s\", %d)", gSystem->ExpandPathName(Form("%s.C", sFunctionName.Data())), sNameFileIn.Data(), iMode)));
  }
}
