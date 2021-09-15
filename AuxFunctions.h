// Settings and auxiliary functions for the task AliAnalysisTaskV0sInJets
// Author: Vit Kucera (vit.kucera@cern.ch)

// upper edges of centrality bins
//const Int_t centBinRanges[] = {10, 20, 40, 60, 80}; // Vit Kucera, initial binning
//const Int_t centBinRanges[] = {10, 20, 40}; // Vit Kucera, restricted binning
//const Int_t centBinRanges[] = {5, 10, 20, 40, 60, 80}; // Iouri, LF analysis
//const Int_t centBinRanges[] = {10, 30, 50, 80}; // Alice Zimmermann
const Int_t centBinRanges[] = {10}; // only central
const Int_t iNCentBins = sizeof(centBinRanges) / sizeof(centBinRanges[0]);

// axis: pT of jets
Double_t dBinsPtJet[] = {0, 10, 20, 30, 40, 100}; // Vit Kucera
//Double_t dBinsPtJet[] = {0, 10, 20, 30, 40, 60, 100}; // Alice Zimmermann
//Double_t dBinsPtJet[] = {0, 20, 60, 100}; // Alice Zimmermann
//Double_t dBinsPtJet[] = {0, 100}; // all jets
const Int_t iNBinsPtJet = sizeof(dBinsPtJet) / sizeof(dBinsPtJet[0]) - 1;
const Int_t iNBinsPtJetInit = int((dBinsPtJet[iNBinsPtJet] - dBinsPtJet[0]) / 5.);

// axis: K0S invariant mass
Int_t iNMassBinsK0s = 200;
Float_t fMassK0sMin = 0.35;
Float_t fMassK0sMax = 0.65;
// axis: Lambda invariant mass
Int_t iNMassBinsLambda = 200;
Float_t fMassLambdaMin = 1.05;
Float_t fMassLambdaMax = 1.25;

// axis: pT of V0 inclusive
//Double_t dBinsPtV0All[] = {0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.6, 2.8, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 9, 10, 12};
//const Int_t iNBinsPtV0All = sizeof(dBinsPtV0All) / sizeof(dBinsPtV0All[0]) - 1; // number of bins of the rebinned V0 pt axis
//const Int_t iNBinsPtV0AllInit = int(10 * (dBinsPtV0All[iNBinsPtV0All] - dBinsPtV0All[0])); // number of bins of the V0 pt axis before rebinning
//Double_t dBinsPtV0AllLF[] = {0.0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0, 10.0, 12.0};
Double_t dBinsPtV0AllLF[] = {0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0, 10.0, 12.0};
const Int_t iNBinsPtV0AllLF = sizeof(dBinsPtV0AllLF) / sizeof(dBinsPtV0AllLF[0]) - 1; // number of bins of the rebinned V0 pt axis

// axis: pT of V0 in jets
//Double_t dBinsPtV0InJet[] = {0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.6, 2.8, 3, 4, 5, 8, 10, 15, 30};
Double_t dBinsPtV0InJet[] = {0, 0.5, 1, 1.5, 2, 3, 4, 5, 7, 10, 12};
//Double_t dBinsPtV0InJet[] = {0, 0.5, 1.5, 2.5, 3.5, 4.5, 6, 8, 12}; // Vojta
const Int_t iNBinsPtV0InJet = sizeof(dBinsPtV0InJet) / sizeof(dBinsPtV0InJet[0]) - 1; // number of bins of the rebinned V0 pt axis
const Int_t iNBinsPtV0InJetInit = int(10 * (dBinsPtV0InJet[iNBinsPtV0InJet] - dBinsPtV0InJet[0])); // number of bins of the V0 pt axis before rebinning

// axis: pT of V0 for determining the pt range to study variations of pt-integrate yield dependent on cut variations
Double_t dBinsPtV0Sys[] = {0, 2, 10, 12};
const Int_t iNBinsPtV0Sys = sizeof(dBinsPtV0Sys) / sizeof(dBinsPtV0Sys[0]) - 1; // number of bins of the rebinned V0 pt axis

// axis: eta of V0s for the efficiency correction
Float_t fRangeEtaV0Max = 0.8;
Double_t dBinsEtaInJet[] = { -fRangeEtaV0Max, -0.7, -0.65, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, fRangeEtaV0Max};
//Double_t dBinsEtaInJet[] = { -fRangeEtaV0Max, fRangeEtaV0Max};
//Double_t dBinsEtaInJet[] = { -fRangeEtaV0Max, 0, fRangeEtaV0Max};
//Double_t dBinsEtaInJet[] = { -fRangeEtaV0Max, -0.7, 0.7, fRangeEtaV0Max};
const Int_t iNBinsEtaInJet = sizeof(dBinsEtaInJet) / sizeof(dBinsEtaInJet[0]) - 1;
// axis: eta of V0s for the efficiency correction of inclusive spectra
Double_t dBinsEtaIncl[] = { -fRangeEtaV0Max, -0.75, -0.7, -0.65, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, fRangeEtaV0Max};
const Int_t iNBinsEtaIncl = sizeof(dBinsEtaIncl) / sizeof(dBinsEtaIncl[0]) - 1;

// Variations of selection cuts for systematics
Double_t dVarSysDCAV[] = {0.2}; // [cm] default: 0.1 (>0.1)
Double_t dVarSysDCAD[] = {0.45}; // [sigma] default: 1 (<1)
Double_t dVarSysCPAK[] = {0.9996}; // [] default: 0.998 (>0.998)
Double_t dVarSysCPAL[] = {0.9994}; // [] default: 0.998 (>0.998)
Double_t dVarSysTauK[] = {2.8}; // [tau] default: 5
Double_t dVarSysTauL[] = {2.8}; // [tau] default: 5
Double_t dVarSysRMin[] = {7.3}; // [cm] default: 5 (>5)
Double_t dVarSysRMax[] = {40}; // [cm] default: 100 (<100)
const Int_t iNVarSysDCAV = sizeof(dVarSysDCAV)/sizeof(dVarSysDCAV[0]);
const Int_t iNVarSysDCAD = sizeof(dVarSysDCAD)/sizeof(dVarSysDCAD[0]);
const Int_t iNVarSysCPA = sizeof(dVarSysCPAK)/sizeof(dVarSysCPAK[0]);
const Int_t iNVarSysTau = sizeof(dVarSysTauK)/sizeof(dVarSysTauK[0]);
const Int_t iNVarSysRMin = sizeof(dVarSysRMin)/sizeof(dVarSysRMin[0]);
const Int_t iNVarSysRMax = sizeof(dVarSysRMax)/sizeof(dVarSysRMax[0]);
Double_t* dArrayVarSysK[] = {dVarSysDCAV, dVarSysDCAD, dVarSysCPAK, dVarSysTauK, dVarSysRMin, dVarSysRMax};
Double_t* dArrayVarSysL[] = {dVarSysDCAV, dVarSysDCAD, dVarSysCPAL, dVarSysTauL, dVarSysRMin, dVarSysRMax};
Double_t** dArrayVarSys[] = {dArrayVarSysK, dArrayVarSysL};
Int_t iArrayNVarSys[] = {iNVarSysDCAV, iNVarSysDCAD, iNVarSysCPA, iNVarSysTau, iNVarSysRMin, iNVarSysRMax};
TString sArrayVarSysName[] = {"DCAV", "DCAD", "CPA", "Tau", "RMin", "RMax"};
Double_t dArrayValCutDefault[] = {0.1, 1, 0.998, 5, 5, 100};
enum cut {kDCAV, kDCAD, kCPA, kTau, kRMin, kRMax};
cut kArrayVarSysType[] = {kDCAV, kDCAD, kCPA, kTau, kRMin, kRMax};
const Int_t iNArrayVarSys = sizeof(dArrayVarSysK)/sizeof(dArrayVarSysK[0]);

Int_t GetCentralityBinIndex(Float_t centrality)
{
// returns index of the centrality bin corresponding to the provided value of centrality
  if(centrality < 0 || centrality > centBinRanges[iNCentBins - 1])
    return -1;
  for(Int_t i = 0; i < iNCentBins; i++)
  {
    if(centrality <= centBinRanges[i])
      return i;
  }
  return -1;
}

Int_t GetCentralityBinEdge(Int_t index)
{
// returns the upper edge of the centrality bin corresponding to the provided value of index
  if(index < 0 || index >= iNCentBins)
    return -1;
  return centBinRanges[index];
}

TString GetCentBinLabel(Int_t index)
{
  TString lowerEdge = ((index == 0) ? "0" : Form("%d", GetCentralityBinEdge(index - 1)));
  TString upperEdge = Form("%d", GetCentralityBinEdge(index));
  return Form("%s-%s %%", lowerEdge.Data(), upperEdge.Data());
}

const Int_t iNPtBins = 5; // V0 pt bins for tuning cuts
Float_t fPtBinRanges[iNPtBins] = {1, 2, 3, 6, 9};

Int_t GetPtBinIndex(Float_t pt)
{
  if(pt < 0 || pt > fPtBinRanges[iNPtBins - 1])
    return -1;
  for(Int_t i = 0; i < iNPtBins; i++)
  {
    if(pt <= fPtBinRanges[i])
      return i;
  }
  return -1;
}

Float_t MassPeakSigmaOld(Float_t pt, Int_t particle)
{
  switch(particle)
  {
    case 0: // K0S
      return 0.0044 + 0.0004 * (pt - 1.);
      break;
    case 1: // Lambda
      return 0.0023 + 0.00034 * (pt - 1.);
      break;
    default:
      return 0;
      break;
  }
}
