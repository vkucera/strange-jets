#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TMath.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"

#include "TStyle.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStopwatch.h"

#include "AuxROOTFunctions.h"
#include "AuxFunctions.h"
#include "BinCounterObject.cpp"

Bool_t bIsPbPb = 1;
TString sEnergy = "#sqrt{#it{s}_{NN}} = 2.76 TeV";
TString sSystem = "Pb#minusPb";
TLatex* labelSystem; // collision system
TLegend* legend;
// canvas size
Int_t iCanHeight = 600;
Int_t iCanWidth = 600;
TString sImageSuf = "png";
//TString sImageSuf = "pdf";

// ================================
Float_t fRadiusJet = 0.2; // R
Float_t fDistanceV0Jet = 0.2; // D
// ================================

Float_t fEtaV0Max = 0.7;
Float_t fEtaJetMax = fEtaV0Max - fDistanceV0Jet;
Float_t fNormalisationJetEta = 2.*fEtaJetMax;
Float_t fNormalisationV0Eta = 2.*fEtaV0Max;
Float_t fNormalisationV0Area = fNormalisationV0Eta * TMath::TwoPi();
//Float_t fRadiusConeExcluded = 0.8;
//Float_t fAreaConeExcluded = TMath::Pi()*fRadiusConeExcluded*fRadiusConeExcluded;
Float_t fAreaCone = TMath::Pi() * fDistanceV0Jet * fDistanceV0Jet;

TString V0Name[3] = {"K0s", "Lambda", "ALambda"};
TString V0Symbol[3] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}"};
TString V0LabelM[3] = {"#it{m}_{inv} K_{S}^{0} (GeV/#it{c}^{2})", "#it{m}_{inv} #Lambda (GeV/#it{c}^{2})", "#it{m}_{inv} #bar{#Lambda} (GeV/#it{c}^{2})"};

// initial parameters for fitting of invariant-mass distribution
Float_t fMassV0[3] = {0.497614, 1.115680, 1.115680}; // [GeV/c^2] initial centre of Gaussian for global fit
Float_t fSigmaV0[3] = {0.006, 0.0017, 0.0017}; // [GeV/c^2] initial sigma of Gaussian for global fit
// ranges of invariant mass for the bin counting
// fixed regions [GeV/c^2]
Float_t fRangeBgOut[3][2] = {{0.38, 0.65}, {1.1, 1.155}, {1.1, 1.155}}; // dBgOutMin, dBgOutMax, maximum range for calculations
Float_t fRangeBgIn[3][2] = {{0.43, 0.57}, {1.105, 1.13}, {1.105, 1.13}}; // dBgInMin, dBgInMax
Float_t fRangeSignal[3][2] = {{0.43, 0.57}, {1.105, 1.13}, {1.105, 1.13}}; // dSigMin, dSigMax
// old settings
//Float_t fRangeBgOut[3][2] = {{0.4, 0.6}, {1.095, 1.16}, {1.095, 1.16}}; // dBgOutMin, dBgOutMax, maximum range for calculations
//Float_t fRangeBgIn[3][2] = {{0.44, 0.55}, {1.105, 1.13}, {1.105, 1.13}}; // dBgInMin, dBgInMax
//Float_t fRangeSignal[3][2] = {{0.44, 0.55}, {1.105, 1.13}, {1.105, 1.13}}; // dSigMin, dSigMax

// variations

//// 1 K0S
//Float_t fRangeBgOut[3][2] = {{0.38, 0.65}, {1.1, 1.155}, {1.1, 1.155}}; // dBgOutMin, dBgOutMax, maximum range for calculations
//Float_t fRangeBgIn[3][2] = {{0.45, 0.55}, {1.105, 1.13}, {1.105, 1.13}}; // dBgInMin, dBgInMax
//Float_t fRangeSignal[3][2] = {{0.45, 0.55}, {1.105, 1.13}, {1.105, 1.13}}; // dSigMin, dSigMax

//// 2 Lambda
//Float_t fRangeBgOut[3][2] = {{0.38, 0.65}, {1.1, 1.155}, {1.1, 1.155}}; // dBgOutMin, dBgOutMax, maximum range for calculations
//Float_t fRangeBgIn[3][2] = {{0.43, 0.57}, {1.110, 1.125}, {1.110, 1.125}}; // dBgInMin, dBgInMax
//Float_t fRangeSignal[3][2] = {{0.43, 0.57}, {1.110, 1.125}, {1.110, 1.125}}; // dSigMin, dSigMax

//// 3 K0S & Lambda
//Float_t fRangeBgOut[3][2] = {{0.38, 0.65}, {1.1, 1.155}, {1.1, 1.155}}; // dBgOutMin, dBgOutMax, maximum range for calculations
//Float_t fRangeBgIn[3][2] = {{0.45, 0.55}, {1.110, 1.125}, {1.110, 1.125}}; // dBgInMin, dBgInMax
//Float_t fRangeSignal[3][2] = {{0.45, 0.55}, {1.110, 1.125}, {1.110, 1.125}}; // dSigMin, dSigMax


// sigma regions:
Double_t dNSigmaBC[3][3] = {{5, 5.5, 20}, {5.5, 5.5, 20}, {3.5, 4, 20}}; // dSigmaSig, dSigmaBgIn, dSigmaBgOut, LF
//Double_t dNSigmaBC[2][3] = {{4,5,16},{4,5,9}}; // Alice

Float_t fFitSigma[3] = {5e-3, 2.5e-3, 2.5e-3};
//      Float_t fFitSigma[2] = {0.004, 0.0015};
// range of signal mass peak, should correspond to analysis window
Float_t fSignalWindow[2][2] = {{fMassV0[0] - 3 * fFitSigma[0], fMassV0[0] + 3 * fFitSigma[0]}, {fMassV0[1] - 3 * fFitSigma[1], fMassV0[1] + 3 * fFitSigma[1]}};

// V0 selection steps
static const Int_t iNCategV0 = 17;
TString categV0[iNCategV0] = {"all"/*0*/, "mass range"/*1*/, "rec. method"/*2*/, "tracks TPC"/*3*/, "track pt"/*4*/, "DCA prim v"/*5*/, "DCA daughters"/*6*/, "CPA"/*7*/, "volume"/*8*/, "track #it{#eta}"/*9*/, "V0 #it{y} & #it{#eta}"/*10*/, "lifetime"/*11*/, "PID"/*12*/, "Arm.-Pod."/*13*/, "inclusive"/*14*/, "in jet event"/*15*/, "in jet"/*16*/};

Double_t dEpsilon = 1e-4;

// Global ranges
// centrality bins
Int_t iCentMin = 0;
//Int_t iCentMax = iNCentBins-1;
Int_t iCentMax = 0;
// jet pT bins
Int_t iJetMin = 1;
//Int_t iJetMax = iNBinsPtJet-1;
Int_t iJetMax = 2;

// axis ranges:
// pT inclusive
Float_t fPtAllXMin = 0;
Float_t fPtAllXMax = 12;
Float_t fPtAllYMin = 1e-5;
Float_t fPtAllYMax = 5e2;
Int_t iBinPtInclFirst = 2;
Int_t iBinPtInclLast = iNBinsPtV0AllLF;
// pT in jets
Float_t fPtInJetsXMin = 0.; // minimum y value
Float_t fPtInJetsXMax = 10.; // maximum y value
Float_t fPtInJetsYMin = 1e-4; // minimum y value
Float_t fPtInJetsYMax = 10; // maximum y value
Int_t iBinPtInJetsFirst = 5;
Int_t iBinPtInJetsLast = iNBinsPtV0InJet - 1;
// pT for systematics
Int_t iBinPtSysFirst = 2;
Int_t iBinPtSysLast = iBinPtSysFirst;
// jet pT
Float_t fPtJetXMin = 0.; // minimum y value
Float_t fPtJetXMax = 100.; // maximum y value
Float_t fPtJetYMin = 1e-8; // minimum y value
Float_t fPtJetYMax = 1e0; // maximum y value

// BC settings
Bool_t bFixedBC = 1; // fixed ranges in bin counting
Int_t iDegreePolSideBands = 2;
Int_t iNRebin = 0; // number of bins to be merged
Bool_t bVerboseBC = 1;
TString sOptionFitSB = "SLRIN"; // option M causes crash of Fit for pol0 in side bands, E crashes with fixed side bands

//=================================================
// Global switches
Bool_t bYieldsOnly = 0; // don't normalize spectra (event number, bin width, jet area)

// Corrections
Bool_t bCorrAll = 0; // switch on all corrections (iCorrEff, iCorrFD, iUEsubtraction)

Int_t iCorrEff = 0; // apply V0 efficiency correction
//{
Int_t iCorrEffInclEta = 0; // apply efficiency reweighting in eta for inclusive particles
Int_t iCorrEffInclusive = 0; // use inclusive rebinned V0 for the efficiency correction of V0s in UE, JC, jets
Int_t iCorrEffPur = 1; // apply purity scaling for the extraction of signal V0s in JC and UE in pt-eta bins for the efficiency weighting
//}
Int_t iCorrFD = 0; // apply feed-down correction
//{
Int_t iCorrFDSel = 1; // choice of FD method: 0 - PYTHIA jets, 1 - calculated inclusive (default, needed for bulk), 2 - LF 2010 0-10 %
//}
Int_t iUEsubtraction = 0; // subtract bulk from spectra in jet cones
//{
Int_t iUESub = 0; // choice of UE method: 0 NJ (default), 1 OC, 2 RC, 3 MCC, 4 PC
//}

Bool_t bOpenPtBins = 1; // open range of jet pT
Bool_t bNormPerArea = 1; // normalize per area unit instead of per jet
Bool_t bEffEtaWindow = 1; // calculate pT-eta efficiency in eta windows instead of eta slices
Bool_t bZeroGenErr = 1; // set all stat errors of generated spectra to zero
Bool_t bInclusiveDensity = 0; // normalize inclusive spectra per area -> density
//=================================================

Int_t iStd = 1;
//{
Int_t iEvents = 0; // event selection
//{
Int_t iEventsBias = 0; // event selection
//}
Int_t iCandVsCut = 0; // candidates ater cuts
Int_t iSignalVsCent = 0; // signal in cent bins
Int_t iSignalVsCut = 0; // signal after cuts
Int_t iPtInclusive = 1; // inclusive spectra
//{
Int_t iPtInclusiveJetBins = 1; // enable inclusive spectra with the jet pt binning (needed for inclusive FD with jet pt binning)
Int_t iPtInclusiveSys = 0; // enable extraction of signal integrated over pt
//}
Int_t iPtEtaInclusive = 0; // inclusive spectra as 2D histogram: pt_V0; eta_V0 (needed for the 2D signal extraction in jets with purity)
//{
Int_t iSwitchPtV0 = 0; // pt_V0 binning: 0 inclusive V0s, 1 V0s in jets (used in iEffEtaPtInclusive, needed for 2D efficiency with jet pt binning)
//}
Int_t iPtInJets = 0; // spectra in jets
Int_t iCorrelations = 0; // angular V0-jet correlations
//}
Int_t iMC = 0; // simulated data (efficiency calculation)
//{
Int_t iEffPtInclusive = 1; // efficiency of inclusive V0s as 1D (pt_V0) histogram
Int_t iEffPtInJets = 0; // efficiency of V0s in jet cones as 1D (pt_V0) histogram and make ratio w.r.t. inclusive V0s
Int_t iEffEtaPtInclusive = 1; // efficiency of inclusive V0s as eta_V0 in pt_V0 slices, store it as 2D histogram: pt_V0; eta_V0, needed for efficiency scaling (required always)
Int_t iEffEtaPtDaughters = 0;
Int_t iEffPtEtaInclusive = 0; // efficiency of inclusive V0s as 1D (pt_V0) histogram in eta_V0 slices
Int_t iEffPtEtaInJets = 0; // efficiency of V0s in jet cones as 1D (pt_V0) histogram in eta_V0 slices, comparison with inclusive V0s
Int_t iMassResol = 0;
//}
Int_t iJetSpectrum = 0; // spectra of jets selected for V0 analysis
//{
Int_t iPtCorrel = 0; // pt correlations leading track, trigger track, jet
//}
//Int_t iCuts = 0;
//{
//Int_t iTestCut = 0;
//}

Double_t GetMeanFeedDownLambdaLF(Double_t dPtMin, Double_t dPtMax);
Double_t GetMeanFeedDownLambdaPYTHIA(Double_t dPtMin, Double_t dPtMax);
Double_t GetMeanValue(TH1* his);
TH1D* GetScaledEfficiency(TH2* hisPtEtaMeasured, TH2* hisPtEtaEffBase, Int_t iSwitchPt = 1);
void DrawResultsLoopSys(TString sNameFileIn, Int_t V0Type = 0, TString sNameFileEff = "", TString sNameFileFD = "", TString sFlag = "", Int_t iMode = 0);

void DrawResults(TString sNameFileIn, Int_t V0Type = 0, TString sNameFileEff = "", TString sNameFileFD = "", TString sFlag = "" /* modifier of the branch name*/, Int_t iMode = 0)
{
  printf("DrawResults: Start\n");
  printf("Flag %s\n", sFlag.Data());

  // counting duration of analysis
  TStopwatch timer;
  timer.Start();

  gStyle->SetCanvasColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetOptStat(0);

  // modes: 0 uncorrected, 1 simulated, 2 corrected
  switch (iMode)
  {
    case 0:
      iStd = 1;
      iMC = 0;
      break;
    case 1:
      iMC = 1;
      break;
    case 2:
      iStd = 1;
      iMC = 0;
      iCorrEff = 1;
      iCorrFD = 1;
      iUEsubtraction = 1;
      break;
    default:
      printf("Error: Wrong mode.\n");
      return;
      break;
  }


  TString branchName = ""; // name of the jet task (branch)

//  TString sYear = "2010";
//  Int_t iFilterBit = 0; // JETAN
//  Int_t iBgSub = 2; // JETAN, type of rhosubtraction
//  iFilterBit = ((sYear == "2011") ? 768 : 272);
//  iFilterBit = 0;

//  branchName = Form("clustersAOD_ANTIKT%02d_B1_Filter%05d_Cut00150_T00_J02",int(10*fRadiusJet),iFilterBit);
//  branchName = Form("clustersAOD_ANTIKT%02d_B%d_Filter%05d_Cut00150_Skip00", int(10 * fRadiusJet), iBgSub, iFilterBit);
//  branchName = Form("Jet_AKTChargedR%03d_PicoTracks_pT0150_pt_scheme_D%02d", int(100 * fRadiusJet), int(10 * fDistanceV0Jet));
  branchName = Form("Jet_AKTChargedR%03d_tracks_pT0150_pt_scheme", int(100 * fRadiusJet));

//  sFlag = "C0005";
//  sFlag = "C0510";
//  sFlag = "Lead10";
//  sFlag = "CPA_T";
//  sFlag = "CPA_TT";
//  sFlag = "Tau_T";
//  sFlag = "Tau_T";
//  sFlag = "DCAD_T";
//  sFlag = "DCAD_TT";
//  sFlag = "Ionut";

  if(bCorrAll)
  {
    iCorrEff = 1;
    iCorrFD = 1;
    iUEsubtraction = 1;

    iPtInclusive = 1;
    iPtInJets = 1;
    iJetSpectrum = 1;

    iPtInclusiveJetBins = 0;
  }

  if(iMC)
  {
    branchName = "";
    iStd = 0;
    iCorrEff = 0;
    iCorrFD = 0;
    if(iEffPtInJets)
      iEffPtInclusive = 1;
    if(iEffPtEtaInJets)
      iEffPtEtaInclusive = 1;
    iJetSpectrum = 0;
  }
  if(sFlag.Length())
  {
    if(branchName.Length())
      branchName += "_";
    branchName += Form("%s", sFlag.Data());
  }
//  branchName = "";

  if(V0Type == 0)
  {
//    bFixedBC = 1;
    iDegreePolSideBands = 2;
    iBinPtInclFirst = 2;
//    iBinPtInclFirst = 1;
    iCorrFD = 0;
  }
  if(V0Type == 1 || V0Type == 2)
  {
//    bFixedBC = 1;
    iDegreePolSideBands = 3;
    iBinPtInclFirst = 5;
//    iBinPtInclFirst = 1;
  }
  if(iPtInJets && iCorrEff && !iCorrEffInclusive)
  {
    iPtEtaInclusive = 1;
    iSwitchPtV0 = 1;
  }

  TString sLabelCollisionText = Form("#splitline{%s,}{%s}", sSystem.Data(), sEnergy.Data());
  TLatex* labelCollision = new TLatex();
  labelCollision->SetTextFont(42);
  labelCollision->SetTextSize(0.04);
  labelCollision->SetTextAlign(23);

  // bin counting ranges
  // fixed range
  Double_t dBgOutMin, dBgOutMax, dBgInMin, dBgInMax, dSigMin, dSigMax;
  dBgOutMin = fRangeBgOut[V0Type][0];
  dBgOutMax = fRangeBgOut[V0Type][1];
  dBgInMin = fRangeBgIn[V0Type][0];
  dBgInMax = fRangeBgIn[V0Type][1];
  dSigMin = fRangeSignal[V0Type][0];
  dSigMax = fRangeSignal[V0Type][1];
  // multiples of sigma
  Double_t dSigmaSig = dNSigmaBC[V0Type][0];
  Double_t dSigmaBgIn = dNSigmaBC[V0Type][1];
  Double_t dSigmaBgOut = dNSigmaBC[V0Type][2];

// function for fitting background
  TF1* funBg = new TF1("bg", "pol3", dBgOutMin, dBgOutMax);
// function for fitting shape of signal peak in range of 3 sigma
//  TF1* funSignal = new TF1("signal","gaus",fSignalWindow[V0Type][0],fSignalWindow[V0Type][1]);
  // array of 3 histograms to store results of bin counting (all, bg, sig)
//  TH1D* hisArrayBC[3];

  printf("Loading file %s ", sNameFileIn.Data());
  TFile* fileHisto = 0;
  fileHisto = new TFile(sNameFileIn.Data(), "READ");
  if(fileHisto->IsZombie())
  {
    printf("failed\nError: Cannot load file %s\n", sNameFileIn.Data());
    return;
  }
  printf("OK\n");

  TFile* fileEff = 0;
  if(iStd && iCorrEff)
  {
    if(!sNameFileEff.Length())
    {
      printf("Error: Name of efficiency file not specified!\n");
      return;
    }
    printf("Loading file %s ", sNameFileEff.Data());
    fileEff = new TFile(sNameFileEff.Data(), "READ");
    if(fileEff->IsZombie())
    {
      printf("failed\nError: Cannot load file %s\n", sNameFileEff.Data());
      return;
    }
    printf("OK\n");
  }

  TFile* fileFD = 0;
  if(iStd && iPtInJets && iCorrFD && iCorrFDSel != 2)
  {
    if(!sNameFileFD.Length())
    {
      printf("Error: Name of feed-down file not specified!\n");
      return;
    }
    printf("Loading file %s ", sNameFileFD.Data());
    fileFD = new TFile(sNameFileFD.Data(), "READ");
    if(fileFD->IsZombie())
    {
      printf("failed\nError: Cannot load file %s\n", sNameFileFD.Data());
      return;
    }
    printf("OK\n");
  }

  TString sPwd = gSystem->pwd();
  TString sDirOutName = sPwd;
//  if(sFlag.Length())
//    sDirOutName += "/" + sFlag;
  gSystem->mkdir(sDirOutName.Data(), 1);
  TString kOutputFileName = Form("%s/Output%s.root", sDirOutName.Data(), V0Name[V0Type].Data());
//  TString kOutputFileName = Form("Output%s.root", V0Name[V0Type].Data());
  printf("Creating new output file %s ", kOutputFileName.Data());
  TFile* fileOutput = new TFile(kOutputFileName.Data(), "RECREATE");
//  TFile* fileOutput = new TFile(kOutputFileName.Data(),"UPDATE");
  if(fileOutput->IsZombie())
  {
    printf("failed\nError: Cannot create output file\n");
    return;
  }
  printf("OK\n");

  TString sNameDirPt = "Spectra";
  TDirectory* dirOutSpectra = fileOutput->mkdir(sNameDirPt.Data());

  printf("Executing block for drawing V0s\n");
  TString dirNameV0 = "V0";
  TString listName = "V0histo";
  if(branchName.Length())
  {
    dirNameV0 += Form("_%s", branchName.Data());
    listName += Form("_%s", branchName.Data());
  }

  TString listStdName = listName + "_Std";
  TList* listStd = GetList(fileHisto, dirNameV0, listStdName.Data());
  if(!listStd)
  {
    printf("Error: Couldn't get list %s\n", listStdName.Data());
    return;
  }
//  TString listCutsName = listName+"_Cuts";
//  TList* listCuts = GetList(fileHisto,dirNameV0,listCutsName.Data());
//  if (!listCuts)
//    {
//      printf("Error: Couldn't get list %s\n",listCutsName.Data());
//      return;
//    }
  TString listMCName = listName + "_MC";
  TList* listMC = 0;
  if(iMC)
  {
    listMC = GetList(fileHisto, dirNameV0, listMCName.Data());
    if(!listMC)
    {
      printf("Error: Couldn't get list %s\n", listMCName.Data());
      return;
    }
  }

  printf("Statistics: Events (cent bins):");
  for(Int_t iCent = 0; iCent < iNCentBins; iCent++)
    printf(" %d", centBinRanges[iCent]);
  printf("\n");

  // loading events per centrality bin
  TString hisEventCentName = "fh1EventCent";
  TH1D* hisEventCent = GetHistogram1D(listStd, hisEventCentName.Data());
  if(!hisEventCent)
    return;
  CheckHistogram(hisEventCent);
  Double_t eventsCent[iNCentBins];
  Double_t eventTotal = 0;
  printf("Loading numbers of events per cent bin\nStatistics: Events:");
  for(Int_t iCent = 0; iCent < iNCentBins; iCent++)
  {
    eventsCent[iCent] = hisEventCent->GetBinContent(iCent + 1);
    printf(" %d", int(eventsCent[iCent]));
    eventTotal += eventsCent[iCent];
  }
  printf("; total: %d\n", int(eventTotal));

  // loading jet pt spectrum for normalisation in pt spectra
  TString hisJetPtName = "fh1PtJet_%d";
  TH1D* hisJetPt[iNCentBins];
  if((iStd && iPtInJets) || (iMC && iEffPtInJets))
    for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
    {
      hisJetPt[iCent] = GetHistogram1D(listStd, Form(hisJetPtName.Data(), iCent));
      if(!hisJetPt[iCent])
        return;
      CheckHistogram(hisJetPt[iCent]);
    }

  TH1D* hisEffPtInclExt[iNCentBins];
  TH1D* hisEffPtInJetsExt[iNCentBins];
  TH1D* hisEffInPCExt[iNCentBins];
  TH2D* hisEffPtEtaInclBase[iNCentBins];
  TH2D* hisEffPtEtaInclBaseInclBins[iNCentBins];
  TH2D* hisEffPtEtaInclBaseJetBins[iNCentBins];
  TString sNameHisEffPtIncl = "fh1EffPtIncl%s_%d";
  TString sNameHisEffPtInJets = "fh1EffPtInJets%s_%d";
  TString sNameHisEffPtInJetsTrue = "fh1EffPtInJets%s_C%d-J%d";
  TString sNameHisEffPtInPC = "fh1EffPtInPC%s_%d";
  TString sNameHisEffPtEtaIncl = "fh1EffPtEtaIncl%s_C%d_E%d";
  TString sNameHisEffPtEtaIncl2D = "fh2EffPtEtaIncl%s_C%d";
  TString sNameHisEffPtEtaInJets = "fh1EffPtEtaInJets%s_C%d_J%d_E%d";
  TString sNameHisGenPtEtaIncl = "fh1GenPtEtaIncl%s_C%d_E%d";
  TString sNameHisGenPtEtaInJets = "fh1GenPtEtaInJets%s_C%d_J%d_E%d";
  TDirectoryFile* dirEff;
  if(iStd && iCorrEff && (iPtInclusive || iPtEtaInclusive || iPtInJets || iCorrelations))
  {
    printf("Loading directory %s ", sNameDirPt.Data());
    dirEff = (TDirectoryFile*)fileEff->Get(sNameDirPt.Data());
    if(!dirEff)
    {
      printf("failed\nError: Loading dir\n");
      return;
    }
    printf("OK\n");

    if(iPtInJets || iCorrelations)
    {
      for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
      {
        // inclusive eff pT-eta for the eta scaling
        TString sNameHisEffPtEtaInclBase = Form(sNameHisEffPtEtaIncl2D.Data(), V0Name[V0Type].Data(), iCent);
        TString sNameHisEffPtEtaInclBaseJet = Form("%s-JetBins", sNameHisEffPtEtaInclBase.Data());
        printf("Loading histogram %s ", sNameHisEffPtEtaInclBaseJet.Data());
        hisEffPtEtaInclBase[iCent] = (TH2D*)dirEff->Get(sNameHisEffPtEtaInclBaseJet.Data());
        if(!hisEffPtEtaInclBase[iCent])
        {
          printf("failed. (Error)\n");
          return;
        }
        printf("OK\n");
      }
    }

    if(iPtEtaInclusive && iCorrEffInclEta)
    {
      for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
      {
        TString sNameHisEffPtEtaInclBase = Form(sNameHisEffPtEtaIncl2D.Data(), V0Name[V0Type].Data(), iCent);
        if(iSwitchPtV0 == 0) // inclusive pt binning
        {
          TString sNameHisEffPtEtaInclBaseIncl = Form("%s-InclBins", sNameHisEffPtEtaInclBase.Data());
          printf("Loading histogram %s ", sNameHisEffPtEtaInclBaseIncl.Data());
          hisEffPtEtaInclBaseInclBins[iCent] = (TH2D*)dirEff->Get(sNameHisEffPtEtaInclBaseIncl.Data());
          if(!hisEffPtEtaInclBaseInclBins[iCent])
          {
            printf("failed. (Error)\n");
            return;
          }
          printf("OK\n");
        }
        else if(iSwitchPtV0 == 1) // in-jet pt binning
        {
          TString sNameHisEffPtEtaInclBaseJet = Form("%s-JetBins", sNameHisEffPtEtaInclBase.Data());
          printf("Loading histogram %s ", sNameHisEffPtEtaInclBaseJet.Data());
          hisEffPtEtaInclBaseJetBins[iCent] = (TH2D*)dirEff->Get(sNameHisEffPtEtaInclBaseJet.Data());
          if(!hisEffPtEtaInclBaseJetBins[iCent])
          {
            printf("failed. (Error)\n");
            return;
          }
          printf("OK\n");
        }
      }
    }

    if(iPtInclusive)
    {
      for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
      {
        TString sNameHisEffPtInclFin = Form(sNameHisEffPtIncl.Data(), V0Name[V0Type].Data(), iCent);
        printf("Loading histogram %s ", sNameHisEffPtInclFin.Data());
        hisEffPtInclExt[iCent] = (TH1D*)dirEff->Get(sNameHisEffPtInclFin.Data());
        if(!hisEffPtInclExt[iCent])
        {
          printf("failed. (Error)\n");
          return;
        }
        printf("OK\n");
      }
    }

    if((iPtInJets && iCorrEffInclusive) || (iPtInclusive && iPtInclusiveJetBins))
    {
      for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
      {
        TString sNameHisEffPtInJetsFin = Form(sNameHisEffPtInJets.Data(), V0Name[V0Type].Data(), iCent);
        printf("Loading histogram %s", sNameHisEffPtInJetsFin.Data());
        hisEffPtInJetsExt[iCent] = (TH1D*)dirEff->Get(sNameHisEffPtInJetsFin.Data());
        if(!hisEffPtInJetsExt[iCent])
        {
          printf("failed. (Error)\n");
          return;
        }
        printf("OK\n");
//                    TString sNameHisEffPtInPCFin = Form(sNameHisEffPtInPC.Data(),V0Name[V0Type].Data(),iCent);
//                    printf("Loading histogram %s\n",sNameHisEffPtInPCFin.Data());
//                    hisEffInPCExt[iCent] = (TH1D*)dirEff->Get(sNameHisEffPtInPCFin.Data());
//                    if (!hisEffInPCExt[iCent])
//                      return;
      }
    }
  }

  // get feed-down data
  TDirectoryFile* dirFD;
  TH1D* hisFDInclJetBins[iNCentBins];
  if(iStd && iPtInJets && iCorrFD && iCorrFDSel != 2)
  {
    // get feed-down directory
    printf("Loading directory %s ", sNameDirPt.Data());
    dirFD = (TDirectoryFile*)fileFD->Get(sNameDirPt.Data());
    if(!dirFD)
    {
      printf("failed (Error)\n");
      return;
    }
    printf("OK\n");

    // get feed down histograms
    for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
    {
      TString sAppendix = "-JetBins";
      TString sNameHisFDInclJetBins = Form("fh1FDFractionIncl_C%d%s", iCent, sAppendix.Data());
      printf("Loading histogram %s ", sNameHisFDInclJetBins.Data());
      hisFDInclJetBins[iCent] = (TH1D*)dirFD->Get(sNameHisFDInclJetBins.Data());
      if(!hisFDInclJetBins[iCent])
      {
        printf("failed. (Error)\n");
        return;
      }
      printf("OK\n");
    }
  }

  gSystem->cd(sDirOutName.Data()); // has to be after all reading from input files

  if(iStd)
  {
    if(iEvents)
    {
      TString sNameHisCentTot = "fh1EventCent";
      TH1D* hisCentrality;
      TCanvas* canCentrality = new TCanvas("canCentrality", "", iCanWidth, iCanHeight);
      hisCentrality = GetHistogram1D(listStd, sNameHisCentTot.Data());
      if(!hisCentrality)
        return;
      fileOutput->cd();
//            hisCentrality->Sumw2(kFALSE);
      hisCentrality->Write();
      canCentrality->cd();
      hisCentrality->Draw("");
      canCentrality->SaveAs(Form("canCentrality.%s", sImageSuf.Data()));
      delete canCentrality;

      TString sNameHisCentTot2 = "fh1EventCent2";
      TString sNameHisCentTot2Jets = "fh1EventCent2Jets";
      TString sNameHisCentTot2NoJets = "fh1EventCent2NoJets";
      TH1D* hisCentrality2;
      TH1D* hisCentrality2Jets;
      TH1D* hisCentrality2NoJets;
      TCanvas* canCentrality2 = new TCanvas("canCentrality2", "", iCanWidth, iCanHeight);
      hisCentrality2 = GetHistogram1D(listStd, sNameHisCentTot2.Data());
      if(!hisCentrality2)
        return;
      fileOutput->cd();
//            hisCentrality2->Sumw2(0);
      hisCentrality2->Write();
      canCentrality2->cd();
      hisCentrality2->Draw("");
      canCentrality2->SaveAs(Form("canCentrality2.%s", sImageSuf.Data()));
      delete canCentrality2;

      if(iEventsBias)
      {
        TCanvas* canCentrality2Bias = new TCanvas("canCentrality2Bias", "", iCanWidth, iCanHeight);
        TCanvas* canCentrality2BiasRatio = new TCanvas("canCentrality2BiasRatio", "", iCanWidth, iCanHeight);
        TMultiGraph* mgrCentBias = new TMultiGraph();
        TMultiGraph* mgrCentBiasRatio = new TMultiGraph();
        hisCentrality2Jets = GetHistogram1D(listStd, sNameHisCentTot2Jets.Data());
        if(!hisCentrality2Jets)
          return;
        hisCentrality2NoJets = GetHistogram1D(listStd, sNameHisCentTot2NoJets.Data());
        if(!hisCentrality2NoJets)
          return;
        hisCentrality2->Scale(1. / hisCentrality2->Integral());
        hisCentrality2Jets->Scale(1. / hisCentrality2Jets->Integral());
        hisCentrality2NoJets->Scale(1. / hisCentrality2NoJets->Integral());
//          hisCentrality2->Sumw2(0);
        TGraphErrors* grCentAll = MakeGraphErrors(hisCentrality2, "all events", iMyColors[0], iMyMarkersFull[0]);
        TGraphErrors* grCentJets = MakeGraphErrors(hisCentrality2Jets, "jet events", iMyColors[1], iMyMarkersFull[0]);
        TGraphErrors* grCentNoJets = MakeGraphErrors(hisCentrality2NoJets, "no-jet events", iMyColors[2], iMyMarkersFull[0]);
        mgrCentBias->Add(grCentAll);
        mgrCentBias->Add(grCentJets);
        mgrCentBias->Add(grCentNoJets);

        TH1D* hisCentRatioJNJ = DivideHistograms1D(hisCentrality2Jets, hisCentrality2NoJets);
        TGraphErrors* grCentRatioJNJ = MakeGraphErrors(hisCentRatioJNJ, "jet events/no-jet events", iMyColors[3], iMyMarkersFull[0]);
        mgrCentBiasRatio->Add(grCentRatioJNJ);
        TH1D* hisCentRatioJAll = DivideHistograms1D(hisCentrality2Jets, hisCentrality2);
        TGraphErrors* grCentRatioJAll = MakeGraphErrors(hisCentRatioJAll, "jet events/all events", iMyColors[1], iMyMarkersFull[0]);
        mgrCentBiasRatio->Add(grCentRatioJAll);
        TH1D* hisCentRatioNJAll = DivideHistograms1D(hisCentrality2NoJets, hisCentrality2);
        TGraphErrors* grCentRatioNJAll = MakeGraphErrors(hisCentRatioNJAll, "no-jet events/all events", iMyColors[2], iMyMarkersFull[0]);
        mgrCentBiasRatio->Add(grCentRatioNJAll);

        canCentrality2Bias->cd();
        canCentrality2Bias->SetLeftMargin(0.15);
        mgrCentBias->SetTitle("centrality distribution;centrality (%);arb. unit");
        mgrCentBias->Draw("AP0");
        mgrCentBias->GetXaxis()->SetLimits(0, 10);
        mgrCentBias->GetYaxis()->SetRangeUser(0.06, 0.12);
        mgrCentBias->GetYaxis()->SetTitleOffset(2);
        legend = canCentrality2Bias->BuildLegend(0.25, 0.15, 0.75, 0.4);
        SetLegend(legend);
        canCentrality2Bias->SaveAs(Form("canCentrality2Bias.%s", sImageSuf.Data()));
        delete canCentrality2Bias;
        delete mgrCentBias;

        TLine* lineOneCent = new TLine(0, 1, 10, 1);
        canCentrality2BiasRatio->cd();
        canCentrality2BiasRatio->SetLeftMargin(0.15);
        mgrCentBiasRatio->SetTitle("centrality distribution, ratios;centrality (%);arb. unit");
        mgrCentBiasRatio->Draw("AP0");
        mgrCentBiasRatio->GetXaxis()->SetLimits(0, 10);
        mgrCentBiasRatio->GetYaxis()->SetRangeUser(0.8, 1.1);
        mgrCentBiasRatio->GetYaxis()->SetTitleOffset(2);
        legend = canCentrality2BiasRatio->BuildLegend(0.25, 0.15, 0.75, 0.4);
        SetLegend(legend);
        lineOneCent->Draw();
        canCentrality2BiasRatio->SaveAs(Form("canCentrality2BiasRatio.%s", sImageSuf.Data()));
        delete canCentrality2BiasRatio;
        delete mgrCentBiasRatio;
        delete lineOneCent;
      }
    }

    if(iCandVsCut)
    {
      TH1D* hisV0CounterCent;
      TString hisV0CounterCentName = Form("fh1V0CounterCent%s_%%d", V0Name[V0Type].Data());
      TGraphErrors* grV0CounterCent;
      TCanvas* canV0CounterAbs = new TCanvas("canV0CounterAbs", "", iCanWidth, iCanHeight);
      TMultiGraph* mgrV0CounterAbs = new TMultiGraph();
      TCanvas* canV0CounterRel = new TCanvas("canV0CounterRel", "", iCanWidth, iCanHeight);
      TMultiGraph* mgrV0CounterRel = new TMultiGraph();
      Float_t fYMax = 1;
      for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
      {
        hisV0CounterCent = GetHistogram1D(listStd, Form(hisV0CounterCentName.Data(), iCent));
        if(!hisV0CounterCent)
          return;
        printf("Statistics: Candidates (c. %d): inclusive: %d, in jet events: %d, in jets: %d\n", iCent, int(hisV0CounterCent->GetBinContent(15)), int(hisV0CounterCent->GetBinContent(16)), int(hisV0CounterCent->GetBinContent(17)));
        if(!eventsCent[iCent])
          continue;
        if(!hisV0CounterCent->GetBinContent(1))
          continue;
        if(iCent == 0)
          fYMax = hisV0CounterCent->GetBinContent(1);
        grV0CounterCent = MakeGraphErrors(hisV0CounterCent, GetCentBinLabel(iCent).Data(), iMyColors[iCent]);
        mgrV0CounterAbs->Add(grV0CounterCent);
        hisV0CounterCent->Scale(1. / hisV0CounterCent->GetBinContent(1));
        grV0CounterCent = MakeGraphErrors(hisV0CounterCent, GetCentBinLabel(iCent).Data(), iMyColors[iCent]);
        mgrV0CounterRel->Add(grV0CounterCent);
      }
      mgrV0CounterAbs->SetTitle(Form("Number of %s candidates after cuts;;counts", V0Symbol[V0Type].Data()));
      mgrV0CounterRel->SetTitle(Form("Number of %s candidates after cuts;;relative counts", V0Symbol[V0Type].Data()));
      canV0CounterAbs->cd();
      canV0CounterAbs->SetLeftMargin(0.15);
      mgrV0CounterAbs->Draw("AP0");
      mgrV0CounterAbs->GetXaxis()->Set(iNCategV0, 0, iNCategV0);
      for(Int_t j = 0; j < iNCategV0; j++)
        mgrV0CounterAbs->GetXaxis()->SetBinLabel(j + 1, categV0[j].Data());
      mgrV0CounterAbs->GetYaxis()->SetTitleOffset(2);
      mgrV0CounterAbs->SetMinimum(1);
      mgrV0CounterAbs->SetMaximum(2 * fYMax);
      canV0CounterAbs->SetLogy();
      legend = canV0CounterAbs->BuildLegend(0.2, 0.15, 0.35, 0.3);
      SetLegend(legend);
      canV0CounterAbs->SaveAs(Form("canCounterAbs%s.%s", V0Name[V0Type].Data(), sImageSuf.Data()));
      delete canV0CounterAbs;

      canV0CounterRel->cd();
      canV0CounterRel->SetLeftMargin(0.15);
      mgrV0CounterRel->Draw("AP0");
      mgrV0CounterRel->GetXaxis()->Set(iNCategV0, 0, iNCategV0);
      for(Int_t j = 0; j < iNCategV0; j++)
        mgrV0CounterRel->GetXaxis()->SetBinLabel(j + 1, categV0[j].Data());
      mgrV0CounterRel->GetYaxis()->SetTitleOffset(2);
      mgrV0CounterRel->SetMaximum(2);
      mgrV0CounterRel->SetMinimum(1e-6);
      canV0CounterRel->SetLogy();
      legend = canV0CounterRel->BuildLegend(0.2, 0.15, 0.35, 0.3);
      SetLegend(legend);
      canV0CounterRel->SaveAs(Form("canCounterRel%s.%s", V0Name[V0Type].Data(), sImageSuf.Data()));
      delete canV0CounterRel;
    }

    if(iSignalVsCent)
    {
      TCanvas* canInvMassCentAll = new TCanvas("canM", "", iCanWidth, iCanHeight);
      TCanvas* canInvMassCentBg = new TCanvas("canMBg", "", iCanWidth, iCanHeight);
      TCanvas* canInvMassCentSB = new TCanvas("canMSB", "", iCanWidth, iCanHeight);
      TCanvas* canInvMassCentCand = new TCanvas("canMCand", "", iCanWidth, iCanHeight);
      TH1D* hisInvMassCent;
      TH1D* hisInvMassCentSB = new TH1D("hisSB", Form("%s signal purity;centrality", V0Symbol[V0Type].Data()), iNCentBins, 0, iNCentBins);
      TH1D* hisInvMassCentCandNew = new TH1D("hisCandNew", Form("%s signal candidates;centrality;counts", V0Symbol[V0Type].Data()), iNCentBins, 0, iNCentBins);
      TGraphErrors* grInvMassCentAll;
//          TGraphErrors* grInvMassCentBg;
      TMultiGraph* mgrInvMassCentAll = new TMultiGraph();
      TMultiGraph* mgrInvMassCentBg = new TMultiGraph();
      TString hisInvMassCentName = Form("fh1V0InvMass%sCent_%%d", V0Name[V0Type].Data());
      for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
      {
        hisInvMassCent = 0;
        if(!eventsCent[iCent])
          continue;
        hisInvMassCent = GetHistogram1D(listStd, Form(hisInvMassCentName.Data(), iCent));
        if(!hisInvMassCent)
          return;

        Bool_t bStatusBC = kTRUE;
        BinCounterObject* binCount = new BinCounterObject(hisInvMassCent, fMassV0[V0Type], fSigmaV0[V0Type]);
        binCount->SetVerbose(bVerboseBC);
        binCount->SetDegreePolMax(iDegreePolSideBands);
        bStatusBC &= binCount->SetXLimits(dBgOutMin, dBgOutMax);
        if(bFixedBC)
          bStatusBC &= binCount->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax);
        else
          bStatusBC &= binCount->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
        binCount->SetFitOptionSideBands(sOptionFitSB.Data());
//                  binCount->SetBackgroundSubtraction(0);
//              binCount->Plot("canBinCounter0");
        bStatusBC &= binCount->EstimateParameters();
//              binCount->Plot("canBinCounter1");
        bStatusBC &= binCount->Run();
        binCount->Plot("canBinCounter2");
        if(!bStatusBC)
        {
          printf("Something wrong with bin counting\n");
          delete binCount;
          continue;
        }
        Double_t dSignalErr = 0;
        Double_t dSignal = binCount->GetSignalAndError(&dSignalErr);
        printf("Signal: %f +- %f\n", dSignal, dSignalErr);
        printf("Mean: %f +- %f\n", binCount->GetMean(), binCount->GetMeanError());
        printf("Sigma: %f +- %f\n", binCount->GetSigma(), binCount->GetSigmaError());
        Double_t dPurityErr = 0;
        Double_t dPurity = binCount->GetPurityAndError(&dPurityErr);
        hisInvMassCentSB->SetBinContent(iCent + 1, dPurity);
        hisInvMassCentSB->SetBinError(iCent + 1, dPurityErr);
        delete binCount;

        hisInvMassCentCandNew->SetBinContent(iCent + 1, dSignal);
        hisInvMassCentCandNew->SetBinError(iCent + 1, dSignalErr);

        grInvMassCentAll = MakeGraphErrors(hisInvMassCent, GetCentBinLabel(iCent), iMyColors[iCent]);
        mgrInvMassCentAll->Add(grInvMassCentAll);

        hisInvMassCentSB->GetXaxis()->SetBinLabel(iCent + 1, GetCentBinLabel(iCent).Data());
        hisInvMassCentCandNew->GetXaxis()->SetBinLabel(iCent + 1, GetCentBinLabel(iCent).Data());
      }
      canInvMassCentAll->cd();
      canInvMassCentAll->SetLogy();
      mgrInvMassCentAll->SetTitle(Form("%s signal + background;#it{m}_{inv} (GeV/#it{c}^{2});counts", V0Symbol[V0Type].Data()));
      mgrInvMassCentAll->Draw("AP0");
      legend = canInvMassCentAll->BuildLegend(0.7, 0.6, 0.85, 0.85);
      SetLegend(legend);
      canInvMassCentAll->SaveAs(Form("canInvMassCentAll%s.%s", V0Name[V0Type].Data(), sImageSuf.Data()));
      delete canInvMassCentAll;
      delete mgrInvMassCentAll;

//              canInvMassCentBg->cd();
//              canInvMassCentBg->SetLogy();
//              mgrInvMassCentBg->SetTitle(Form("%s background;#it{m}_{inv} (GeV/#it{c}^{2});counts",V0Symbol[V0Type].Data()));
//              mgrInvMassCentBg->Draw("AP0");
//              legend = canInvMassCentBg->BuildLegend(0.7,0.6,0.85,0.85);
//              SetLegend(legend);
//              canInvMassCentBg->SaveAs(Form("canInvMassCentBg%s.%s",V0Name[V0Type].Data(),sImageSuf.Data()));
      delete canInvMassCentBg;
      delete mgrInvMassCentBg;

      canInvMassCentSB->cd();
      hisInvMassCentSB->Draw();
      canInvMassCentSB->SaveAs(Form("canInvMassCentSB%s.%s", V0Name[V0Type].Data(), sImageSuf.Data()));
      delete canInvMassCentSB;

      canInvMassCentCand->cd();
      hisInvMassCentCandNew->Draw();
      canInvMassCentCand->SaveAs(Form("canInvMassCentCand%s.%s", V0Name[V0Type].Data(), sImageSuf.Data()));
      delete canInvMassCentCand;
    }

    if(iSignalVsCut)
    {
      TCanvas* canInvMassCutAll = new TCanvas("can", "", iCanWidth, iCanHeight);
      TCanvas* canInvMassCutBg = new TCanvas("canBg", "", iCanWidth, iCanHeight);
      TCanvas* canInvMassCutSignal = new TCanvas("canSig", "", iCanWidth, iCanHeight);
      TCanvas* canInvMassCutSB = new TCanvas("canSB", "", iCanWidth, iCanHeight);
      TCanvas* canInvMassCutCand = new TCanvas("canCand", "", iCanWidth, iCanHeight);
      TH1D* hisInvMassCut;
      TH1D* hisInvMassCutSB = new TH1D("hisSB", Form("%s signal purity", V0Symbol[V0Type].Data()), iNCategV0, 0, iNCategV0);
      TH1D* hisInvMassCutCand = new TH1D("hisCand", Form("%s signal candidates", V0Symbol[V0Type].Data()), iNCategV0, 0, iNCategV0);
      TGraphErrors* grInvMassCutAll;
      TGraphErrors* grInvMassCutBg;
      TGraphErrors* grInvMassCutSignal;
      TMultiGraph* mgrInvMassCutAll = new TMultiGraph();
      TMultiGraph* mgrInvMassCutBg = new TMultiGraph();
      TMultiGraph* mgrInvMassCutSignal = new TMultiGraph();
      TString hisInvMassCutName = Form("fh1V0InvMass%sAll_%%d", V0Name[V0Type].Data());
      TH1D* hisInvMassCutBg = 0;
      for(Int_t iCut = 0; iCut < iNCategV0; iCut++)
      {
        hisInvMassCut = GetHistogram1D(listStd, Form(hisInvMassCutName.Data(), iCut));
        if(!hisInvMassCut)
          return;
        if(hisInvMassCut->GetEntries() < 100)
          continue;
        grInvMassCutAll = new TGraphErrors(hisInvMassCut);
        grInvMassCutAll->SetLineColor(iMyColors[iCut % iNMyColors]);
        grInvMassCutAll->SetMarkerColor(iMyColors[iCut % iNMyColors]);
        grInvMassCutAll->SetTitle(categV0[iCut].Data());
        mgrInvMassCutAll->Add(grInvMassCutAll);

        // copy of the histogram with excluded peak
        hisInvMassCutBg = new TH1D("hisInvMassCutBg", "", hisInvMassCut->GetNbinsX(), hisInvMassCut->GetBinLowEdge(1), hisInvMassCut->GetBinLowEdge(hisInvMassCut->GetNbinsX() + 1));
        TH1D* hisInvMassCutSignal = new TH1D("hisInvMassCutSignal", "", hisInvMassCut->GetNbinsX(), hisInvMassCut->GetBinLowEdge(1), hisInvMassCut->GetBinLowEdge(hisInvMassCut->GetNbinsX() + 1));
        for(Int_t i = 1; i <= hisInvMassCutBg->GetNbinsX(); i++)
        {
          if((hisInvMassCutBg->GetBinLowEdge(i) > fRangeSignal[V0Type][0]) && (hisInvMassCutBg->GetBinLowEdge(i) < fRangeSignal[V0Type][1]))
            continue;
          hisInvMassCutBg->SetBinContent(i, hisInvMassCut->GetBinContent(i));
          hisInvMassCutBg->SetBinError(i, hisInvMassCut->GetBinError(i));
        }
        // fit background
        if(!(categV0[iCut].Contains("peak")))
        {
          printf("InvMassCutBg: Fitting background: cut=%d, with %.0f entries\n", iCut, hisInvMassCutBg->Integral());
          hisInvMassCutBg->Fit(funBg, "R");
//          funBg = hisInvMassCutBg->GetFunction(funBg);
        }
        grInvMassCutBg = new TGraphErrors(hisInvMassCutBg);
        grInvMassCutBg->SetLineColor(iMyColors[iCut % iNMyColors]);
        grInvMassCutBg->SetMarkerColor(iMyColors[iCut % iNMyColors]);
        grInvMassCutBg->SetTitle(categV0[iCut].Data());
        mgrInvMassCutBg->Add(grInvMassCutBg);
        // get and plot the fit function
        if(funBg)
        {
          TGraph* grFitResultBg = new TGraph(funBg);
          grFitResultBg->SetLineColor(iMyColors[iCut % iNMyColors]);
          grFitResultBg->SetLineStyle(2);
          mgrInvMassCutBg->Add(grFitResultBg);
          mgrInvMassCutAll->Add(grFitResultBg);
        }

        Double_t fIntegralAllBC = 0; // sum of bin content under the peak
        Double_t fIntegralSignalBC = 0; // sum of bin content in the peak only
        Double_t fBinDiff; // signal height in a bin
        Double_t fIntegralAllBCError2 = 0; // square error of the sum under peak
        for(Int_t i = 1; i <= hisInvMassCut->GetNbinsX(); i++)
        {
          if((hisInvMassCut->GetBinLowEdge(i) < fRangeSignal[V0Type][0]) || (hisInvMassCut->GetBinLowEdge(i) > fRangeSignal[V0Type][1]))
            continue;
          if(categV0[iCut].Contains("peak") && ((hisInvMassCut->GetBinLowEdge(i) < fSignalWindow[V0Type][0]) || (hisInvMassCut->GetBinLowEdge(i) > fSignalWindow[V0Type][1])))
            continue;
          fIntegralAllBC += hisInvMassCut->GetBinContent(i);
          fIntegralAllBCError2 += (hisInvMassCut->GetBinError(i)) * (hisInvMassCut->GetBinError(i));
          if(funBg)
            fBinDiff = hisInvMassCut->GetBinContent(i) - funBg->Eval(hisInvMassCut->GetBinCenter(i));
          else
            fBinDiff = hisInvMassCut->GetBinContent(i);
          fIntegralSignalBC += fBinDiff;
          hisInvMassCutSignal->SetBinContent(i, fBinDiff);
          hisInvMassCutSignal->SetBinError(i, hisInvMassCut->GetBinError(i));
        }
        Double_t fSB = fIntegralSignalBC / fIntegralAllBC;
        Double_t fSBError = (fIntegralAllBC - fIntegralSignalBC) / (fIntegralAllBC * fIntegralAllBC) * sqrt(fIntegralAllBCError2);
        printf("Signal purity for %s: %f +- %f\n", categV0[iCut].Data(), fSB, fSBError);
        hisInvMassCutSB->SetBinContent(iCut + 1, fSB);
        hisInvMassCutSB->SetBinError(iCut + 1, fSBError);
        hisInvMassCutSB->GetXaxis()->SetBinLabel(iCut + 1, categV0[iCut].Data());
        hisInvMassCutCand->SetBinContent(iCut + 1, fIntegralSignalBC);
        hisInvMassCutCand->SetBinError(iCut + 1, sqrt(fIntegralAllBCError2));
        hisInvMassCutCand->GetXaxis()->SetBinLabel(iCut + 1, categV0[iCut].Data());
        // plot the extracted signal
        grInvMassCutSignal = new TGraphErrors(hisInvMassCutSignal);
        grInvMassCutSignal->SetTitle(categV0[iCut].Data());
        grInvMassCutSignal->SetLineColor(iMyColors[iCut % iNMyColors]);
        grInvMassCutSignal->SetMarkerColor(iMyColors[iCut % iNMyColors]);
        mgrInvMassCutSignal->Add(grInvMassCutSignal);
      }
      canInvMassCutAll->cd();
      canInvMassCutAll->SetLogy();
      mgrInvMassCutAll->SetMinimum(1);
      mgrInvMassCutAll->SetTitle(Form("%s signal + background;#it{m}_{inv} (GeV/#it{c}^{2});counts", V0Symbol[V0Type].Data()));
      mgrInvMassCutAll->Draw("apl");
      legend = canInvMassCutAll->BuildLegend(0.7, 0.2, 0.85, 0.85);
      SetLegend(legend);
      canInvMassCutAll->SaveAs(Form("canInvMassCutAll%s.%s", V0Name[V0Type].Data(), sImageSuf.Data()));
      delete canInvMassCutAll;

      canInvMassCutBg->cd();
      canInvMassCutBg->SetLogy();
//              mgrInvMassCutBg->SetMinimum(1);
//              mgrInvMassCutBg->SetMinimum(8e4);
//              mgrInvMassCutBg->SetMaximum(1e6);
      mgrInvMassCutBg->SetTitle(Form("%s background;#it{m}_{inv} (GeV/#it{c}^{2});counts", V0Symbol[V0Type].Data()));
      mgrInvMassCutBg->Draw("apl");
      legend = canInvMassCutBg->BuildLegend(0.7, 0.2, 0.85, 0.85);
      SetLegend(legend);
      canInvMassCutBg->SaveAs(Form("canInvMassCutBg%s.%s", V0Name[V0Type].Data(), sImageSuf.Data()));
      delete canInvMassCutBg;

      canInvMassCutSignal->cd();
      canInvMassCutSignal->SetLogy();
      mgrInvMassCutSignal->SetMinimum(1);
      mgrInvMassCutSignal->SetTitle(Form("%s extracted signal;#it{m}_{inv} (GeV/#it{c}^{2});counts", V0Symbol[V0Type].Data()));
      mgrInvMassCutSignal->Draw("apl");
      legend = canInvMassCutSignal->BuildLegend(0.7, 0.2, 0.85, 0.85);
      SetLegend(legend);
      canInvMassCutSignal->SaveAs(Form("canInvMassCutSignal%s.%s", V0Name[V0Type].Data(), sImageSuf.Data()));
      delete canInvMassCutSignal;

      canInvMassCutSB->cd();
      hisInvMassCutSB->Draw();
      canInvMassCutSB->SaveAs(Form("canInvMassCutSB%s.%s", V0Name[V0Type].Data(), sImageSuf.Data()));
      delete canInvMassCutSB;

      canInvMassCutCand->cd();
      canInvMassCutCand->SetLogy();
      hisInvMassCutCand->Draw();
      canInvMassCutCand->SaveAs(Form("canInvMassCutCand%s.%s", V0Name[V0Type].Data(), sImageSuf.Data()));
      delete canInvMassCutCand;
    }

    TCanvas* canPtDensityCompare = new TCanvas("canPtDensityCompare", "", iCanWidth, iCanHeight);
    TMultiGraph* mgrPtDensityCompare = new TMultiGraph();

    if(iPtInclusive)  // inclusive V0 pt spectra
    {
      for(Int_t iBinning = 0; iBinning < 3; iBinning++)
      {
        if(iBinning == 1 && !iPtInclusiveJetBins)
          continue;
        if(iBinning == 2 && !iPtInclusiveSys)
          continue;
        TCanvas* canPtInclusiveAll[iNCentBins];
        TCanvas* canPtInclusiveBg[iNCentBins];
        TCanvas* canPtInclusiveSBAll;
        TCanvas* canPtInclusiveMean;
        TCanvas* canPtInclusiveSigma;
        TCanvas* canPtInclusiveSpectrum = new TCanvas("canPtInclusiveSpectrum", "", iCanWidth, iCanHeight); // pt spectrum, bin counting
        TMultiGraph* mgrPtInclusiveSBAll;
        TMultiGraph* mgrPtInclusiveMean;
        TMultiGraph* mgrPtInclusiveSigma;
        TMultiGraph* mgrPtInclusiveSpectrum = new TMultiGraph();

        TGraphErrors* grPtInclusiveAll;
        TGraphErrors* grPtInclusiveSpectrum;

        canPtInclusiveSBAll = new TCanvas(Form("canSBAll"), "", iCanWidth, iCanHeight);
        canPtInclusiveMean = new TCanvas(Form("canMean"), "", iCanWidth, iCanHeight);
        canPtInclusiveSigma = new TCanvas(Form("canSigma"), "", iCanWidth, iCanHeight);
        mgrPtInclusiveSBAll = new TMultiGraph();
        mgrPtInclusiveMean = new TMultiGraph();
        mgrPtInclusiveSigma = new TMultiGraph();

        Int_t iBinPtInclFirstTmp = iBinPtInclFirst;
        Int_t iBinPtInclLastTmp = iBinPtInclLast;
        Int_t iNBinsPtV0InclTmp = iNBinsPtV0AllLF;
        Double_t* dBinsPtV0InclTmp = dBinsPtV0AllLF;
        TString sAppendix = "";

        if(iBinning == 1)
        {
          iBinPtInclFirstTmp = iBinPtInJetsFirst;
          iBinPtInclLastTmp = iBinPtInJetsLast;
          iNBinsPtV0InclTmp = iNBinsPtV0InJet;
          dBinsPtV0InclTmp = dBinsPtV0InJet;
          sAppendix = "-JetBins";
        }
        if(iBinning == 2)
        {
          iBinPtInclFirstTmp = iBinPtSysFirst;
          iBinPtInclLastTmp = iBinPtSysFirst;
          iNBinsPtV0InclTmp = iNBinsPtV0Sys;
          dBinsPtV0InclTmp = dBinsPtV0Sys;
          sAppendix = "-Sys";
        }

        for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
        {
          printf("Cent %d: Start\n", iCent);
          if(!eventsCent[iCent])
            continue;
          canPtInclusiveAll[iCent] = new TCanvas(Form("canAll%d", iCent), "", iCanWidth, iCanHeight);
          canPtInclusiveBg[iCent] = new TCanvas(Form("canBg%d", iCent), "", iCanWidth, iCanHeight);

          TH1D* hisPtInclusiveSBBC = new TH1D(Form("hisSBBC%d", iCent), Form("%s signal purity, %s;#it{p}_{T} (GeV/c)", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()), iNBinsPtV0InclTmp, dBinsPtV0InclTmp);
          TH1D* hisPtInclusiveSpectrum = new TH1D(Form("hisCand%d", iCent), Form("%s spectrum, %s;#it{p}_{T} (GeV/c)", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()), iNBinsPtV0InclTmp, dBinsPtV0InclTmp);
          TMultiGraph* mgrPtInclusiveAll = new TMultiGraph();
          TMultiGraph* mgrPtInclusiveBg = new TMultiGraph();
          TString spaPtInclusiveName = Form("fhnV0Inclusive%s_C%%d", V0Name[V0Type].Data());

          THnSparseD* spaPtIncl2D = GetSparseD(listStd, Form(spaPtInclusiveName.Data(), iCent));
          if(!spaPtIncl2D)
            return;
          spaPtIncl2D->GetAxis(0)->SetTitle(V0LabelM[V0Type].Data());

          TH1D* hisPtInclusiveMean = new TH1D(Form("hisMean%d", iCent), "Mean", iNBinsPtV0InclTmp, dBinsPtV0InclTmp);
          TH1D* hisPtInclusiveSigma = new TH1D(Form("hisSigma%d", iCent), "Sigma", iNBinsPtV0InclTmp, dBinsPtV0InclTmp);
          TH1D* hisPtInclusiveSigmaExp = new TH1D(Form("hisSigmaExp%d", iCent), "Expected Sigma", iNBinsPtV0InclTmp, dBinsPtV0InclTmp); // expected sigma from formula

          for(Int_t iPt = iBinPtInclFirstTmp; iPt <= iBinPtInclLastTmp; iPt++) // should start from 2
          {
            printf("Pt %d: Start\n", iPt);
            // make projection of mass spectrum in selected pt range
            Int_t iPtBinFirst = spaPtIncl2D->GetAxis(1)->FindBin(dBinsPtV0InclTmp[iPt - 1] + dEpsilon);
            Int_t iPtBinLast = spaPtIncl2D->GetAxis(1)->FindBin(dBinsPtV0InclTmp[iPt] - dEpsilon);
            spaPtIncl2D->GetAxis(1)->SetRange(iPtBinFirst, iPtBinLast);

            TH1D* hisMass = (TH1D*)spaPtIncl2D->Projection(0, "e");
            if(!hisMass)
              return;
            if(iNRebin > 1)
              hisMass->Rebin(iNRebin);
            hisMass->SetName(Form("PtInclusive%s-C%d-%d_B%d", V0Name[V0Type].Data(), iCent, iPt, iBinning));
            Int_t iEntries = hisMass->Integral();
            printf("PtInclusiveBg: Making projection %.1f-%.1f, bins: %d-%d, entries: %d\n", dBinsPtV0InclTmp[iPt - 1], dBinsPtV0InclTmp[iPt], iPtBinFirst, iPtBinLast, iEntries);

            grPtInclusiveAll = MakeGraphErrors(hisMass, Form("%.1f-%.1f", dBinsPtV0InclTmp[iPt - 1], dBinsPtV0InclTmp[iPt]), iMyColors[(iPt - 1) % iNMyColors]);
            mgrPtInclusiveAll->Add(grPtInclusiveAll);

            printf("Pt %d: BC Start\n", iPt);
            Bool_t bStatusBC = kTRUE;
            BinCounterObject* binCount = new BinCounterObject(hisMass, fMassV0[V0Type], fSigmaV0[V0Type]);
            binCount->SetVerbose(bVerboseBC);
            binCount->SetDegreePolMax(iDegreePolSideBands);
            bStatusBC &= binCount->SetXLimits(dBgOutMin, dBgOutMax);
            if(bFixedBC)
            {
              bStatusBC &= binCount->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax); // normal settings
              // settings for counting in the entire range
//              bStatusBC &= binCount->FixRegions(dBgOutMin,dBgOutMax,dBgInMin,dBgInMax,dBgOutMin,dBgOutMax);
            }
            else
              bStatusBC &= binCount->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
            binCount->SetFitOptionSideBands(sOptionFitSB.Data());
//            binCount->SetBackgroundSubtraction(0); // set 0 if you want just counts
            TString sTitle = Form("%s, inclusive, c. %s, pT: %.1f-%.1f GeV/c", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtV0InclTmp[iPt - 1], dBinsPtV0InclTmp[iPt]);
//            binCount->Plot("canBinCounter0",Form("%s, step %d",sTitle.Data(),0));
            bStatusBC &= binCount->EstimateParameters();
//            binCount->Plot("canBinCounter1",Form("%s, step %d",sTitle.Data(),1));
            bStatusBC &= binCount->Run();
            binCount->Plot("canBinCounter2", Form("%s, step %d", sTitle.Data(), 2));
            printf("Pt %d: BC End\n", iPt);
            if(!bStatusBC)
            {
              printf("Something wrong with bin counting\n");
              delete binCount;
              continue;
            }
            Double_t dSignalErr = 0;
            Double_t dSignal = binCount->GetSignalAndError(&dSignalErr);
//            Double_t dBgErr = 0;
//            Double_t dBg = binCount->GetBgAndError(&dBgErr);
            Double_t dPurityErr = 0;
            Double_t dPurity = binCount->GetPurityAndError(&dPurityErr);
            hisPtInclusiveSBBC->SetBinContent(iPt, dPurity);
            hisPtInclusiveSBBC->SetBinError(iPt, dPurityErr);
            printf("Signal: %f +- %f\n", dSignal, dSignalErr);
            printf("Mean: %f +- %f\n", binCount->GetMean(), binCount->GetMeanError());
            printf("Sigma: %f +- %f\n", binCount->GetSigma(), binCount->GetSigmaError());

            // feed-down correction
            if(iCorrFD && (V0Type == 1 || V0Type == 2))
            {
              printf("Correcting %s for feed-down, pT bin %d\n", hisPtInclusiveSpectrum->GetName(), iPt);
              Double_t dFD;
              dFD = GetMeanFeedDownLambdaLF(dBinsPtV0InclTmp[iPt - 1], dBinsPtV0InclTmp[iPt]);
              dSignal *= (1 - dFD);
              dSignalErr *= (1 - dFD);
            }

            hisPtInclusiveSpectrum->SetBinContent(iPt, dSignal);
            hisPtInclusiveSpectrum->SetBinError(iPt, dSignalErr);
            hisPtInclusiveMean->SetBinContent(iPt, binCount->GetMean());
            hisPtInclusiveMean->SetBinError(iPt, binCount->GetMeanError());
            hisPtInclusiveSigma->SetBinContent(iPt, binCount->GetSigma());
            hisPtInclusiveSigma->SetBinError(iPt, binCount->GetSigmaError());
            TF1* funFitResultBg = binCount->GetFunctionBackground();
            TGraph* grFitResultBg = new TGraph(funFitResultBg);
            grFitResultBg->SetTitle(Form("%.1f-%.1f", dBinsPtV0InclTmp[iPt - 1], dBinsPtV0InclTmp[iPt]));
            grFitResultBg->SetLineColor(iMyColors[(iPt - 1) % iNMyColors]);
            grFitResultBg->SetLineStyle(2);
            grFitResultBg->SetFillStyle(0);
            mgrPtInclusiveBg->Add(grFitResultBg, "l");
            delete binCount;

            Double_t fSignalSigma = MassPeakSigmaOld((dBinsPtV0InclTmp[iPt - 1] + dBinsPtV0InclTmp[iPt]) / 2, V0Type);
            hisPtInclusiveSigmaExp->SetBinContent(iPt, fSignalSigma);
            hisPtInclusiveSigmaExp->SetBinError(iPt, 0);
            printf("Pt %d: End\n", iPt);
          }

          canPtInclusiveAll[iCent]->cd();
          canPtInclusiveAll[iCent]->SetLogy();
          mgrPtInclusiveAll->SetTitle(Form("%s signal + background, %s;m;counts", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
          mgrPtInclusiveAll->SetMinimum(1);
          mgrPtInclusiveAll->Draw("AP0");
          legend = canPtInclusiveAll[iCent]->BuildLegend(0.7, 0.2, 0.85, 0.85);
          SetLegend(legend);
          canPtInclusiveAll[iCent]->SaveAs(Form("canPtInclusiveAll%s-C%d%s.%s", V0Name[V0Type].Data(), iCent, sAppendix.Data(), sImageSuf.Data()));
          delete canPtInclusiveAll[iCent];

          canPtInclusiveBg[iCent]->cd();
          canPtInclusiveBg[iCent]->SetLogy();
          mgrPtInclusiveBg->SetTitle(Form("%s background, %s", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
          mgrPtInclusiveBg->Draw("AP0");
          legend = canPtInclusiveBg[iCent]->BuildLegend(0.7, 0.2, 0.85, 0.85);
          SetLegend(legend);
          canPtInclusiveBg[iCent]->SaveAs(Form("canPtInclusiveBg%s-C%d%s.%s", V0Name[V0Type].Data(), iCent, sAppendix.Data(), sImageSuf.Data()));
          delete canPtInclusiveBg[iCent];

          TGraphErrors* grSBBCAll = MakeGraphErrors(hisPtInclusiveSBBC, Form("%s BC", GetCentBinLabel(iCent).Data()), iMyColors[iCent], iMyMarkersEmpty[iCent]);
          mgrPtInclusiveSBAll->Add(grSBBCAll);
          delete hisPtInclusiveSBBC;
          TGraphErrors* grMean = MakeGraphErrors(hisPtInclusiveMean, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersFull[iCent]);
          mgrPtInclusiveMean->Add(grMean);
          delete hisPtInclusiveMean;
          TGraphErrors* grSigma = MakeGraphErrors(hisPtInclusiveSigma, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersFull[iCent]);
          mgrPtInclusiveSigma->Add(grSigma);
          delete hisPtInclusiveSigma;
//          TGraphErrors* grSigmaExp = MakeGraphErrors(hisPtInclusiveSigmaExp,"expected",iMyColors[iCent],iMyMarkersEmpty[iCent]);
//          mgrPtInclusiveSigma->Add(grSigmaExp);
          delete hisPtInclusiveSigmaExp;

          dirOutSpectra->cd();
          // Event normalization
          if(!bYieldsOnly)
            hisPtInclusiveSpectrum->Scale(1. / eventsCent[iCent], "width");
          if(bInclusiveDensity)
            hisPtInclusiveSpectrum->Scale(1. / fNormalisationV0Area);
          if(!iCorrFD)
            hisPtInclusiveSpectrum->Write(Form("fh1PtInclusive%s_C%d_Raw%s", V0Name[V0Type].Data(), iCent, sAppendix.Data()));
          // Efficiency correction
          if(iCorrEff)
          {
            printf("Correcting %s for efficiency\n", hisPtInclusiveSpectrum->GetName());
            if(iBinning == 1)
              hisPtInclusiveSpectrum = (TH1D*)DivideHistograms1D(hisPtInclusiveSpectrum, hisEffPtInJetsExt[iCent]);
            else
              hisPtInclusiveSpectrum = (TH1D*)DivideHistograms1D(hisPtInclusiveSpectrum, hisEffPtInclExt[iCent]);
            if(!hisPtInclusiveSpectrum)
              return;
          }
          hisPtInclusiveSpectrum->SetTitle(Form("%s: #it{p}_{T} spectrum, inclusive, c. %s;#it{p}_{T}^{%s} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{%s}}{d#it{p}_{T}^{%s}} (#it{c} GeV^{#minus1})", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data()));
          hisPtInclusiveSpectrum->Write(Form("fh1PtInclusive%s_C%d%s", V0Name[V0Type].Data(), iCent, sAppendix.Data()));
          grPtInclusiveSpectrum = MakeGraphErrors(hisPtInclusiveSpectrum, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersFull[iCent]);
          mgrPtInclusiveSpectrum->Add(grPtInclusiveSpectrum);
          delete hisPtInclusiveSpectrum;
          printf("Cent %d: End\n", iCent);
        }

        canPtInclusiveMean->cd();
        canPtInclusiveMean->SetLeftMargin(0.15);
        mgrPtInclusiveMean->SetTitle(Form("%s: Mean #it{m}_{inv};#it{p}_{T} (GeV/c);mean (GeV/c^{2})", V0Symbol[V0Type].Data()));
        if(V0Type == 0)
        {
          mgrPtInclusiveMean->SetMinimum(0.4976);
          mgrPtInclusiveMean->SetMaximum(0.4997);
        }
        else
        {
          mgrPtInclusiveMean->SetMinimum(1.11557);
          mgrPtInclusiveMean->SetMaximum(1.1165);
        }
        mgrPtInclusiveMean->Draw("AP0");
        mgrPtInclusiveMean->GetYaxis()->SetTitleOffset(2);
        mgrPtInclusiveMean->GetXaxis()->SetLimits(fPtAllXMin, fPtAllXMax);
        legend = canPtInclusiveMean->BuildLegend(0.2, 0.6, 0.35, 0.85);
        SetLegend(legend);
        canPtInclusiveMean->SaveAs(Form("canPtInclusiveMean%s%s.%s", V0Name[V0Type].Data(), sAppendix.Data(), sImageSuf.Data()));
        delete canPtInclusiveMean;
        delete mgrPtInclusiveMean;

        canPtInclusiveSigma->cd();
        canPtInclusiveSigma->SetLeftMargin(0.15);
        mgrPtInclusiveSigma->SetTitle(Form("%s: Sigma #it{m}_{inv};#it{p}_{T} (GeV/c);#sigma (GeV/c^{2})", V0Symbol[V0Type].Data()));
        if(V0Type == 0)
        {
          mgrPtInclusiveSigma->SetMinimum(0.0035);
          mgrPtInclusiveSigma->SetMaximum(0.011);
        }
        else
        {
          mgrPtInclusiveSigma->SetMinimum(0.0013);
          mgrPtInclusiveSigma->SetMaximum(0.0033);
        }
        mgrPtInclusiveSigma->Draw("AP0");
        mgrPtInclusiveSigma->GetYaxis()->SetTitleOffset(2);
        mgrPtInclusiveSigma->GetXaxis()->SetLimits(fPtAllXMin, fPtAllXMax);
        legend = canPtInclusiveSigma->BuildLegend(0.7, 0.15, 0.85, 0.45);
        SetLegend(legend);
        canPtInclusiveSigma->SaveAs(Form("canPtInclusiveSigma%s%s.%s", V0Name[V0Type].Data(), sAppendix.Data(), sImageSuf.Data()));
        delete canPtInclusiveSigma;
        delete mgrPtInclusiveSigma;

        canPtInclusiveSBAll->cd();
        canPtInclusiveSBAll->SetLeftMargin(0.15);
        mgrPtInclusiveSBAll->SetTitle(Form("%s signal purity;#it{p}_{T} (GeV/c);purity", V0Symbol[V0Type].Data()));
        mgrPtInclusiveSBAll->SetMinimum(0.);
        mgrPtInclusiveSBAll->SetMaximum(1.1);
        mgrPtInclusiveSBAll->Draw("AP0");
        mgrPtInclusiveSBAll->GetYaxis()->SetTitleOffset(2);
        mgrPtInclusiveSBAll->GetXaxis()->SetLimits(fPtAllXMin, fPtAllXMax);
        legend = canPtInclusiveSBAll->BuildLegend(0.7, 0.2, 0.85, 0.45);
        SetLegend(legend);
        canPtInclusiveSBAll->SaveAs(Form("canPtInclusivePurity%s%s.%s", V0Name[V0Type].Data(), sAppendix.Data(), sImageSuf.Data()));
        delete canPtInclusiveSBAll;
        delete mgrPtInclusiveSBAll;
        canPtInclusiveSpectrum->cd();
        canPtInclusiveSpectrum->SetLeftMargin(0.16);
        mgrPtInclusiveSpectrum->SetTitle(Form("%s: #it{p}_{T} spectrum, inclusive;#it{p}_{T}^{h} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{h}}{d#it{p}_{T}^{h}} (#it{c} GeV^{#minus1})", V0Symbol[V0Type].Data()));
        if(!bYieldsOnly)
        {
          mgrPtInclusiveSpectrum->SetMinimum(fPtAllYMin);
          mgrPtInclusiveSpectrum->SetMaximum(fPtAllYMax);
        }
        mgrPtInclusiveSpectrum->Draw("AP0");
        mgrPtInclusiveSpectrum->GetYaxis()->SetTitleOffset(2);
        mgrPtInclusiveSpectrum->GetXaxis()->SetLimits(fPtAllXMin, fPtAllXMax);
        if(!bYieldsOnly)
          canPtInclusiveSpectrum->SetLogy();
        legend = canPtInclusiveSpectrum->BuildLegend(0.7, 0.6, 0.85, 0.85);
        SetLegend(legend);
        labelSystem = labelCollision->DrawLatex(fPtAllXMin + (fPtAllXMax - fPtAllXMin) / 5., 10 * fPtAllYMin, sLabelCollisionText.Data());
        canPtInclusiveSpectrum->SaveAs(Form("canPtInclusiveSpectrum%s%s.%s", V0Name[V0Type].Data(), sAppendix.Data(), sImageSuf.Data()));
        delete canPtInclusiveSpectrum;
        delete mgrPtInclusiveSpectrum;
      }
    }

    if(iPtEtaInclusive)  // inclusive V0 pt spectra
    {
      Int_t iNBinsPtV0Tmp, iBinPtFirstTmp, iBinPtLastTmp;
      Double_t* dBinsPtV0Tmp;
      TH2D* hisEffPtEtaInclBaseTmp = 0;
      TH1D* hisEffPtScaled = 0;
      if(iSwitchPtV0 == 0) // inclusive pt bins
      {
        iNBinsPtV0Tmp = iNBinsPtV0AllLF;
        dBinsPtV0Tmp = dBinsPtV0AllLF;
        iBinPtFirstTmp = iBinPtInclFirst;
        iBinPtLastTmp = iBinPtInclLast;
      }
      if(iSwitchPtV0 == 1) // in-jet pt bins
      {
        iNBinsPtV0Tmp = iNBinsPtV0InJet;
        dBinsPtV0Tmp = dBinsPtV0InJet;
        iBinPtFirstTmp = iBinPtInJetsFirst;
        iBinPtLastTmp = iBinPtInJetsLast;
      }

      for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
      {
        printf("Cent %d: Start\n", iCent);
        if(!eventsCent[iCent])
          continue;
        TH2D* hisPtEtaInclusivePurity = new TH2D(Form("hisPtEtaInclPurity_C%d", iCent), Form("%s signal purity, %s;#it{p}_{T} (GeV/c);#eta", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()), iNBinsPtV0Tmp, dBinsPtV0Tmp, iNBinsEtaInJet, dBinsEtaInJet);
        TH2D* hisPtEtaInclusiveSigma = new TH2D(Form("hisPtEtaInclSigma_C%d", iCent), Form("%s signal sigma, %s;#it{p}_{T} (GeV/c);#eta", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()), iNBinsPtV0Tmp, dBinsPtV0Tmp, iNBinsEtaInJet, dBinsEtaInJet);
        TH2D* hisPtEtaInclusiveSpectrum = new TH2D(Form("hisPtEtaInclSignal_C%d", iCent), Form("%s signal, %s;#it{p}_{T} (GeV/c);#eta", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()), iNBinsPtV0Tmp, dBinsPtV0Tmp, iNBinsEtaInJet, dBinsEtaInJet);
        TH2D* hisPtEtaInclusiveSpectrumRaw = new TH2D(Form("hisPtEtaInclCandidates_C%d", iCent), Form("%s candidates, %s;#it{p}_{T} (GeV/c);#eta", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()), iNBinsPtV0Tmp, dBinsPtV0Tmp, iNBinsEtaInJet, dBinsEtaInJet);
//              hisPtEtaInclusiveSpectrumRaw->Sumw2();
        TString spaPtInclusiveName = Form("fhnV0Inclusive%s_C%%d", V0Name[V0Type].Data());

        if(iCorrEff && iCorrEffInclEta)
        {
          if(iSwitchPtV0 == 0)
            hisEffPtEtaInclBaseTmp = hisEffPtEtaInclBaseInclBins[iCent];
          if(iSwitchPtV0 == 1)
            hisEffPtEtaInclBaseTmp = hisEffPtEtaInclBaseJetBins[iCent];
        }

        THnSparseD* spaPtEtaIncl = GetSparseD(listStd, Form(spaPtInclusiveName.Data(), iCent));
        if(!spaPtEtaIncl)
          return;
        if(bFixedBC)
        {
          Int_t iBinMFirst = spaPtEtaIncl->GetAxis(0)->FindBin(dSigMin + dEpsilon);
          Int_t iBinMLast = spaPtEtaIncl->GetAxis(0)->FindBin(dSigMax - dEpsilon);
          spaPtEtaIncl->GetAxis(0)->SetRange(iBinMFirst, iBinMLast);
        }
        spaPtEtaIncl->GetAxis(0)->SetTitle(V0LabelM[V0Type].Data());
        spaPtEtaIncl->GetAxis(1)->SetRange(0, 0);
        TH2D* hisPtEtaInclProj = (TH2D*)spaPtEtaIncl->Projection(2, 1, "e");
        if(!hisPtEtaInclProj)
          return;
        if(!RebinHistogram2D(hisPtEtaInclProj, hisPtEtaInclusiveSpectrumRaw))
          return;
        dirOutSpectra->cd();
        hisPtEtaInclusiveSpectrumRaw->Write(Form("fh2PtEtaIncl%s_C%d-Candidates", V0Name[V0Type].Data(), iCent));
        spaPtEtaIncl->GetAxis(0)->SetRange(0, 0);

        for(Int_t iPt = iBinPtFirstTmp; iPt <= iBinPtLastTmp; iPt++) // should start from 2
//              for (Int_t iPt=iBinPtInclFirst; iPt<=iBinPtInclFirst+5; iPt++) // should start from 2
        {
          printf("Pt %d: Start\n", iPt);
          // make projection of mass spectrum in selected pt range
          Int_t iPtBinFirst = spaPtEtaIncl->GetAxis(1)->FindBin(dBinsPtV0Tmp[iPt - 1] + dEpsilon);
          Int_t iPtBinLast = spaPtEtaIncl->GetAxis(1)->FindBin(dBinsPtV0Tmp[iPt] - dEpsilon);
          spaPtEtaIncl->GetAxis(1)->SetRange(iPtBinFirst, iPtBinLast);
          for(Int_t iEta = 0; iEta < iNBinsEtaInJet; iEta++)
          {
            Int_t iEtaPtFirstBin = spaPtEtaIncl->GetAxis(2)->FindBin(dBinsEtaInJet[iEta] + dEpsilon);
            Int_t iEtaPtLastBin = spaPtEtaIncl->GetAxis(2)->FindBin(dBinsEtaInJet[iEta + 1] - dEpsilon);
            spaPtEtaIncl->GetAxis(2)->SetRange(iEtaPtFirstBin, iEtaPtLastBin);

            TH1D* hisMass = (TH1D*)spaPtEtaIncl->Projection(0, "e");
            if(!hisMass)
              return;
            if(iNRebin > 1)
              hisMass->Rebin(iNRebin);
            hisMass->SetName(Form("PtEtaInclusive%s-C%d-Pt%d-E%d", V0Name[V0Type].Data(), iCent, iPt, iEta));
            Int_t iEntries = hisMass->Integral();
            printf("PtInclusiveBg: Making projection %.1f-%.1f, bins: %d-%d, entries: %d\n", dBinsPtV0Tmp[iPt - 1], dBinsPtV0Tmp[iPt], iPtBinFirst, iPtBinLast, iEntries);

            printf("Pt %d, Eta: %d: BC Start\n", iPt, iEta);
            Bool_t bStatusBC = kTRUE;
            BinCounterObject* binCount = new BinCounterObject(hisMass, fMassV0[V0Type], fSigmaV0[V0Type]);
            binCount->SetVerbose(bVerboseBC);
            binCount->SetDegreePolMax(iDegreePolSideBands);
            bStatusBC &= binCount->SetXLimits(dBgOutMin, dBgOutMax);
            if(bFixedBC)
            {
              bStatusBC &= binCount->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax); // normal settings
              // settings for counting in the entire range
//                          bStatusBC &= binCount->FixRegions(dBgOutMin,dBgOutMax,dBgInMin,dBgInMax,dBgOutMin,dBgOutMax);
            }
            else
              bStatusBC &= binCount->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
            binCount->SetFitOptionSideBands(sOptionFitSB.Data());
//                      binCount->SetBackgroundSubtraction(0); // set 0 if you want just counts
            TString sTitle = Form("%s, inclusive, c. %s, pT: %g-%g GeV/c, #eta: %g - %g", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtV0Tmp[iPt - 1], dBinsPtV0Tmp[iPt], dBinsEtaInJet[iEta], dBinsEtaInJet[iEta + 1]);
//                      binCount->Plot("canBinCounter0",Form("%s, step %d",sTitle.Data(),0));
            bStatusBC &= binCount->EstimateParameters();
//                      binCount->Plot("canBinCounter1",Form("%s, step %d",sTitle.Data(),1));
            bStatusBC &= binCount->Run();
            binCount->Plot("canBinCounter2", Form("%s, step %d", sTitle.Data(), 2));
            printf("Pt %d, Eta %d: BC End\n", iPt, iEta);
            if(!bStatusBC)
            {
              printf("Something wrong with bin counting\n");
              delete binCount;
              continue;
            }
            Double_t dSignalErr = 0;
            Double_t dSignal = binCount->GetSignalAndError(&dSignalErr);
            Double_t dPurityErr = 0;
            Double_t dPurity = binCount->GetPurityAndError(&dPurityErr);
            hisPtEtaInclusiveSpectrum->SetBinContent(iPt, iEta + 1, dSignal);
            hisPtEtaInclusiveSpectrum->SetBinError(iPt, iEta + 1, dSignalErr);
            hisPtEtaInclusivePurity->SetBinContent(iPt, iEta + 1, dPurity);
            hisPtEtaInclusivePurity->SetBinError(iPt, iEta + 1, dPurityErr);
            hisPtEtaInclusiveSigma->SetBinContent(iPt, iEta + 1, binCount->GetSigma());
            hisPtEtaInclusiveSigma->SetBinError(iPt, iEta + 1, binCount->GetSigmaError());
//                      printf("Signal: %f +- %f\n",dSignal,dSignalErr);
//                      printf("Mean: %f +- %f\n",binCount->GetMean(),binCount->GetMeanError());
//                      printf("Sigma: %f +- %f\n",binCount->GetSigma(),binCount->GetSigmaError());
            delete binCount;
            printf("Pt %d: End\n", iPt);
          }
        }
        hisPtEtaInclusiveSpectrum->Write(Form("fh2PtEtaIncl%s_C%d-Signal", V0Name[V0Type].Data(), iCent));
        hisPtEtaInclusivePurity->Write(Form("fh2PtEtaIncl%s_C%d-Purity", V0Name[V0Type].Data(), iCent));
        hisPtEtaInclusiveSigma->Write(Form("fh2PtEtaIncl%s_C%d-Sigma", V0Name[V0Type].Data(), iCent));
        TH2D* hisPtEtaInclSigVsRaw = DivideHistograms2D(hisPtEtaInclusiveSpectrum, hisPtEtaInclusiveSpectrumRaw);
        if(!hisPtEtaInclSigVsRaw)
          return;
        hisPtEtaInclSigVsRaw->Write(Form("fh2PtEtaIncl%s_C%d-RatioSigVsCand", V0Name[V0Type].Data(), iCent));
        TH2D* hisPtEtaInclPurVsRatio = DivideHistograms2D(hisPtEtaInclusivePurity, hisPtEtaInclSigVsRaw);
        if(!hisPtEtaInclPurVsRatio)
          return;
        hisPtEtaInclPurVsRatio->Write(Form("fh2PtEtaIncl%s_C%d-PurVsRatio", V0Name[V0Type].Data(), iCent));
        TH1D* hisPtInclusiveSpectrumProj = hisPtEtaInclusiveSpectrum->ProjectionX(Form("%s_px", hisPtEtaInclusiveSpectrum->GetName()), 0, -1, "e");
        hisPtInclusiveSpectrumProj->Scale(1. / eventsCent[iCent], "width");
        hisPtInclusiveSpectrumProj->Write(Form("fh1PtInclusiveProj%s_C%d", V0Name[V0Type].Data(), iCent));
        hisPtEtaInclusiveSpectrum->Scale(1. / eventsCent[iCent]);
        hisPtEtaInclusiveSpectrum->Write(Form("fh2PtEtaIncl%s_C%d-Spectrum", V0Name[V0Type].Data(), iCent));
        if(iPtInclusive && iSwitchPtV0 == 0)
        {
          TH1D* hisPtInclusiveSpectrumTrue = (TH1D*)dirOutSpectra->Get(Form("fh1PtInclusive%s_C%d_Raw", V0Name[V0Type].Data(), iCent));
          if(!hisPtInclusiveSpectrumTrue)
          {
            printf("Loading raw spectrum failed. (Error)\n");
            return;
          }
          TH1D* hisRatioPtInclusive = DivideHistograms1D(hisPtInclusiveSpectrumProj, hisPtInclusiveSpectrumTrue);
          if(!hisRatioPtInclusive)
            return;
          hisRatioPtInclusive->Write(Form("fh1PtInclusiveRatio%s_C%d", V0Name[V0Type].Data(), iCent));

          if(iCorrEff && iCorrEffInclEta)
          {
            hisEffPtScaled = GetScaledEfficiency(hisPtEtaInclusiveSpectrum, hisEffPtEtaInclBaseTmp, 0);
            hisEffPtScaled->Write(Form(sNameHisEffPtIncl.Data(), V0Name[V0Type].Data(), iCent));
            TH1D* hisPtInclusiveSpectrum = DivideHistograms1D(hisPtInclusiveSpectrumTrue, hisEffPtScaled);
            if(!hisPtInclusiveSpectrum)
              return;
            hisPtInclusiveSpectrum->Write(Form("fh1PtInclusive%s_C%d_Eta", V0Name[V0Type].Data(), iCent));
            TH2D* hisPtEtaInclusiveSpectrumCorrEff = DivideHistograms2D(hisPtEtaInclusiveSpectrum, hisEffPtEtaInclBaseTmp);
            if(!hisPtEtaInclusiveSpectrumCorrEff)
              return;
            hisPtEtaInclusiveSpectrumCorrEff->Write(Form("fh2PtEtaIncl%s_C%d-SpectrumCorrEff", V0Name[V0Type].Data(), iCent));
          }
        }
        printf("Cent %d: End\n", iCent);
      }
    }

    if(iPtInJets)
    {
      printf("In jets\n");
      TString spaPtInJetsName = Form("fhnV0InJet%s_%%d", V0Name[V0Type].Data());
      TString spaPtInPCName = Form("fhnV0InPerp%s_%%d", V0Name[V0Type].Data());
      TString spaPtInRCName = Form("fhnV0InRnd%s_%%d", V0Name[V0Type].Data());
      TString spaPtInMCCName = Form("fhnV0InMed%s_%%d", V0Name[V0Type].Data());
      TString spaPtOutJetsName = Form("fhnV0OutJet%s_%%d", V0Name[V0Type].Data());
      TString spaPtNoJetsName = Form("fhnV0NoJet%s_%%d", V0Name[V0Type].Data());
      TString hisPtInRCNormName = "fh1NRndConeCent";
      TString hisPtInMCCNormName = "fh1NMedConeCent";
      TString hisPtInOCNormName = "fh1AreaExcluded";
      TString hisPtInNJNormName = "fh1EventCounterCutCent_%d";

      TString sLabelsAxesSpectra = "";
      if(bNormPerArea)
        sLabelsAxesSpectra = "#it{p}_{T}^{h} (GeV/#it{c});#frac{1}{#Delta#it{#eta} #Delta#it{#phi}} #frac{d#it{N}}{d#it{p}_{T}} (#it{c} GeV^{#minus1})";
      else
        sLabelsAxesSpectra = "#it{p}_{T}^{h} (GeV/#it{c});#frac{1}{#it{N}_{jet}} #frac{d#it{N}}{d#it{p}_{T}} (#it{c} GeV^{#minus1})";

      TH1D* hisPtInRCNorm = GetHistogram1D(listStd, hisPtInRCNormName.Data());
      if(!hisPtInRCNorm)
        return;
      TH1D* hisPtInMCCNorm = GetHistogram1D(listStd, hisPtInMCCNormName.Data());
      if(!hisPtInMCCNorm)
        return;

      printf("Statistics: Jets (pt bins):");
      for(Int_t iJet = iJetMin; iJet <= iJetMax; iJet++)
        printf(" %.1f", dBinsPtJet[iJet]);
      printf(", open bins: %d\n", bOpenPtBins);

      // bulk methods
      // random cones
      TCanvas* canPtInRC = new TCanvas("canPtInRC", "", iCanWidth, iCanHeight);
      TMultiGraph* mgrPtInRC = new TMultiGraph();
      // outside cones
      TCanvas* canPtInOC = new TCanvas("canPtInOC", "", iCanWidth, iCanHeight);
      TMultiGraph* mgrPtInOC = new TMultiGraph();
      // no jet events
      TCanvas* canPtInNJ = new TCanvas("canPtInNJ", "", iCanWidth, iCanHeight);
      TMultiGraph* mgrPtInNJ = new TMultiGraph();

      Double_t dArrayJetNorm[iNBinsPtJet];
      Double_t dArrayJetNormErr[iNBinsPtJet];
      Int_t iNBinsPtJetNorm = hisJetPt[0]->GetXaxis()->GetNbins();

      for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
      {
        printf("Cent %d: Start\n", iCent);
        if(!eventsCent[iCent])
          continue;

        // jet spectrum for normalisation
//                  TH1D* hisJetPtNorm = (TH1D*)hisJetPt[iCent]->Clone(Form("hisJetPtNorm%d",iCent));
//                  hisJetPtNorm = (TH1D*)hisJetPtNorm->Rebin(iNBinsPtJet,Form("hisJetPtNorm%d",iCent),dBinsPtJet);
        for(Int_t i = 0; i < iNBinsPtJet; i++)
        {
          Int_t iBinFirst = hisJetPt[iCent]->GetXaxis()->FindBin(dBinsPtJet[i] + dEpsilon);
          Int_t iBinLast;
          if(bOpenPtBins)
            iBinLast = iNBinsPtJetNorm + 1;
          else
            iBinLast = hisJetPt[iCent]->GetXaxis()->FindBin(dBinsPtJet[i + 1] - dEpsilon);
          dArrayJetNorm[i] = hisJetPt[iCent]->IntegralAndError(iBinFirst, iBinLast, dArrayJetNormErr[i]);
          dArrayJetNormErr[i] = 0;
        }

        printf("Statistics: Jets (Cent %d):", iCent);
        for(Int_t iJet = iJetMin; iJet <= iJetMax; iJet++)
          printf(" %d", int(dArrayJetNorm[iJet]));
        printf("\n");

        // V0s in jets (jet cones minus bulk)
        TCanvas* canPtInJets = new TCanvas(Form("canPtInJets-%d", iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrPtInJets = new TMultiGraph();
        // V0s in jet cones
        TCanvas* canPtInJC = new TCanvas(Form("canPtInJC-%d", iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrPtInJC = new TMultiGraph();
        // bulk methods dependent on jet pT
        // perpendicular cones
        TCanvas* canPtInPC = new TCanvas(Form("canPtInPC-%d", iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrPtInPC = new TMultiGraph();
        // compare bulk methods in each centrality bin
        TCanvas* canPtInBulk = new TCanvas("canPtInBulk", "", iCanWidth, iCanHeight);
        TMultiGraph* mgrPtInBulk = new TMultiGraph();
        TCanvas* canPtInBulkCompare = new TCanvas("canPtInBulkCompare", "", iCanWidth, iCanHeight);
        TMultiGraph* mgrPtInBulkCompare = new TMultiGraph();

        TCanvas* canEffCompare = new TCanvas("canEffCompare", "", iCanWidth, iCanHeight);
        TMultiGraph* mgrEffCompare = new TMultiGraph();

        TCanvas* canPurityInJC = new TCanvas("canPurityInJC", "", iCanWidth, iCanHeight);
        TMultiGraph* mgrPurityInJC = new TMultiGraph();

        // V0 in jet cones
        THnSparseD* spaPtInJets = GetSparseD(listStd, Form(spaPtInJetsName.Data(), iCent));
        if(!spaPtInJets)
          return;
        // V0 in perp. cones
        THnSparseD* spaPtInPC = GetSparseD(listStd, Form(spaPtInPCName.Data(), iCent));
        if(!spaPtInPC)
          return;
//              printf("PC %d tot integral: %g\n",iCent,spaPtInPC->ComputeIntegral());
        // V0 in rnd. cones
        THnSparseD* spaPtInRC = GetSparseD(listStd, Form(spaPtInRCName.Data(), iCent));
        if(!spaPtInRC)
          return;
        // V0 in med. cones
        THnSparseD* spaPtInMCC = GetSparseD(listStd, Form(spaPtInMCCName.Data(), iCent));
        if(!spaPtInMCC)
          return;
        // V0 in outside jet cones
        THnSparseD* spaPtOutJet = GetSparseD(listStd, Form(spaPtOutJetsName.Data(), iCent));
        if(!spaPtOutJet)
          return;
        // V0 in in no jet events
        THnSparseD* spaPtNoJet = GetSparseD(listStd, Form(spaPtNoJetsName.Data(), iCent));
        if(!spaPtNoJet)
          return;

        TCanvas* canPtInRCAll = new TCanvas(Form("canInRCAll%s-%d", V0Name[V0Type].Data(), iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrPtInRCAll = new TMultiGraph();
        TCanvas* canPtInRCBg = new TCanvas(Form("canInRCBg%s-%d", V0Name[V0Type].Data(), iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrPtInRCBg = new TMultiGraph();

        TCanvas* canPtInNJAll = new TCanvas(Form("canInNJAll%s-%d", V0Name[V0Type].Data(), iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrPtInNJAll = new TMultiGraph();
        TCanvas* canPtInNJBg = new TCanvas(Form("canInNJBg%s-%d", V0Name[V0Type].Data(), iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrPtInNJBg = new TMultiGraph();

        TCanvas* canPtInOCAll = new TCanvas(Form("canInOCAll%s-%d", V0Name[V0Type].Data(), iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrPtInOCAll = new TMultiGraph();
        TCanvas* canPtInOCBg = new TCanvas(Form("canInOCBg%s-%d", V0Name[V0Type].Data(), iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrPtInOCBg = new TMultiGraph();

        // pt spectrum in rnd. cones
        TH1D* hisPtInRCSpectrum = new TH1D(Form("hisPtInRCSpectrum-%d", iCent), Form("%s: #it{p}_{T} spectrum, in Rnd. cones, c. %s;#it{p}_{T}^{%s} (GeV/#it{c});#frac{1}{#it{N}_{jet}} #frac{d#it{N}^{%s}}{d#it{p}_{T}^{%s}} (#it{c} GeV^{#minus1})", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data()), iNBinsPtV0InJet, dBinsPtV0InJet);
        hisPtInRCSpectrum->Sumw2();
        TH2D* hisPtEtaInRCSpectrum = new TH2D(Form("hisPtEtaInRCSpectrum-%d", iCent), Form("%s: #it{p}_{T}-#eta spectrum, in Rnd. cones, c. %s;#it{p}_{T}^{%s} (GeV/#it{c});#eta", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), V0Symbol[V0Type].Data()), iNBinsPtV0InJet, dBinsPtV0InJet, iNBinsEtaInJet, dBinsEtaInJet);
        TH1D* hisEffInRC;
        // pt spectrum in med. cones
        TH1D* hisPtInMCCSpectrum = new TH1D(Form("hisPtInMCCSpectrum-%d", iCent), Form("%s: #it{p}_{T} spectrum, in M-C cones, c. %s;#it{p}_{T}^{%s} (GeV/#it{c});#frac{1}{#it{N}_{jet}} #frac{d#it{N}^{%s}}{d#it{p}_{T}^{%s}} (#it{c} GeV^{#minus1})", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data()), iNBinsPtV0InJet, dBinsPtV0InJet);
        hisPtInMCCSpectrum->Sumw2();
        TH2D* hisPtEtaInMCCSpectrum = new TH2D(Form("hisPtEtaInMCCSpectrum-%d", iCent), Form("%s: #it{p}_{T} spectrum, in M-C cones, c. %s;#it{p}_{T}^{%s} (GeV/#it{c});#eta", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), V0Symbol[V0Type].Data()), iNBinsPtV0InJet, dBinsPtV0InJet, iNBinsEtaInJet, dBinsEtaInJet);
        TH1D* hisEffInMCC;
        // pt spectrum in no jet events
        TH1D* hisPtInNJSpectrum = new TH1D(Form("hisPtInNJSpectrum-%d", iCent), Form("%s: #it{p}_{T} spectrum, in no jet events, c. %s;#it{p}_{T}^{%s} (GeV/#it{c});#frac{1}{#it{N}_{jet}} #frac{d#it{N}^{%s}}{d#it{p}_{T}^{%s}} (#it{c} GeV^{#minus1})", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data()), iNBinsPtV0InJet, dBinsPtV0InJet);
        hisPtInNJSpectrum->Sumw2();
        TH2D* hisPtEtaInNJSpectrum = new TH2D(Form("hisPtEtaInNJSpectrum-%d", iCent), Form("%s: #it{p}_{T} spectrum, in no jet events, c. %s;#it{p}_{T}^{%s} (GeV/#it{c});#eta", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), V0Symbol[V0Type].Data()), iNBinsPtV0InJet, dBinsPtV0InJet, iNBinsEtaInJet, dBinsEtaInJet);
        TH1D* hisEffInNJ;
        // pt spectrum outside jet cones
        TH1D* hisPtInOCSpectrum = new TH1D(Form("hisPtInOCSpectrum-%d", iCent), Form("%s: #it{p}_{T} spectrum, outside jet cones, c. %s;#it{p}_{T}^{%s} (GeV/#it{c});#frac{1}{#it{N}_{jet}} #frac{d#it{N}^{%s}}{d#it{p}_{T}^{%s}} (#it{c} GeV^{#minus1})", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data()), iNBinsPtV0InJet, dBinsPtV0InJet);
        hisPtInOCSpectrum->Sumw2();
        TH2D* hisPtEtaInOCSpectrum = new TH2D(Form("hisPtEtaInOCSpectrum-%d", iCent), Form("%s: #it{p}_{T} spectrum, outside jet cones, c. %s;#it{p}_{T}^{%s} (GeV/#it{c});#eta", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), V0Symbol[V0Type].Data()), iNBinsPtV0InJet, dBinsPtV0InJet, iNBinsEtaInJet, dBinsEtaInJet);
        TH1D* hisEffInOC;

        TH1D* hisPtInNJNorm = GetHistogram1D(listStd, Form(hisPtInNJNormName.Data(), iCent));
        if(!hisPtInNJNorm)
          return;
        TH1D* hisPtInOCNorm = GetHistogram1D(listStd, hisPtInOCNormName.Data());
        if(!hisPtInOCNorm)
          return;
        Float_t fNEventsOK = hisPtInNJNorm->GetBinContent(9); // 9 after 2016-03-18, 3 before
        Float_t fNEventsWJets = hisPtInNJNorm->GetBinContent(12); // 12 after 2016-03, 6 before
        Float_t fNEventsNoJets = fNEventsOK - fNEventsWJets;
        Float_t fAreaExcluded = hisPtInOCNorm->GetBinContent(iCent + 1);
        Float_t fNConesRnd = hisPtInRCNorm->GetBinContent(iCent + 1);
        Float_t fNConesMed = hisPtInMCCNorm->GetBinContent(iCent + 1);
        printf("Number of events: EvtHisto: %.0f, EvtOK: %.0f, EvtWJets: %.0f\n", eventsCent[iCent], fNEventsOK, fNEventsWJets);

        TH2D* hisPtEtaPurity = 0;
        Int_t iBinMFirst = 0, iBinMLast = -1;
        if(iCorrEff && !iCorrEffInclusive)
        {
          dirOutSpectra->cd();
          hisPtEtaPurity = (TH2D*)dirOutSpectra->Get(Form("fh2PtEtaIncl%s_C%d-Purity", V0Name[V0Type].Data(), iCent));
          if(!hisPtEtaPurity)
          {
            printf("Loading purity failed. (Error)\n");
            return;
          }
          if(iCorrEffPur)
          {
            iBinMFirst = spaPtInRC->GetAxis(0)->FindBin(dSigMin + dEpsilon);
            iBinMLast = spaPtInRC->GetAxis(0)->FindBin(dSigMax - dEpsilon);
            spaPtInRC->GetAxis(0)->SetRange(iBinMFirst, iBinMLast);
            spaPtNoJet->GetAxis(0)->SetRange(iBinMFirst, iBinMLast);
            spaPtInMCC->GetAxis(0)->SetRange(iBinMFirst, iBinMLast);
            spaPtOutJet->GetAxis(0)->SetRange(iBinMFirst, iBinMLast);
          }
          // RC
          TH2D* hisPtEtaInRCSpectrumRaw = (TH2D*)spaPtInRC->Projection(2, 1, "e");
          if(!RebinHistogram2D(hisPtEtaInRCSpectrumRaw, hisPtEtaInRCSpectrum))
            return;
          TH2D* hisPtEtaInRCSpectrumSig = MultiplyHistograms2D(hisPtEtaInRCSpectrum, hisPtEtaPurity);
          if(iCorrEffPur)
            hisEffInRC = GetScaledEfficiency(hisPtEtaInRCSpectrumSig, hisEffPtEtaInclBase[iCent], 1);
          else
            hisEffInRC = GetScaledEfficiency(hisPtEtaInRCSpectrum, hisEffPtEtaInclBase[iCent], 1);
          if(!hisEffInRC)
            return;
          hisPtEtaInRCSpectrumRaw->Write(Form("RC-Raw-%d", iCent));
          hisPtEtaInRCSpectrum->Write(Form("RC-RawRebin-%d", iCent));
          hisPtEtaInRCSpectrumSig->Write(Form("RC-Sig-%d", iCent));
          hisEffInRC->Write(Form("RC-Eff-%d", iCent));
          // NJ
          TH2D* hisPtEtaInNJSpectrumRaw = (TH2D*)spaPtNoJet->Projection(2, 1, "e");
          if(!RebinHistogram2D(hisPtEtaInNJSpectrumRaw, hisPtEtaInNJSpectrum))
            return;
          TH2D* hisPtEtaInNJSpectrumSig = MultiplyHistograms2D(hisPtEtaInNJSpectrum, hisPtEtaPurity);
          if(iCorrEffPur)
            hisEffInNJ = GetScaledEfficiency(hisPtEtaInNJSpectrumSig, hisEffPtEtaInclBase[iCent], 1);
          else
            hisEffInNJ = GetScaledEfficiency(hisPtEtaInNJSpectrum, hisEffPtEtaInclBase[iCent], 1);
          if(!hisEffInNJ)
            return;
          hisPtEtaInNJSpectrumRaw->Write(Form("NJ-Raw-%d", iCent));
          hisPtEtaInNJSpectrum->Write(Form("NJ-RawRebin-%d", iCent));
          hisPtEtaInNJSpectrumSig->Write(Form("NJ-Sig-%d", iCent));
          hisEffInNJ->Write(Form("NJ-Eff-%d", iCent));
          // MCC
          TH2D* hisPtEtaInMCCSpectrumRaw = (TH2D*)spaPtInMCC->Projection(2, 1, "e");
          if(!RebinHistogram2D(hisPtEtaInMCCSpectrumRaw, hisPtEtaInMCCSpectrum))
            return;
          TH2D* hisPtEtaInMCCSpectrumSig = MultiplyHistograms2D(hisPtEtaInMCCSpectrum, hisPtEtaPurity);
          if(iCorrEffPur)
            hisEffInMCC = GetScaledEfficiency(hisPtEtaInMCCSpectrumSig, hisEffPtEtaInclBase[iCent], 1);
          else
            hisEffInMCC = GetScaledEfficiency(hisPtEtaInMCCSpectrum, hisEffPtEtaInclBase[iCent], 1);
          if(!hisEffInMCC)
            return;
          // OC
          TH2D* hisPtEtaInOCSpectrumRaw = (TH2D*)spaPtOutJet->Projection(2, 1, "e");
          if(!RebinHistogram2D(hisPtEtaInOCSpectrumRaw, hisPtEtaInOCSpectrum))
            return;
          TH2D* hisPtEtaInOCSpectrumSig = MultiplyHistograms2D(hisPtEtaInOCSpectrum, hisPtEtaPurity);
          if(iCorrEffPur)
            hisEffInOC = GetScaledEfficiency(hisPtEtaInOCSpectrumSig, hisEffPtEtaInclBase[iCent], 1);
          else
            hisEffInOC = GetScaledEfficiency(hisPtEtaInOCSpectrum, hisEffPtEtaInclBase[iCent], 1);
          if(!hisEffInOC)
            return;

          spaPtInRC->GetAxis(0)->SetRange(0, 0);
          spaPtNoJet->GetAxis(0)->SetRange(0, 0);
          spaPtInMCC->GetAxis(0)->SetRange(0, 0);
          spaPtOutJet->GetAxis(0)->SetRange(0, 0);
        }

        for(Int_t iV0 = iBinPtInJetsFirst - 1; iV0 <= iBinPtInJetsLast - 1; iV0++) // y (V0 pt bins)
        {
          Int_t iV0PtFirstBinRnd = spaPtInRC->GetAxis(1)->FindBin(dBinsPtV0InJet[iV0] + dEpsilon);
          Int_t iV0PtLastBinRnd = spaPtInRC->GetAxis(1)->FindBin(dBinsPtV0InJet[iV0 + 1] - dEpsilon);
          spaPtInRC->GetAxis(1)->SetRange(iV0PtFirstBinRnd, iV0PtLastBinRnd);

          // get invariant mass spectrum for a range of bins in y
          printf("iInRC: Making projection V0: %.1f-%.1f, bins: V0 %d-%d\n", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1], iV0PtFirstBinRnd, iV0PtLastBinRnd);

          // bin counting
          TH1D* hisMassInRC = (TH1D*)spaPtInRC->Projection(0, "e");
          if(!hisMassInRC)
            return;
          if(iNRebin > 1)
            hisMassInRC->Rebin(iNRebin);
          hisMassInRC->SetName(Form("PtInRC%s_C%d-%d", V0Name[V0Type].Data(), iCent, iV0));

          TGraphErrors* grPtInRCAll = MakeGraphErrors(hisMassInRC, Form("%.1f-%.1f", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1]), iMyColors[iV0 % iNMyColors], iMyMarkersEmpty[iV0 % iNMyMarkersEmpty], 0.5);
          mgrPtInRCAll->Add(grPtInRCAll);

          if(iUEsubtraction)
          {
            // V0 in random cones
            Bool_t bStatusBCRnd = kTRUE;
            BinCounterObject* binCountRnd = new BinCounterObject(hisMassInRC, fMassV0[V0Type], fSigmaV0[V0Type]);
            binCountRnd->SetVerbose(bVerboseBC);
            binCountRnd->SetDegreePolMax(iDegreePolSideBands);
            bStatusBCRnd &= binCountRnd->SetXLimits(dBgOutMin, dBgOutMax);
            if(bFixedBC)
              bStatusBCRnd &= binCountRnd->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax);
            else
              bStatusBCRnd &= binCountRnd->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
            binCountRnd->SetFitOptionSideBands(sOptionFitSB.Data());
//                  binCountRnd->Plot("canBinCounter0");
            bStatusBCRnd &= binCountRnd->EstimateParameters();
//                  binCountRnd->Plot("canBinCounter1");
            bStatusBCRnd &= binCountRnd->Run();
            binCountRnd->Plot("canBinCounter2");
            if(!bStatusBCRnd)
            {
              printf("Something wrong with bin counting\n");
            }
            else
            {
              Double_t dSignalRndErr = 0;
              Double_t dSignalRnd = binCountRnd->GetSignalAndError(&dSignalRndErr);
              hisPtInRCSpectrum->SetBinContent(iV0 + 1, dSignalRnd);
              hisPtInRCSpectrum->SetBinError(iV0 + 1, dSignalRndErr);
              printf("Signal: %f +- %f\n", dSignalRnd, dSignalRndErr);
              printf("Mean: %f +- %f\n", binCountRnd->GetMean(), binCountRnd->GetMeanError());
              printf("Sigma: %f +- %f\n", binCountRnd->GetSigma(), binCountRnd->GetSigmaError());
              TF1* funFitBgInRC = binCountRnd->GetFunctionSideBands();
              TGraph* grFunPtInRCBg = new TGraph(funFitBgInRC);
              grFunPtInRCBg->SetTitle(Form("%.1f-%.1f", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1]));
              grFunPtInRCBg->SetLineColor(iMyColors[iV0 % iNMyColors]);
              grFunPtInRCBg->SetLineStyle(2);
              grFunPtInRCBg->SetFillStyle(0);
              mgrPtInRCBg->Add(grFunPtInRCBg, "l");
            }
            delete binCountRnd;
          }

          // get invariant mass spectrum for a range of bins in y
          printf("iInMCC: Making projection V0: %.1f-%.1f, bins: V0 %d-%d\n", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1], iV0PtFirstBinRnd, iV0PtLastBinRnd);

          // bin counting
          spaPtInMCC->GetAxis(1)->SetRange(iV0PtFirstBinRnd, iV0PtLastBinRnd);
          TH1D* hisMassInMCC = (TH1D*)spaPtInMCC->Projection(0, "e");
          if(!hisMassInMCC)
            return;
          if(iNRebin > 1)
            hisMassInMCC->Rebin(iNRebin);
          hisMassInMCC->SetName(Form("PtInMCC%s_C%d-%d", V0Name[V0Type].Data(), iCent, iV0));

          if(iUEsubtraction && bIsPbPb)
          {
            // V0 in median-cluster cones
            Bool_t bStatusBCMed = kTRUE;
            BinCounterObject* binCountMed = new BinCounterObject(hisMassInMCC, fMassV0[V0Type], fSigmaV0[V0Type]);
            binCountMed->SetVerbose(bVerboseBC);
            binCountMed->SetDegreePolMax(iDegreePolSideBands);
            bStatusBCMed &= binCountMed->SetXLimits(dBgOutMin, dBgOutMax);
            if(bFixedBC)
              bStatusBCMed &= binCountMed->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax);
            else
              bStatusBCMed &= binCountMed->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
            binCountMed->SetFitOptionSideBands(sOptionFitSB.Data());
//                  binCountMed->Plot("canBinCounter0");
            bStatusBCMed &= binCountMed->EstimateParameters();
//                  binCountMed->Plot("canBinCounter1");
            bStatusBCMed &= binCountMed->Run();
            binCountMed->Plot("canBinCounter2");
            if(!bStatusBCMed)
            {
              printf("Something wrong with bin counting\n");
            }
            else
            {
              Double_t dSignalMedErr = 0;
              Double_t dSignalMed = binCountMed->GetSignalAndError(&dSignalMedErr);
              hisPtInMCCSpectrum->SetBinContent(iV0 + 1, dSignalMed);
              hisPtInMCCSpectrum->SetBinError(iV0 + 1, dSignalMedErr);
            }
            delete binCountMed;
          }

          // get invariant mass spectrum for a range of bins in y
          printf("iNoJet: Making projection V0: %.1f-%.1f, bins: V0 %d-%d\n", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1], iV0PtFirstBinRnd, iV0PtLastBinRnd);

          // bin counting
          spaPtNoJet->GetAxis(1)->SetRange(iV0PtFirstBinRnd, iV0PtLastBinRnd);
          TH1D* hisMassInNJ = (TH1D*)spaPtNoJet->Projection(0, "e");
          if(!hisMassInNJ)
            return;
          if(iNRebin > 1)
            hisMassInNJ->Rebin(iNRebin);
          hisMassInNJ->SetName(Form("PtInNJ%s_C%d-%d", V0Name[V0Type].Data(), iCent, iV0));

          TGraphErrors* grPtInNJAll = MakeGraphErrors(hisMassInNJ, Form("%.1f-%.1f", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1]), iMyColors[iV0 % iNMyColors], iMyMarkersEmpty[iV0 % iNMyMarkersEmpty], 0.5);
          mgrPtInNJAll->Add(grPtInNJAll);

          if(iUEsubtraction)
          {
            // V0 in jet-less events
            Bool_t bStatusBCInNJ = kTRUE;
            BinCounterObject* binCountInNJ = new BinCounterObject(hisMassInNJ, fMassV0[V0Type], fSigmaV0[V0Type]);
            binCountInNJ->SetVerbose(bVerboseBC);
            binCountInNJ->SetDegreePolMax(iDegreePolSideBands);
            bStatusBCInNJ &= binCountInNJ->SetXLimits(dBgOutMin, dBgOutMax);
            if(bFixedBC)
              bStatusBCInNJ &= binCountInNJ->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax);
            else
              bStatusBCInNJ &= binCountInNJ->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
            binCountInNJ->SetFitOptionSideBands(sOptionFitSB.Data());
//                  binCountInNJ->Plot("canBinCounter0");
            bStatusBCInNJ &= binCountInNJ->EstimateParameters();
//                  binCountInNJ->Plot("canBinCounter1");
            bStatusBCInNJ &= binCountInNJ->Run();
            binCountInNJ->Plot("canBinCounter2");
            if(!bStatusBCInNJ)
            {
              printf("Something wrong with bin counting\n");
            }
            else
            {
              Double_t dSignalInNJErr = 0;
              Double_t dSignalInNJ = binCountInNJ->GetSignalAndError(&dSignalInNJErr);
              hisPtInNJSpectrum->SetBinContent(iV0 + 1, dSignalInNJ);
              hisPtInNJSpectrum->SetBinError(iV0 + 1, dSignalInNJErr);
              printf("Signal: %f +- %f\n", dSignalInNJ, dSignalInNJErr);
              printf("Mean: %f +- %f\n", binCountInNJ->GetMean(), binCountInNJ->GetMeanError());
              printf("Sigma: %f +- %f\n", binCountInNJ->GetSigma(), binCountInNJ->GetSigmaError());
              TF1* funFitBgInNJ = binCountInNJ->GetFunctionSideBands();
              TGraph* grFunPtInNJBg = new TGraph(funFitBgInNJ);
              grFunPtInNJBg->SetTitle(Form("%.1f-%.1f", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1]));
              grFunPtInNJBg->SetLineColor(iMyColors[iV0 % iNMyColors]);
              grFunPtInNJBg->SetLineStyle(2);
              grFunPtInNJBg->SetFillStyle(0);
              mgrPtInNJBg->Add(grFunPtInNJBg, "l");
            }
            delete binCountInNJ;
          }

          // get invariant mass spectrum for a range of bins in y
          printf("iOutsideJet: Making projection V0: %.1f-%.1f, bins: V0 %d-%d\n", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1], iV0PtFirstBinRnd, iV0PtLastBinRnd);

          // bin counting
          spaPtOutJet->GetAxis(1)->SetRange(iV0PtFirstBinRnd, iV0PtLastBinRnd);
          TH1D* hisMassInOC = (TH1D*)spaPtOutJet->Projection(0, "e");
          if(!hisMassInOC)
            return;
          if(iNRebin > 1)
            hisMassInOC->Rebin(iNRebin);
          hisMassInOC->SetName(Form("PtInOC%s_C%d-%d", V0Name[V0Type].Data(), iCent, iV0));

          TGraphErrors* grPtInOCAll = MakeGraphErrors(hisMassInOC, Form("%.1f-%.1f", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1]), iMyColors[iV0 % iNMyColors], iMyMarkersEmpty[iV0 % iNMyMarkersEmpty], 0.5);
          mgrPtInOCAll->Add(grPtInOCAll);

          if(iUEsubtraction)
          {
            // V0 outside cones
            Bool_t bStatusBCInOC = kTRUE;
            BinCounterObject* binCountInOC = new BinCounterObject(hisMassInOC, fMassV0[V0Type], fSigmaV0[V0Type]);
            binCountInOC->SetVerbose(bVerboseBC);
            binCountInOC->SetDegreePolMax(iDegreePolSideBands);
            bStatusBCInOC &= binCountInOC->SetXLimits(dBgOutMin, dBgOutMax);
            if(bFixedBC)
              bStatusBCInOC &= binCountInOC->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax);
            else
              bStatusBCInOC &= binCountInOC->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
            binCountInOC->SetFitOptionSideBands(sOptionFitSB.Data());
//                  binCountInOC->Plot("canBinCounter0");
            bStatusBCInOC &= binCountInOC->EstimateParameters();
//                  binCountInOC->Plot("canBinCounter1");
            bStatusBCInOC &= binCountInOC->Run();
            binCountInOC->Plot("canBinCounter2");
            if(!bStatusBCInOC)
            {
              printf("Something wrong with bin counting\n");
            }
            else
            {
              Double_t dSignalInOCErr = 0;
              Double_t dSignalInOC = binCountInOC->GetSignalAndError(&dSignalInOCErr);
              hisPtInOCSpectrum->SetBinContent(iV0 + 1, dSignalInOC);
              hisPtInOCSpectrum->SetBinError(iV0 + 1, dSignalInOCErr);
              printf("Signal: %f +- %f\n", dSignalInOC, dSignalInOCErr);
              printf("Mean: %f +- %f\n", binCountInOC->GetMean(), binCountInOC->GetMeanError());
              printf("Sigma: %f +- %f\n", binCountInOC->GetSigma(), binCountInOC->GetSigmaError());
              TF1* funFitBgInOC = binCountInOC->GetFunctionSideBands();
              TGraph* grFunPtInOCBg = new TGraph(funFitBgInOC);
              grFunPtInOCBg->SetTitle(Form("%.1f-%.1f", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1]));
              grFunPtInOCBg->SetLineColor(iMyColors[iV0 % iNMyColors]);
              grFunPtInOCBg->SetLineStyle(2);
              grFunPtInOCBg->SetFillStyle(0);
              mgrPtInOCBg->Add(grFunPtInOCBg, "l");
            }
            delete binCountInOC;
          }
        }
        if(iCorrEff)
        {
          printf("Correcting %s for efficiency\n", hisPtInRCSpectrum->GetName());
          if(iCorrEffInclusive)
            hisPtInRCSpectrum = (TH1D*)DivideHistograms1D(hisPtInRCSpectrum, hisEffPtInJetsExt[iCent]);
          else
            hisPtInRCSpectrum = (TH1D*)DivideHistograms1D(hisPtInRCSpectrum, hisEffInRC);
          if(!hisPtInRCSpectrum)
            return;
          if(bIsPbPb)
          {
            printf("Correcting %s for efficiency\n", hisPtInMCCSpectrum->GetName());
            if(iCorrEffInclusive)
              hisPtInMCCSpectrum = (TH1D*)DivideHistograms1D(hisPtInMCCSpectrum, hisEffPtInJetsExt[iCent]);
            else
              hisPtInMCCSpectrum = (TH1D*)DivideHistograms1D(hisPtInMCCSpectrum, hisEffInMCC);
            if(!hisPtInMCCSpectrum)
              return;
          }
          printf("Correcting %s for efficiency\n", hisPtInNJSpectrum->GetName());
          if(iCorrEffInclusive)
            hisPtInNJSpectrum = (TH1D*)DivideHistograms1D(hisPtInNJSpectrum, hisEffPtInJetsExt[iCent]);
          else
            hisPtInNJSpectrum = (TH1D*)DivideHistograms1D(hisPtInNJSpectrum, hisEffInNJ);
          if(!hisPtInNJSpectrum)
            return;
          printf("Correcting %s for efficiency\n", hisPtInOCSpectrum->GetName());
          if(iCorrEffInclusive)
            hisPtInOCSpectrum = (TH1D*)DivideHistograms1D(hisPtInOCSpectrum, hisEffPtInJetsExt[iCent]);
          else
            hisPtInOCSpectrum = (TH1D*)DivideHistograms1D(hisPtInOCSpectrum, hisEffInOC);
          if(!hisPtInOCSpectrum)
            return;
        }
        dirOutSpectra->cd();

        printf("Normalization of UE spectra\n");

        if(iUEsubtraction)
        {
          // normalize per area unit
          Float_t fNormRC = fNConesRnd;
          if(bNormPerArea)
            fNormRC *= fAreaCone;
          hisPtInRCSpectrum = DivideHistogram(hisPtInRCSpectrum, fNormRC, 0, 1);
          if(!hisPtInRCSpectrum)
            return;
          hisPtInRCSpectrum->Write(Form("fh1PtInRC%s_C%d", V0Name[V0Type].Data(), iCent));

          // normalize per area unit
          Float_t fNormMCC = fNConesMed;
          if(bNormPerArea)
            fNormMCC *= fAreaCone;
          if(bIsPbPb)
            hisPtInMCCSpectrum = DivideHistogram(hisPtInMCCSpectrum, fNormMCC, 0, 1);
          if(!hisPtInMCCSpectrum)
            return;
          hisPtInMCCSpectrum->Write(Form("fh1PtInMCC%s_C%d", V0Name[V0Type].Data(), iCent));

          // normalize per area unit
          Float_t fNormNJ = fNEventsNoJets * fNormalisationV0Area;
          if(!fNormNJ)
            return;
          if(!bNormPerArea)
            fNormNJ /= fAreaCone;
          if(fNormNJ)
            hisPtInNJSpectrum = DivideHistogram(hisPtInNJSpectrum, fNormNJ, 0, 1);
          if(!hisPtInNJSpectrum)
            return;
          hisPtInNJSpectrum->Write(Form("fh1PtInNJ%s_C%d", V0Name[V0Type].Data(), iCent));

          // normalize per area unit
          Float_t fNormOC = fNEventsWJets * fNormalisationV0Area - fAreaExcluded;
          if(!bNormPerArea)
            fNormOC /= fAreaCone;
          hisPtInOCSpectrum = DivideHistogram(hisPtInOCSpectrum, fNormOC, 0, 1);
          if(!hisPtInOCSpectrum)
            return;
          hisPtInOCSpectrum->Write(Form("fh1PtInOC%s_C%d", V0Name[V0Type].Data(), iCent));
        }

        // spectrum for UE subtraction
        TH1D* hisPtInBulkSpectrumSub;

        // reference spectrum for comparison
        TH1D* hisPtInBulkSpectrumRef;
        TString sBulkRef = "NJ";
        hisPtInBulkSpectrumRef = hisPtInNJSpectrum;
        // spectrum for UE subtraction
        const Int_t iNUE = 5;
        TH1D* hisPtInBulkSpectrum[iNUE];
        TH1D* hisPtInJetsSpectrumUE[iNUE];
        TString sPtInBulkName[iNUE];
        hisPtInBulkSpectrum[0] = hisPtInNJSpectrum;
        sPtInBulkName[0] = "NJ";
        hisPtInBulkSpectrum[1] = hisPtInOCSpectrum;
        sPtInBulkName[1] = "OC";
        hisPtInBulkSpectrum[2] = hisPtInRCSpectrum;
        sPtInBulkName[2] = "RC";
        hisPtInBulkSpectrum[3] = hisPtInMCCSpectrum;
        sPtInBulkName[3] = "MCC";

        TGraphErrors* grPtInRCSpectrum = MakeGraphErrors(hisPtInRCSpectrum, Form("RC %s", GetCentBinLabel(iCent).Data()), iMyColors[iCent], iMyMarkersEmpty[iCent]);
        mgrPtInRC->Add(grPtInRCSpectrum);
        TGraphErrors* grPtInRCSpectrumP = MakeGraphErrors(hisPtInRCSpectrum, Form("RC"), iMyColors[3], iMyMarkersEmpty[iCent]);
        mgrPtInBulk->Add(grPtInRCSpectrumP);
        TH1D* hisPtInRCSpectrumComp = DivideHistograms1D(hisPtInRCSpectrum, hisPtInBulkSpectrumRef);
        TGraphErrors* grPtInRCSpectrumComp = MakeGraphErrors(hisPtInRCSpectrumComp, Form("RC/%s", sBulkRef.Data()), iMyColors[3], iMyMarkersEmpty[iCent]);
        mgrPtInBulkCompare->Add(grPtInRCSpectrumComp);

        TGraphErrors* grPtInMCCSpectrumP = MakeGraphErrors(hisPtInMCCSpectrum, Form("MCC"), iMyColors[2], iMyMarkersEmpty[iCent]);
        mgrPtInBulk->Add(grPtInMCCSpectrumP);
        TH1D* hisPtInMCCSpectrumComp = DivideHistograms1D(hisPtInMCCSpectrum, hisPtInBulkSpectrumRef);
        TGraphErrors* grPtInMCCSpectrumComp = MakeGraphErrors(hisPtInMCCSpectrumComp, Form("MCC/%s", sBulkRef.Data()), iMyColors[2], iMyMarkersEmpty[iCent]);
        mgrPtInBulkCompare->Add(grPtInMCCSpectrumComp);

        TGraphErrors* grPtInNJSpectrum = MakeGraphErrors(hisPtInNJSpectrum, Form("NJ"), iMyColors[iCent], iMyMarkersEmpty[iCent]);
        mgrPtInNJ->Add(grPtInNJSpectrum);
        TGraphErrors* grPtInNJSpectrumP = MakeGraphErrors(hisPtInNJSpectrum, Form("NJ"), iMyColors[0], iMyMarkersEmpty[iCent]);
        mgrPtInBulk->Add(grPtInNJSpectrumP);
        TH1D* hisPtInNJSpectrumComp = DivideHistograms1D(hisPtInNJSpectrum, hisPtInBulkSpectrumRef);
        TGraphErrors* grPtInNJSpectrumComp = MakeGraphErrors(hisPtInNJSpectrumComp, Form("NJ/%s", sBulkRef.Data()), iMyColors[0], iMyMarkersEmpty[iCent]);
        mgrPtInBulkCompare->Add(grPtInNJSpectrumComp);

        TGraphErrors* grPtInOCSpectrum = MakeGraphErrors(hisPtInOCSpectrum, Form("OC"), iMyColors[iCent], iMyMarkersEmpty[iCent]);
        mgrPtInOC->Add(grPtInOCSpectrum);
        TGraphErrors* grPtInOCSpectrumP = MakeGraphErrors(hisPtInOCSpectrum, Form("OC"), iMyColors[4], iMyMarkersEmpty[iCent]);
        mgrPtInBulk->Add(grPtInOCSpectrumP);
        TH1D* hisPtInOCSpectrumComp = DivideHistograms1D(hisPtInOCSpectrum, hisPtInBulkSpectrumRef);
        TGraphErrors* grPtInOCSpectrumComp = MakeGraphErrors(hisPtInOCSpectrumComp, Form("OC/%s", sBulkRef.Data()), iMyColors[4], iMyMarkersEmpty[iCent]);
        mgrPtInBulkCompare->Add(grPtInOCSpectrumComp);
//              printf("Fitting OC\n");
//              TFitResultPtr result = hisPtInOCSpectrumComp->Fit("pol1","SLI","",2,10);
//              if (!BinCounterObject::IsFitOK(result))
//                return;

        Int_t iJetSel = 1;
        for(Int_t iJet = iJetMin; iJet <= iJetMax; iJet++)  // z (jet pt bins)
//              for (Int_t iJet = iJetSel; iJet <= iJetSel; iJet++) // z (jet pt bins)
        {
          printf("Jets in pt range (rebin): %f +- %f\n", dArrayJetNorm[iJet], dArrayJetNormErr[iJet]);
          if(!dArrayJetNorm[iJet])
            continue;

          printf("Check jets OK\n");
          TCanvas* canPtInJetsAll = new TCanvas(Form("canInJetsAll%s-%d-%d", V0Name[V0Type].Data(), iCent, iJet), "", iCanWidth, iCanHeight);
          TMultiGraph* mgrPtInJetsAll = new TMultiGraph();
          TCanvas* canPtInJetsBg = new TCanvas(Form("canInJetsBg%s-%d-%d", V0Name[V0Type].Data(), iCent, iJet), "", iCanWidth, iCanHeight);
          TMultiGraph* mgrPtInJetsBg = new TMultiGraph();

          TCanvas* canPtInPCAll = new TCanvas(Form("canInPCAll%s-%d-%d", V0Name[V0Type].Data(), iCent, iJet), "", iCanWidth, iCanHeight);
          TMultiGraph* mgrPtInPCAll = new TMultiGraph();
          TCanvas* canPtInPCBg = new TCanvas(Form("canInPCBg%s-%d-%d", V0Name[V0Type].Data(), iCent, iJet), "", iCanWidth, iCanHeight);
          TMultiGraph* mgrPtInPCBg = new TMultiGraph();

          TCanvas* canBulkOverJC = new TCanvas(Form("canBulkOverJC_%d", iCent), "", iCanWidth, iCanHeight);
          TMultiGraph* mgrBulkOverJC = new TMultiGraph();
          TCanvas* canBulkOverSignal = new TCanvas(Form("canBulkOverSignal_%d", iCent), "", iCanWidth, iCanHeight);
          TMultiGraph* mgrBulkOverSignal = new TMultiGraph();

          Int_t iJetPtFirstBin = spaPtInJets->GetAxis(3)->FindBin(dBinsPtJet[iJet] + dEpsilon);
          Int_t iJetPtLastBin;
          if(bOpenPtBins)
            iJetPtLastBin = spaPtInJets->GetAxis(3)->GetNbins() + 1;
          else
            iJetPtLastBin = spaPtInJets->GetAxis(3)->FindBin(dBinsPtJet[iJet + 1] - dEpsilon);
          spaPtInJets->GetAxis(3)->SetRange(iJetPtFirstBin, iJetPtLastBin);
          spaPtInPC->GetAxis(3)->SetRange(iJetPtFirstBin, iJetPtLastBin);

          // pt spectrum in jets
          TH1D* hisPtInJetsSpectrum = new TH1D(Form("hisPtInJetsSpectrum-%d-%d", iCent, iJet), Form("%s: #it{p}_{T} spectrum, in jets, c. %s;#it{p}_{T}^{%s} (GeV/#it{c});#frac{1}{#it{N}_{jet}} #frac{d#it{N}^{%s}}{d#it{p}_{T}^{%s}} (#it{c} GeV^{#minus1})", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data()), iNBinsPtV0InJet, dBinsPtV0InJet);
          hisPtInJetsSpectrum->Sumw2();
          TH2D* hisPtEtaInJetsSpectrum = new TH2D(Form("hisPtEtaInJetsSpectrum-%d-%d", iCent, iJet), Form("%s: #it{p}_{T} spectrum, in jets, c. %s;#it{p}_{T}^{%s} (GeV/#it{c});#frac{1}{#it{N}_{jet}} #frac{d#it{N}^{%s}}{d#it{p}_{T}^{%s}} (#it{c} GeV^{#minus1})", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data()), iNBinsPtV0InJet, dBinsPtV0InJet, iNBinsEtaInJet, dBinsEtaInJet);
          TH1D* hisEffInJC;
          // pt spectrum in perp. cones
          TH1D* hisPtInPCSpectrum = new TH1D(Form("hisPtInPCSpectrum-%d-%d", iCent, iJet), Form("%s: #it{p}_{T} spectrum, in perp. cones, c. %s;#it{p}_{T}^{%s} (GeV/#it{c});#frac{1}{#it{N}_{jet}} #frac{d#it{N}^{%s}}{d#it{p}_{T}^{%s}} (#it{c} GeV^{#minus1})", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data()), iNBinsPtV0InJet, dBinsPtV0InJet);
          hisPtInPCSpectrum->Sumw2();
          TH2D* hisPtEtaInPCSpectrum = new TH2D(Form("hisPtEtaInPCSpectrum-%d-%d", iCent, iJet), Form("%s: #it{p}_{T} spectrum, in perp. cones, c. %s;#it{p}_{T}^{%s} (GeV/#it{c});#frac{1}{#it{N}_{jet}} #frac{d#it{N}^{%s}}{d#it{p}_{T}^{%s}} (#it{c} GeV^{#minus1})", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data(), V0Symbol[V0Type].Data()), iNBinsPtV0InJet, dBinsPtV0InJet, iNBinsEtaInJet, dBinsEtaInJet);
          TH1D* hisEffInPC;

          TH1D* hisPurityInJC = new TH1D(Form("hisPurityInJC-%d-%d", iCent, iJet), Form("%s signal purity in JC, %s;#it{p}_{T} (GeV/c)", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()), iNBinsPtV0InJet, dBinsPtV0InJet);

          if(iCorrEff && !iCorrEffInclusive)
          {
//                  dirOutSpectra->cd();
            spaPtInPC->GetAxis(1)->SetRange(0, 0);
            if(iCorrEffPur)
              spaPtInPC->GetAxis(0)->SetRange(iBinMFirst, iBinMLast);
            TH2D* hisPtEtaInPCSpectrumRaw = (TH2D*)spaPtInPC->Projection(2, 1, "e");
            if(!RebinHistogram2D(hisPtEtaInPCSpectrumRaw, hisPtEtaInPCSpectrum))
              return;
            TH2D* hisPtEtaInPCSpectrumSig = MultiplyHistograms2D(hisPtEtaInPCSpectrum, hisPtEtaPurity);
            if(iCorrEffPur)
              hisEffInPC = GetScaledEfficiency(hisPtEtaInPCSpectrumSig, hisEffPtEtaInclBase[iCent], 1);
            else
              hisEffInPC = GetScaledEfficiency(hisPtEtaInPCSpectrum, hisEffPtEtaInclBase[iCent], 1);
            if(!hisEffInPC)
              return;
            spaPtInPC->GetAxis(0)->SetRange(0, 0);

            spaPtInJets->GetAxis(1)->SetRange(0, 0);
            if(iCorrEffPur)
              spaPtInJets->GetAxis(0)->SetRange(iBinMFirst, iBinMLast);
            TH2D* hisPtEtaInJCSpectrumRaw = (TH2D*)spaPtInJets->Projection(2, 1, "e");
            if(!RebinHistogram2D(hisPtEtaInJCSpectrumRaw, hisPtEtaInJetsSpectrum))
              return;
            TH2D* hisPtEtaInJetsSpectrumSig = MultiplyHistograms2D(hisPtEtaInJetsSpectrum, hisPtEtaPurity);
            if(iCorrEffPur)
              hisEffInJC = GetScaledEfficiency(hisPtEtaInJetsSpectrumSig, hisEffPtEtaInclBase[iCent], 1);
            else
              hisEffInJC = GetScaledEfficiency(hisPtEtaInJetsSpectrum, hisEffPtEtaInclBase[iCent], 1);
            if(!hisEffInJC)
              return;
            spaPtInJets->GetAxis(0)->SetRange(0, 0);

            hisPtEtaInJCSpectrumRaw->Write(Form("JC-Raw-%d-%d", iCent, iJet));
            hisPtEtaInJetsSpectrum->Write(Form("JC-RawRebin-%d-%d", iCent, iJet));
            hisPtEtaInJetsSpectrumSig->Write(Form("JC-Sig-%d-%d", iCent, iJet));
            hisEffInJC->Write(Form("JC-Eff-%d-%d", iCent, iJet));

            TH1D* hisPtEtaInJetsSpectrumSigProj = (TH1D*)hisPtEtaInJetsSpectrumSig->ProjectionX(Form("%s_px", hisPtEtaInJetsSpectrumSig->GetName()), 0, -1, "e");
            hisPtEtaInJetsSpectrumSigProj->Write(Form("JC-SigProj-%d-%d", iCent, iJet));
          }

          for(Int_t iV0 = iBinPtInJetsFirst - 1; iV0 <= iBinPtInJetsLast - 1; iV0++) // y (V0 pt bins)
          {
            Int_t iV0PtFirstBin = spaPtInJets->GetAxis(1)->FindBin(dBinsPtV0InJet[iV0] + dEpsilon);
            Int_t iV0PtLastBin = spaPtInJets->GetAxis(1)->FindBin(dBinsPtV0InJet[iV0 + 1] - dEpsilon);
            spaPtInJets->GetAxis(1)->SetRange(iV0PtFirstBin, iV0PtLastBin);
            spaPtInPC->GetAxis(1)->SetRange(iV0PtFirstBin, iV0PtLastBin);

            // get invariant mass spectrum for a range of bins in y, z
            if(bOpenPtBins)
              printf("iPtInJets: Making projection V0: %.1f-%.1f, jet >%.0f, bins: V0 %d-%d, jet %d-%d\n", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1], dBinsPtJet[iJet], iV0PtFirstBin, iV0PtLastBin, iJetPtFirstBin, iJetPtLastBin);
            else
              printf("iPtInJets: Making projection V0: %.1f-%.1f, jet %.0f-%.0f, bins: V0 %d-%d, jet %d-%d\n", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1], dBinsPtJet[iJet], dBinsPtJet[iJet + 1], iV0PtFirstBin, iV0PtLastBin, iJetPtFirstBin, iJetPtLastBin);

            // bin counting
            TH1D* hisMassInJets = (TH1D*)spaPtInJets->Projection(0, "E");
            if(!hisMassInJets)
              return;
            if(iNRebin > 1)
              hisMassInJets->Rebin(iNRebin);
            hisMassInJets->SetName(Form("PtInJets%s_C%d-J%d-%d", V0Name[V0Type].Data(), iCent, iJet, iV0));
            TH1D* hisMassInPC = (TH1D*)spaPtInPC->Projection(0, "E");
            if(!hisMassInPC)
              return;
            if(iNRebin > 1)
              hisMassInPC->Rebin(iNRebin);
            hisMassInPC->SetName(Form("PtInPC%s_C%d-J%d-%d", V0Name[V0Type].Data(), iCent, iJet, iV0));

            TGraphErrors* grPtInJetsAll = MakeGraphErrors(hisMassInJets, Form("%.1f-%.1f", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1]), iMyColors[iV0 % iNMyColors], iMyMarkersEmpty[iV0 % iNMyMarkersEmpty], 0.5);
            mgrPtInJetsAll->Add(grPtInJetsAll);
            TGraphErrors* grPtInPCAll = MakeGraphErrors(hisMassInPC, Form("%.1f-%.1f", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1]), iMyColors[iV0 % iNMyColors], iMyMarkersEmpty[iV0 % iNMyMarkersEmpty], 0.5);
            mgrPtInPCAll->Add(grPtInPCAll);

            // V0 in jet cones
            Bool_t bStatusBC = kTRUE;
            BinCounterObject* binCount = new BinCounterObject(hisMassInJets, fMassV0[V0Type], fSigmaV0[V0Type]);
            binCount->SetVerbose(bVerboseBC);
            binCount->SetDegreePolMax(iDegreePolSideBands);
            bStatusBC &= binCount->SetXLimits(dBgOutMin, dBgOutMax);
            if(bFixedBC)
              bStatusBC &= binCount->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax);
            else
              bStatusBC &= binCount->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
            binCount->SetFitOptionSideBands(sOptionFitSB.Data());

            TString sTitle = Form("%s, in JC, c. %s, pT: %.1f-%.1f GeV/c", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1]);
//                      binCount->Plot("canBinCounter0",Form("%s, step %d",sTitle.Data(),0));
            bStatusBC &= binCount->EstimateParameters();
//                      binCount->Plot("canBinCounter1",Form("%s, step %d",sTitle.Data(),1));
            bStatusBC &= binCount->Run();
            binCount->Plot("canBinCounter2", Form("%s, step %d", sTitle.Data(), 2));
            if(!bStatusBC)
            {
              printf("Something wrong with bin counting\n");
            }
            else
            {
              Double_t dSignalErr = 0;
              Double_t dSignal = binCount->GetSignalAndError(&dSignalErr);
              hisPtInJetsSpectrum->SetBinContent(iV0 + 1, dSignal);
              hisPtInJetsSpectrum->SetBinError(iV0 + 1, dSignalErr);
              printf("Signal: %f +- %f\n", dSignal, dSignalErr);
              printf("Mean: %f +- %f\n", binCount->GetMean(), binCount->GetMeanError());
              printf("Sigma: %f +- %f\n", binCount->GetSigma(), binCount->GetSigmaError());
              Double_t dPurityErr = 0;
              Double_t dPurity = binCount->GetPurityAndError(&dPurityErr);
              hisPurityInJC->SetBinContent(iV0 + 1, dPurity);
              hisPurityInJC->SetBinError(iV0 + 1, dPurityErr);
              TF1* funFitBgInJets = binCount->GetFunctionSideBands();
              TGraph* grFunPtInJetsBg = new TGraph(funFitBgInJets);
              grFunPtInJetsBg->SetTitle(Form("%.1f-%.1f", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1]));
              grFunPtInJetsBg->SetLineColor(iMyColors[iV0 % iNMyColors]);
              grFunPtInJetsBg->SetLineStyle(2);
              grFunPtInJetsBg->SetFillStyle(0);
              mgrPtInJetsBg->Add(grFunPtInJetsBg, "l");
            }
            delete binCount;

            if(iUEsubtraction)
            {
              // V0 in perpendicular cones
              Bool_t bStatusBCPerp = kTRUE;
              BinCounterObject* binCountPerp = new BinCounterObject(hisMassInPC, fMassV0[V0Type], fSigmaV0[V0Type]);
              binCountPerp->SetVerbose(bVerboseBC);
              binCountPerp->SetDegreePolMax(iDegreePolSideBands);
              bStatusBCPerp &= binCountPerp->SetXLimits(dBgOutMin, dBgOutMax);
              if(bFixedBC)
                bStatusBCPerp &= binCountPerp->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax);
              else
                bStatusBCPerp &= binCountPerp->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
              binCountPerp->SetFitOptionSideBands(sOptionFitSB.Data());
//                      binCountPerp->Plot("canBinCounter0");
              bStatusBCPerp &= binCountPerp->EstimateParameters();
//                      binCountPerp->Plot("canBinCounter1");
              bStatusBCPerp &= binCountPerp->Run();
              binCountPerp->Plot("canBinCounter2");
              if(!bStatusBCPerp)
              {
                printf("Something wrong with bin counting\n");
              }
              else
              {
                Double_t dSignalPerpErr = 0;
                Double_t dSignalPerp = binCountPerp->GetSignalAndError(&dSignalPerpErr);
                hisPtInPCSpectrum->SetBinContent(iV0 + 1, dSignalPerp);
                hisPtInPCSpectrum->SetBinError(iV0 + 1, dSignalPerpErr);
                printf("Signal: %f +- %f\n", dSignalPerp, dSignalPerpErr);
                printf("Mean: %f +- %f\n", binCountPerp->GetMean(), binCountPerp->GetMeanError());
                printf("Sigma: %f +- %f\n", binCountPerp->GetSigma(), binCountPerp->GetSigmaError());
                TF1* funFitBgInPC = binCountPerp->GetFunctionSideBands();
                TGraph* grFunPtInPCBg = new TGraph(funFitBgInPC);
                grFunPtInPCBg->SetTitle(Form("%.1f-%.1f", dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1]));
                grFunPtInPCBg->SetLineColor(iMyColors[iV0 % iNMyColors]);
                grFunPtInPCBg->SetLineStyle(2);
                grFunPtInPCBg->SetFillStyle(0);
                mgrPtInPCBg->Add(grFunPtInPCBg, "l");
              }
              delete binCountPerp;
            }
          }
          TH1D* hisPtInJCSpectrum = (TH1D*)hisPtInJetsSpectrum->Clone(Form("%sNoSub", hisPtInJetsSpectrum->GetName()));
          if(iCorrEff && !iCorrEffInclusive)
          {
            TH1D* hisPtInJCSpectrumProj = (TH1D*)dirOutSpectra->Get(Form("JC-SigProj-%d-%d", iCent, iJet));
            if(!hisPtInJCSpectrumProj)
            {
              printf("Loading JC-SigProj failed. (Error)\n");
              return;
            }
            TH1D* hisPtInJCRatio = DivideHistograms1D(hisPtInJCSpectrumProj, hisPtInJCSpectrum);
            hisPtInJCRatio->Write(Form("JC-SigRatio-%d-%d", iCent, iJet));
          }
          // scale to number of jets (2 PC per 1 jet)
          hisPtInPCSpectrum->Scale(0.5);
          if(!bYieldsOnly)
          {
            // normalize per area unit and by bin width
            Float_t fNormJC = dArrayJetNorm[iJet];
            Float_t fNormJCErr = dArrayJetNormErr[iJet];
            if(bNormPerArea)
            {
              fNormJC *= fAreaCone;
              fNormJCErr *= fAreaCone;
            }
            hisPtInJCSpectrum = DivideHistogram(hisPtInJCSpectrum, fNormJC, fNormJCErr, 1);
            if(!hisPtInJCSpectrum)
              return;
            hisPtInJetsSpectrum = DivideHistogram(hisPtInJetsSpectrum, fNormJC, fNormJCErr, 1);
            if(!hisPtInJetsSpectrum)
              return;
            hisPtInPCSpectrum = DivideHistogram(hisPtInPCSpectrum, fNormJC, fNormJCErr, 1);
            if(!hisPtInPCSpectrum)
              return;
          }
          dirOutSpectra->cd();
          hisPtInJetsSpectrum->Write(Form("fh1PtInJets%s_C%d-J%d_Raw", V0Name[V0Type].Data(), iCent, iJet));
          // efficiency correction
          if(iCorrEff)
          {
            printf("Correcting %s for efficiency\n", hisPtInJCSpectrum->GetName());
            if(iCorrEffInclusive)
              hisPtInJCSpectrum = DivideHistograms1D(hisPtInJCSpectrum, hisEffPtInJetsExt[iCent]);
            else
              hisPtInJCSpectrum = DivideHistograms1D(hisPtInJCSpectrum, hisEffInJC);
            if(!hisPtInJCSpectrum)
              return;
            printf("Correcting %s for efficiency\n", hisPtInJetsSpectrum->GetName());
            if(iCorrEffInclusive)
              hisPtInJetsSpectrum = DivideHistograms1D(hisPtInJetsSpectrum, hisEffPtInJetsExt[iCent]);
            else
              hisPtInJetsSpectrum = DivideHistograms1D(hisPtInJetsSpectrum, hisEffInJC);
            if(!hisPtInJetsSpectrum)
              return;
            printf("Correcting %s for efficiency\n", hisPtInPCSpectrum->GetName());
            if(iCorrEffInclusive)
              hisPtInPCSpectrum = DivideHistograms1D(hisPtInPCSpectrum, hisEffPtInJetsExt[iCent]);
            else
              hisPtInPCSpectrum = DivideHistograms1D(hisPtInPCSpectrum, hisEffInPC);
            if(!hisPtInPCSpectrum)
              return;
          }

          hisPtInBulkSpectrum[4] = hisPtInPCSpectrum;
          sPtInBulkName[4] = "PC";
          for(Int_t iUE = 0; iUE < iNUE; iUE++)
          {
            TH1D* hisPtBulkMinusRef = (TH1D*)hisPtInBulkSpectrum[iUE]->Clone();
            hisPtBulkMinusRef->Add(hisPtInBulkSpectrumRef, -1);
            hisPtBulkMinusRef->Divide(hisPtInBulkSpectrumRef);
            TH1D* hisPtBulkOverRef = DivideHistograms1D(hisPtInBulkSpectrum[iUE], hisPtInBulkSpectrumRef);
            hisPtBulkOverRef->Write(Form("fh1PtBulkOverRef%s_C%d-J%d-%s", V0Name[V0Type].Data(), iCent, iJet, sPtInBulkName[iUE].Data()));
            TH1D* hisPtBulkOverJC = DivideHistograms1D(hisPtInBulkSpectrum[iUE], hisPtInJetsSpectrum);
            hisPtBulkOverJC->Write(Form("fh1PtBulkOverJC%s_C%d-J%d-%s", V0Name[V0Type].Data(), iCent, iJet, sPtInBulkName[iUE].Data()));
            TGraphErrors* grPtBulkOverJC = MakeGraphErrors(hisPtBulkOverJC, sPtInBulkName[iUE].Data(), iMyColors[iUE], iMyMarkersFull[iUE]);
            mgrBulkOverJC->Add(grPtBulkOverJC);
            hisPtInJetsSpectrumUE[iUE] = (TH1D*)hisPtInJetsSpectrum->Clone();
            // subtract V0s in underlying event
            hisPtInJetsSpectrumUE[iUE]->Add(hisPtInBulkSpectrum[iUE], -1);
            TH1D* hisPtBulkOverSignal = DivideHistograms1D(hisPtInBulkSpectrum[iUE], hisPtInJetsSpectrumUE[iUE]);
            hisPtBulkOverSignal->Write(Form("fh1PtBulkOverSignal%s_C%d-J%d-%s", V0Name[V0Type].Data(), iCent, iJet, sPtInBulkName[iUE].Data()));
            TGraphErrors* grPtBulkOverSignal = MakeGraphErrors(hisPtBulkOverSignal, sPtInBulkName[iUE].Data(), iMyColors[iUE], iMyMarkersFull[iUE]);
            mgrBulkOverSignal->Add(grPtBulkOverSignal);
            hisPtBulkMinusRef->Multiply(hisPtBulkOverSignal);
            hisPtBulkMinusRef->Write(Form("fh1PtBulkMinusRef%s_C%d-J%d-%s", V0Name[V0Type].Data(), iCent, iJet, sPtInBulkName[iUE].Data()));
          }

          hisPtInBulkSpectrumSub = hisPtInBulkSpectrum[iUESub];
          if(iUEsubtraction != 0)
          {
            // subtract V0s in underlying event
            hisPtInJetsSpectrum->Add(hisPtInBulkSpectrumSub, -1);
          }

          // feed-down correction
          if(iCorrFD && (V0Type == 1 || V0Type == 2))
          {
            for(Int_t iV0 = iBinPtInJetsFirst - 1; iV0 <= iBinPtInJetsLast - 1; iV0++) // y (V0 pt bins)
            {
              printf("Correcting %s for feed-down, pT bin %d\n", hisPtInJetsSpectrum->GetName(), iV0 + 1);
              Double_t dFD;
              if(iCorrFDSel == 0) // PYTHIA 0-10 %
                dFD = GetMeanFeedDownLambdaPYTHIA(dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1]);
              if(iCorrFDSel == 1) // calculated inclusive
                dFD = hisFDInclJetBins[iCent]->GetBinContent(iV0 + 1);
              if(iCorrFDSel == 2) // LF 2010 0-10 %
                dFD = GetMeanFeedDownLambdaLF(dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1]);
              printf("Feed-down bin: %d, value: %g\n", iV0 + 1, dFD);
              hisPtInJetsSpectrum->SetBinContent(iV0 + 1, (1 - dFD)*hisPtInJetsSpectrum->GetBinContent(iV0 + 1));
              hisPtInJetsSpectrum->SetBinError(iV0 + 1, (1 - dFD)*hisPtInJetsSpectrum->GetBinError(iV0 + 1));
              for(Int_t iUE = 0; iUE < iNUE; iUE++)
              {
                hisPtInJetsSpectrumUE[iUE]->SetBinContent(iV0 + 1, (1 - dFD)*hisPtInJetsSpectrumUE[iUE]->GetBinContent(iV0 + 1));
                hisPtInJetsSpectrumUE[iUE]->SetBinError(iV0 + 1, (1 - dFD)*hisPtInJetsSpectrumUE[iUE]->GetBinError(iV0 + 1));
              }
            }
          }

          dirOutSpectra->cd();
          hisPtInJetsSpectrum->Write(Form("fh1PtInJets%s_C%d-J%d", V0Name[V0Type].Data(), iCent, iJet));
          for(Int_t iUE = 0; iUE < iNUE; iUE++)
            hisPtInJetsSpectrumUE[iUE]->Write(Form("fh1PtInJets%s_C%d-J%d-UE%d", V0Name[V0Type].Data(), iCent, iJet, iUE));
          hisPtInPCSpectrum->Write(Form("fh1PtInPC%s_C%d-J%d", V0Name[V0Type].Data(), iCent, iJet));
          hisPtInBulkSpectrumSub->Write(Form("fh1PtInBulk%s_C%d-J%d-Sub", V0Name[V0Type].Data(), iCent, iJet));
          hisPurityInJC->Write(Form("fh1PurityInJC%s_C%d-J%d", V0Name[V0Type].Data(), iCent, iJet));

          TGraphErrors* grPtInJCSpectrum;
          if(bOpenPtBins)
            grPtInJCSpectrum = MakeGraphErrors(hisPtInJCSpectrum, Form("#it{p}_{T}^{jet} > %1.f GeV/#it{c}", dBinsPtJet[iJet]), iMyColors[iJet], iMyMarkersFull[iJet]);
          else
            grPtInJCSpectrum = MakeGraphErrors(hisPtInJCSpectrum, Form("#it{p}_{T}^{jet}: %1.f-%1.f GeV/#it{c}", dBinsPtJet[iJet], dBinsPtJet[iJet + 1]), iMyColors[iJet], iMyMarkersFull[iJet]);
          mgrPtInJC->Add(grPtInJCSpectrum);
          delete hisPtInJCSpectrum;
          TGraphErrors* grPtInJetsSpectrum;
          if(bOpenPtBins)
            grPtInJetsSpectrum = MakeGraphErrors(hisPtInJetsSpectrum, Form("#it{p}_{T}^{jet} > %1.f GeV/#it{c}", dBinsPtJet[iJet]), iMyColors[iJet], iMyMarkersFull[iJet]);
          else
            grPtInJetsSpectrum = MakeGraphErrors(hisPtInJetsSpectrum, Form("#it{p}_{T}^{jet}: %1.f-%1.f GeV/#it{c}", dBinsPtJet[iJet], dBinsPtJet[iJet + 1]), iMyColors[iJet], iMyMarkersFull[iJet]);
          mgrPtInJets->Add(grPtInJetsSpectrum);
          delete hisPtInJetsSpectrum;
          TGraphErrors* grPurityInJC;
          if(bOpenPtBins)
            grPurityInJC = MakeGraphErrors(hisPurityInJC, Form("#it{p}_{T}^{jet} > %1.f GeV/#it{c}", dBinsPtJet[iJet]), iMyColors[iJet], iMyMarkersFull[iJet]);
          else
            grPurityInJC = MakeGraphErrors(hisPurityInJC, Form("#it{p}_{T}^{jet}: %1.f-%1.f GeV/#it{c}", dBinsPtJet[iJet], dBinsPtJet[iJet + 1]), iMyColors[iJet], iMyMarkersFull[iJet]);
          mgrPurityInJC->Add(grPurityInJC);
          delete hisPurityInJC;
          // bulk
          TGraphErrors* grPtInPCSpectrum;
          if(bOpenPtBins)
            grPtInPCSpectrum = MakeGraphErrors(hisPtInPCSpectrum, Form("PC #it{p}_{T}^{jet} > %1.f GeV/#it{c}", dBinsPtJet[iJet]), iMyColors[1], iMyMarkersEmpty[iJet]);
          else
            grPtInPCSpectrum = MakeGraphErrors(hisPtInPCSpectrum, Form("PC #it{p}_{T^{jet}}: %1.f-%1.f GeV/#it{c}", dBinsPtJet[iJet], dBinsPtJet[iJet + 1]), iMyColors[1], iMyMarkersEmpty[iJet]);
          mgrPtInPC->Add(grPtInPCSpectrum);
          TGraphErrors* grPtInPCSpectrumP;
          if(bOpenPtBins)
            grPtInPCSpectrumP = MakeGraphErrors(hisPtInPCSpectrum, Form("PC #it{p}_{T}^{jet} > %1.f GeV/#it{c}", dBinsPtJet[iJet]), iMyColors[1], iMyMarkersEmpty[iJet]);
          else
            grPtInPCSpectrumP = MakeGraphErrors(hisPtInPCSpectrum, Form("PC #it{p}_{T}: %1.f-%1.f GeV/#it{c}", dBinsPtJet[iJet], dBinsPtJet[iJet + 1]), iMyColors[1], iMyMarkersEmpty[iJet]);
          if(iJet == iJetSel)
            mgrPtInBulk->Add(grPtInPCSpectrumP);
//                      TH1D* hisPtInPCSpectrumComp = CompareHistograms1D(hisPtInPCSpectrum,hisPtInBulkSpectrumRef);
          TH1D* hisPtInPCSpectrumComp = DivideHistograms1D(hisPtInPCSpectrum, hisPtInBulkSpectrumRef);
          TGraphErrors* grPtInPCSpectrumComp = MakeGraphErrors(hisPtInPCSpectrumComp, Form("PC(%1.f GeV/#it{c})/%s", dBinsPtJet[iJet], sBulkRef.Data()), iMyColors[1], iMyMarkersEmpty[iCent]);
          if(iJet == iJetSel)
            mgrPtInBulkCompare->Add(grPtInPCSpectrumComp);
          delete hisPtInPCSpectrum;

          canPtInJetsAll->cd();
          if(bOpenPtBins)
            mgrPtInJetsAll->SetTitle(Form("%s sig+bg in JC, %s, #it{p}_{T}^{jet} > %.0f GeV/c", V0Name[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet]));
          else
            mgrPtInJetsAll->SetTitle(Form("%s sig+bg in JC, %s, #it{p}_{T}^{jet}: %.0f-%.0f GeV/c", V0Name[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet], dBinsPtJet[iJet + 1]));
          mgrPtInJetsAll->Draw("AP0");
          canPtInJetsAll->SetLogy();
          legend = canPtInJetsAll->BuildLegend(0.7, 0.6, 0.85, 0.85);
          SetLegend(legend, 0.02);
          canPtInJetsAll->SaveAs(Form("canPtInJCAll%s_C%d-J%d.%s", V0Name[V0Type].Data(), iCent, iJet, sImageSuf.Data()));
          delete canPtInJetsAll;

          canPtInJetsBg->cd();
          if(bOpenPtBins)
            mgrPtInJetsBg->SetTitle(Form("%s bg in JC, %s, #it{p}_{T}^{jet} > %.0f GeV/c", V0Name[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet]));
          else
            mgrPtInJetsBg->SetTitle(Form("%s bg in JC, %s, #it{p}_{T}^{jet}: %.0f-%.0f GeV/c", V0Name[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet], dBinsPtJet[iJet + 1]));
          mgrPtInJetsBg->Draw("AP0");
          canPtInJetsBg->SetLogy();
          legend = canPtInJetsBg->BuildLegend(0.7, 0.6, 0.85, 0.85);
          SetLegend(legend, 0.02);
          canPtInJetsBg->SaveAs(Form("canPtInJCBg%s_C%d-J%d.%s", V0Name[V0Type].Data(), iCent, iJet, sImageSuf.Data()));
          delete canPtInJetsBg;

          canPtInPCAll->cd();
          if(bOpenPtBins)
            mgrPtInPCAll->SetTitle(Form("%s sig+bg in PC, %s, #it{p}_{T}^{jet} > %.0f GeV/c", V0Name[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet]));
          else
            mgrPtInPCAll->SetTitle(Form("%s sig+bg in PC, %s, #it{p}_{T}^{jet}: %.0f-%.0f GeV/c", V0Name[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet], dBinsPtJet[iJet + 1]));
          mgrPtInPCAll->Draw("AP0");
          canPtInPCAll->SetLogy();
          legend = canPtInPCAll->BuildLegend(0.7, 0.6, 0.85, 0.85);
          SetLegend(legend, 0.02);
          canPtInPCAll->SaveAs(Form("canPtInPCAll%s_C%d-J%d.%s", V0Name[V0Type].Data(), iCent, iJet, sImageSuf.Data()));
          delete canPtInPCAll;

          canPtInPCBg->cd();
          if(bOpenPtBins)
            mgrPtInPCBg->SetTitle(Form("%s bg in PC, %s, #it{p}_{T}^{jet} > %.0f GeV/c", V0Name[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet]));
          else
            mgrPtInPCBg->SetTitle(Form("%s bg in PC, %s, #it{p}_{T}^{jet}: %.0f-%.0f GeV/c", V0Name[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet], dBinsPtJet[iJet + 1]));
          mgrPtInPCBg->Draw("AP0");
          canPtInPCBg->SetLogy();
          legend = canPtInPCBg->BuildLegend(0.7, 0.6, 0.85, 0.85);
          SetLegend(legend, 0.02);
          canPtInPCBg->SaveAs(Form("canPtInPCBg%s_C%d-J%d.%s", V0Name[V0Type].Data(), iCent, iJet, sImageSuf.Data()));
          delete canPtInPCBg;

          canBulkOverJC->cd();
          canBulkOverJC->SetLeftMargin(0.16);
          mgrBulkOverJC->SetTitle(Form("%s: bulk/JC ratio, c. %s, #it{p}_{T}^{jet} > %g GeV/c;#it{p}_{T}^{h} (GeV/#it{c});UE/JC", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet]));
          mgrBulkOverJC->SetMaximum(1.2);
          mgrBulkOverJC->SetMinimum(0);
          mgrBulkOverJC->Draw("AP0");
//                  canBulkOverJC->SetLogy();
          mgrBulkOverJC->GetYaxis()->SetTitleOffset(2);
          mgrBulkOverJC->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
          legend = canBulkOverJC->BuildLegend(0.6, 0.6, 0.85, 0.85);
          SetLegend(legend);
          //              labelSystem=labelCollision->DrawLatex(fPtInJetsXMin+(fPtInJetsXMax-fPtInJetsXMin)/5.,0.4,sLabelCollisionText.Data());
          canBulkOverJC->SaveAs(Form("canBulkOverJC%s_C%d_J%d.%s", V0Name[V0Type].Data(), iCent, iJet, sImageSuf.Data()));
          delete canBulkOverJC;
          delete mgrBulkOverJC;

          canBulkOverSignal->cd();
          canBulkOverSignal->SetLeftMargin(0.16);
          mgrBulkOverSignal->SetTitle(Form("%s: bulk/signal ratio, c. %s, #it{p}_{T}^{jet} > %g GeV/c;#it{p}_{T}^{h} (GeV/#it{c});UE/(JC - UE)", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet]));
          mgrBulkOverSignal->SetMaximum(50);
          mgrBulkOverSignal->SetMinimum(1e-2);
          mgrBulkOverSignal->Draw("AP0");
          canBulkOverSignal->SetLogy();
          mgrBulkOverSignal->GetYaxis()->SetTitleOffset(2);
          mgrBulkOverSignal->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
          legend = canBulkOverSignal->BuildLegend(0.6, 0.6, 0.85, 0.85);
          SetLegend(legend);
          //              labelSystem=labelCollision->DrawLatex(fPtInJetsXMin+(fPtInJetsXMax-fPtInJetsXMin)/5.,0.4,sLabelCollisionText.Data());
          canBulkOverSignal->SaveAs(Form("canBulkOverSignal%s_C%d_J%d.%s", V0Name[V0Type].Data(), iCent, iJet, sImageSuf.Data()));
          delete canBulkOverSignal;
          delete mgrBulkOverSignal;
        }

        if(iUESub != 4) // don't want PC to be considered for the final UE spectrum because UE should not depend on jet pt and PC do
        {
          if(iCorrFD && (V0Type == 1 || V0Type == 2))
          {
            for(Int_t iV0 = iBinPtInJetsFirst - 1; iV0 <= iBinPtInJetsLast - 1; iV0++) // y (V0 pt bins)
            {
              printf("Correcting %s for feed-down, pT bin %d\n", hisPtInBulkSpectrumSub->GetName(), iV0 + 1);
              Double_t dFD;
              if(iCorrFDSel == 2) // LF 2010 0-10 %
                dFD = GetMeanFeedDownLambdaLF(dBinsPtV0InJet[iV0], dBinsPtV0InJet[iV0 + 1]);
              else // calculated inclusive
                dFD = hisFDInclJetBins[iCent]->GetBinContent(iV0 + 1);
              hisPtInBulkSpectrumSub->SetBinContent(iV0 + 1, (1 - dFD)*hisPtInBulkSpectrumSub->GetBinContent(iV0 + 1));
              hisPtInBulkSpectrumSub->SetBinError(iV0 + 1, (1 - dFD)*hisPtInBulkSpectrumSub->GetBinError(iV0 + 1));
            }
          }
          hisPtInBulkSpectrumSub->Write(Form("fh1PtInBulk%s_C%d", V0Name[V0Type].Data(), iCent));
        }
        delete hisPtInNJSpectrum;
        delete hisPtInRCSpectrum;
        delete hisPtInMCCSpectrum;
        delete hisPtInOCSpectrum;

        canPtInJC->cd();
        canPtInJC->SetLeftMargin(0.16);
//                  mgrPtInJC->SetTitle(Form("%s: #it{p}_{T} spectrum, in JC, c. %s;#it{p}_{T}^{h} (GeV/#it{c});#frac{1}{#it{N}_{jet}} #frac{d#it{N}^{h}}{d#it{p}_{T}^{h}} (#it{c} GeV^{#minus1})",V0Symbol[V0Type].Data(),GetCentBinLabel(iCent).Data()));
        mgrPtInJC->SetTitle(Form("%s: #it{p}_{T} spectrum, in JC, c. %s;%s", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), sLabelsAxesSpectra.Data()));
        mgrPtInJC->SetMaximum(fPtInJetsYMax);
        mgrPtInJC->SetMinimum(fPtInJetsYMin);
        mgrPtInJC->Draw("AP0");
        mgrPtInJC->GetYaxis()->SetTitleOffset(2);
        mgrPtInJC->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
        canPtInJC->SetLogy();
        legend = canPtInJC->BuildLegend(0.5, 0.6, 0.85, 0.85);
        SetLegend(legend);
        labelSystem = labelCollision->DrawLatex(fPtInJetsXMin + (fPtInJetsXMax - fPtInJetsXMin) / 5., 10 * fPtInJetsYMin, sLabelCollisionText.Data());
        canPtInJC->SaveAs(Form("canPtInJCSpectrum%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canPtInJC;
        delete mgrPtInJC;

        canPtInJets->cd();
        canPtInJets->SetLeftMargin(0.16);
        mgrPtInJets->SetTitle(Form("%s: #it{p}_{T} spectrum, in jets, c. %s;%s", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), sLabelsAxesSpectra.Data()));
        mgrPtInJets->SetMaximum(fPtInJetsYMax);
        mgrPtInJets->SetMinimum(fPtInJetsYMin);
        mgrPtInJets->Draw("AP0");
        mgrPtInJets->GetYaxis()->SetTitleOffset(2);
        mgrPtInJets->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
        canPtInJets->SetLogy();
        legend = canPtInJets->BuildLegend(0.5, 0.6, 0.85, 0.85);
        SetLegend(legend);
        labelSystem = labelCollision->DrawLatex(fPtInJetsXMin + (fPtInJetsXMax - fPtInJetsXMin) / 5., 10 * fPtInJetsYMin, sLabelCollisionText.Data());
        canPtInJets->SaveAs(Form("canPtInJetsSpectrum%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canPtInJets;
        delete mgrPtInJets;

        canPurityInJC->cd();
        canPurityInJC->SetLeftMargin(0.16);
        mgrPurityInJC->SetTitle(Form("%s: signal purity in JC, c. %s;#it{p}_{T} (GeV/#it{c});purity", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrPurityInJC->SetMaximum(1.1);
        mgrPurityInJC->SetMinimum(0);
        mgrPurityInJC->Draw("AP0");
        mgrPurityInJC->GetYaxis()->SetTitleOffset(2);
        mgrPurityInJC->GetXaxis()->SetLimits(fPtAllXMin, fPtAllXMax);
        legend = canPurityInJC->BuildLegend(0.5, 0.15, 0.85, 0.35);
        SetLegend(legend);
        labelSystem = labelCollision->DrawLatex(fPtAllXMin + (fPtAllXMax - fPtAllXMin) / 5., 0.2, sLabelCollisionText.Data());
        canPurityInJC->SaveAs(Form("canPtInJCPurity%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canPurityInJC;
        delete mgrPurityInJC;

        canPtInPC->cd();
        canPtInPC->SetLeftMargin(0.16);
        mgrPtInPC->SetTitle(Form("%s: #it{p}_{T} spectrum, in PC, c. %s;%s", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), sLabelsAxesSpectra.Data()));
        mgrPtInPC->SetMaximum(fPtInJetsYMax);
        mgrPtInPC->SetMinimum(fPtInJetsYMin);
        mgrPtInPC->Draw("AP0");
        mgrPtInPC->GetYaxis()->SetTitleOffset(2);
        mgrPtInPC->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
        canPtInPC->SetLogy();
        legend = canPtInPC->BuildLegend(0.5, 0.6, 0.85, 0.85);
        SetLegend(legend);
        labelSystem = labelCollision->DrawLatex(fPtInJetsXMin + (fPtInJetsXMax - fPtInJetsXMin) / 5., 10 * fPtInJetsYMin, sLabelCollisionText.Data());
        canPtInPC->SaveAs(Form("canPtInPCSpectrum%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canPtInPC;
        delete mgrPtInPC;

        canPtInBulk->cd();
        canPtInBulk->SetLeftMargin(0.16);
//                  mgrPtInBulk->SetTitle(Form("%s: #it{p}_{T} spectrum, in bulk, c. %s;#it{p}_{T}^{h} (GeV/#it{c});#frac{1}{#it{N}_{jet}} #frac{d#it{N}^{h}}{d#it{p}_{T}^{h}} (#it{c} GeV^{#minus1})",V0Symbol[V0Type].Data(),GetCentBinLabel(iCent).Data()));
        mgrPtInBulk->SetTitle(Form("%s: #it{p}_{T} spectrum, in UE, c. %s;%s", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), sLabelsAxesSpectra.Data()));
        mgrPtInBulk->SetMaximum(fPtInJetsYMax);
        mgrPtInBulk->SetMinimum(fPtInJetsYMin);
        mgrPtInBulk->Draw("AP0");
        mgrPtInBulk->GetYaxis()->SetTitleOffset(2);
        mgrPtInBulk->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
        canPtInBulk->SetLogy();
        legend = canPtInBulk->BuildLegend(0.5, 0.6, 0.85, 0.85);
        SetLegend(legend);
        labelSystem = labelCollision->DrawLatex(fPtInJetsXMin + (fPtInJetsXMax - fPtInJetsXMin) / 5., 10 * fPtInJetsYMin, sLabelCollisionText.Data());
        canPtInBulk->SaveAs(Form("canPtInBulkSpectrum%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canPtInBulk;
        delete mgrPtInBulk;

        canPtInBulkCompare->cd();
        canPtInBulkCompare->SetLeftMargin(0.16);
        mgrPtInBulkCompare->SetTitle(Form("%s: #it{p}_{T} spectra comparison, in UE, c. %s;#it{p}_{T}^{h} (GeV/#it{c});ratio", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrPtInBulkCompare->SetMaximum(1.3);
        mgrPtInBulkCompare->SetMinimum(0.7);
        mgrPtInBulkCompare->Draw("AP0");
        mgrPtInBulkCompare->GetYaxis()->SetTitleOffset(2);
        mgrPtInBulkCompare->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
        legend = canPtInBulkCompare->BuildLegend(0.5, 0.15, 0.8, 0.3);
        SetLegend(legend);
        labelSystem = labelCollision->DrawLatex(fPtInJetsXMin + (fPtInJetsXMax - fPtInJetsXMin) / 5., 0.4, sLabelCollisionText.Data());
        canPtInBulkCompare->SaveAs(Form("canPtInBulkCompare%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canPtInBulkCompare;
        delete mgrPtInBulkCompare;

        canPtInRCAll->cd();
        mgrPtInRCAll->SetTitle(Form("%s sig+bg in RC, %s", V0Name[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrPtInRCAll->Draw("AP0");
        canPtInRCAll->SetLogy();
        legend = canPtInRCAll->BuildLegend(0.7, 0.6, 0.85, 0.85);
        SetLegend(legend, 0.02);
        canPtInRCAll->SaveAs(Form("canPtInRCAll%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canPtInRCAll;

        canPtInRCBg->cd();
        mgrPtInRCBg->SetTitle(Form("%s bg in RC, %s", V0Name[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrPtInRCBg->Draw("AP0");
        canPtInRCBg->SetLogy();
        legend = canPtInRCBg->BuildLegend(0.7, 0.6, 0.85, 0.85);
        SetLegend(legend, 0.02);
        canPtInRCBg->SaveAs(Form("canPtInRCBg%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canPtInRCBg;
      }
      canPtInRC->cd();
      canPtInRC->SetLeftMargin(0.16);
      mgrPtInRC->SetTitle(Form("%s: #it{p}_{T} spectrum, in RC;#it{p}_{T}^{h} (GeV/#it{c});%s", V0Symbol[V0Type].Data(), sLabelsAxesSpectra.Data()));
      mgrPtInRC->SetMaximum(fPtInJetsYMax);
      mgrPtInRC->SetMinimum(fPtInJetsYMin);
      mgrPtInRC->Draw("AP0");
      mgrPtInRC->GetYaxis()->SetTitleOffset(2);
      mgrPtInRC->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
      canPtInRC->SetLogy();
      legend = canPtInRC->BuildLegend(0.5, 0.6, 0.85, 0.85);
      SetLegend(legend);
      labelSystem = labelCollision->DrawLatex(fPtInJetsXMin + (fPtInJetsXMax - fPtInJetsXMin) / 5., 10 * fPtInJetsYMin, sLabelCollisionText.Data());
      canPtInRC->SaveAs(Form("canPtInRCSpectrum%s.%s", V0Name[V0Type].Data(), sImageSuf.Data()));
      delete canPtInRC;
      delete mgrPtInRC;
    }

    if(iCorrelations && V0Type != 2)
    {
      printf("Correlations V0-jet\n");
      TString sNameSpaCorrelSE = Form("fhnV0CorrelSE%s_%%d", V0Name[V0Type].Data()); // m_V0; pt_V0; eta_V0; pt_jet; delta-phi_V0-jet; (delta-eta_V0-jet)
      THnSparseD* spaCorrelSE; // same-event correlations
      TString sNameSpaCorrelME = Form("fhnV0CorrelME%s_%%d", V0Name[V0Type].Data());
      THnSparseD* spaCorrelME; // mixed-event correlations

      for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
      {
        printf("Cent %d: Start\n", iCent);
        if(!eventsCent[iCent])
          continue;

        spaCorrelSE = GetSparseD(listStd, Form(sNameSpaCorrelSE.Data(), iCent));
        if(!spaCorrelSE)
          continue;
        spaCorrelME = GetSparseD(listStd, Form(sNameSpaCorrelME.Data(), iCent));
        if(!spaCorrelME)
          continue;

        // restrict pt_jet, axis 3
        for(Int_t iJet = iJetMin; iJet <= iJetMax; iJet++)
        {
          Double_t dPtJetMin = dBinsPtJet[iJet];
          Double_t dPtJetMax = dBinsPtJet[iJet + 1];
          Int_t iPtJetBinFirst = spaCorrelSE->GetAxis(3)->FindBin(dPtJetMin + dEpsilon);
          Int_t iPtJetBinLast = spaCorrelSE->GetAxis(3)->FindBin(dPtJetMax - dEpsilon);
          if(bOpenPtBins)
            iPtJetBinLast = spaCorrelSE->GetAxis(3)->GetNbins() + 1;
          dPtJetMin = spaCorrelSE->GetAxis(3)->GetBinLowEdge(iPtJetBinFirst);
          dPtJetMax = spaCorrelSE->GetAxis(3)->GetBinUpEdge(iPtJetBinLast);
          printf("Jet pt range: %g-%g\n", dPtJetMin, dPtJetMax);
          spaCorrelSE->GetAxis(3)->SetRange(iPtJetBinFirst, iPtJetBinLast);
          spaCorrelME->GetAxis(3)->SetRange(iPtJetBinFirst, iPtJetBinLast);

          // histogram for pt spectrum
          TH1D* hisPtCorrelYield = new TH1D(Form("hisPtCorrelYield-C%d-J%d", iCent, iJet), Form("Yield of %s in jets %g-%g GeV/c, cent %s;GeV/c;N_{V0-jet}", V0Symbol[V0Type].Data(), dPtJetMin, dPtJetMax, GetCentBinLabel(iCent).Data()), iNBinsPtV0InJet, dBinsPtV0InJet);
          TH1D* hisPtCorrelSE = new TH1D(Form("hisPtCorrelSE-C%d-J%d", iCent, iJet), Form("pt spectrum of %s in correlation with jets in SE %g-%g GeV/c, cent %s;GeV/c;N_{V0-jet}", V0Symbol[V0Type].Data(), dPtJetMin, dPtJetMax, GetCentBinLabel(iCent).Data()), iNBinsPtV0InJet, dBinsPtV0InJet);
          TH1D* hisPtCorrelME = new TH1D(Form("hisPtCorrelME-C%d-J%d", iCent, iJet), Form("pt spectrum of %s in correlation with jets in ME %g-%g GeV/c, cent %s;GeV/c;N_{V0-jet}", V0Symbol[V0Type].Data(), dPtJetMin, dPtJetMax, GetCentBinLabel(iCent).Data()), iNBinsPtV0InJet, dBinsPtV0InJet);

          TCanvas* canDPhiME = new TCanvas(Form("canDPhiME-C%d-J%d", iCent, iJet), "", iCanWidth, iCanHeight);
          TMultiGraph* mgrDPhiME = new TMultiGraph();
          TGraphErrors* grDPhiME;
          TCanvas* canCF = new TCanvas(Form("canCF-C%d-J%d", iCent, iJet), "", iCanWidth, iCanHeight);
          TMultiGraph* mgrCF = new TMultiGraph();
          TGraphErrors* grCF;

          // restrict pt_V0, axis 1
          for(Int_t iV0 = iBinPtInJetsFirst - 1; iV0 <= iBinPtInJetsLast - 1; iV0++)
          {
            Int_t iPtV0BinFirst = spaCorrelSE->GetAxis(1)->FindBin(dBinsPtV0InJet[iV0] + dEpsilon);
            Int_t iPtV0BinLast = spaCorrelSE->GetAxis(1)->FindBin(dBinsPtV0InJet[iV0 + 1] - dEpsilon);
            Double_t dPtV0Min = spaCorrelSE->GetAxis(1)->GetBinLowEdge(iPtV0BinFirst);
            Double_t dPtV0Max = spaCorrelSE->GetAxis(1)->GetBinUpEdge(iPtV0BinLast);
            printf("V0 pt range: %g-%g\n", dPtV0Min, dPtV0Max);
            spaCorrelSE->GetAxis(1)->SetRange(iPtV0BinFirst, iPtV0BinLast);
            spaCorrelSE->GetAxis(2)->SetRange();
            spaCorrelSE->GetAxis(4)->SetRange();
            spaCorrelME->GetAxis(1)->SetRange(iPtV0BinFirst, iPtV0BinLast);
            spaCorrelME->GetAxis(2)->SetRange();
            spaCorrelME->GetAxis(4)->SetRange();

            TH1D* hisDPhiSE = (TH1D*)spaCorrelSE->Projection(4, "e");
            dirOutSpectra->cd();
            hisDPhiSE->Write(Form("hisDPhiSE%s-C%d-J%d-V%d-Raw", V0Name[V0Type].Data(), iCent, iJet, iV0)); // final delta-phi distribution for given pt V0
            hisDPhiSE->Reset();
            TH1D* hisDPhiSETmp;
            hisDPhiSETmp = (TH1D*)hisDPhiSE->Clone(Form("hisDPhiSETmp-C%d-J%d-V%d", iCent, iJet, iV0)); // temporary delta-phi distribution for loop over eta

            TH1D* hisDPhiME = (TH1D*)spaCorrelME->Projection(4, "e");
            dirOutSpectra->cd();
            hisDPhiME->Write(Form("hisDPhiME%s-C%d-J%d-V%d-Raw", V0Name[V0Type].Data(), iCent, iJet, iV0)); // final delta-phi distribution for given pt V0
            hisDPhiME->Reset();
            TH1D* hisDPhiMETmp = (TH1D*)hisDPhiME->Clone(Form("hisDPhiMETmp-C%d-J%d-V%d", iCent, iJet, iV0)); // temporary delta-phi distribution for loop over eta

            // restrict eta_V0, axis 2, one efficiency and one purity for all delta-phi_V0-jet
            for(Int_t iEta = 0; iEta < iNBinsEtaInJet; iEta++)
            {
              Int_t iEtaV0BinFirst = spaCorrelSE->GetAxis(2)->FindBin(dBinsEtaInJet[iEta] + dEpsilon);
              Int_t iEtaV0BinLast = spaCorrelSE->GetAxis(2)->FindBin(dBinsEtaInJet[iEta + 1] - dEpsilon);
              Double_t dEtaV0Min = spaCorrelSE->GetAxis(2)->GetBinLowEdge(iEtaV0BinFirst);
              Double_t dEtaV0Max = spaCorrelSE->GetAxis(2)->GetBinUpEdge(iEtaV0BinLast);
              printf("V0 eta range %g - %g, Eta bins: %d-%d\n", dEtaV0Min, dEtaV0Max, iEtaV0BinFirst, iEtaV0BinLast);
              spaCorrelSE->GetAxis(2)->SetRange(iEtaV0BinFirst, iEtaV0BinLast);
              spaCorrelME->GetAxis(2)->SetRange(iEtaV0BinFirst, iEtaV0BinLast);
              spaCorrelSE->GetAxis(4)->SetRange(); // reset delta-phi_V0-jet range
              spaCorrelME->GetAxis(4)->SetRange();
              spaCorrelSE->GetAxis(0)->SetRange();
              spaCorrelME->GetAxis(0)->SetRange();

              hisDPhiSETmp->Reset(); // reset delta-phi distribution for each eta bin
              hisDPhiMETmp->Reset(); // reset delta-phi distribution for each eta bin

              BinCounterObject* binCount;
              TH1D* hisMass;
              Bool_t bStatusBC;
              Double_t dIntegral;
              Double_t dSignalErr;
              Double_t dSignal;
              Double_t dPuritySEErr;
              Double_t dPuritySE;
              Double_t dPurityMEErr;
              Double_t dPurityME;
              TString sTitle;

              // signal extraction, axis 0, integrate over delta-phi_V0-jet
              // SE
              hisMass = (TH1D*)spaCorrelSE->Projection(0, "e");
              if(!hisMass)
                return;
              if(iNRebin > 1)
                hisMass->Rebin(iNRebin);
              dIntegral = hisMass->Integral();
              hisMass->SetName(Form("CorrelSE%s-C%d-J%d-V%d-E%d", V0Name[V0Type].Data(), iCent, iJet, iV0, iEta));
              printf("%s: entries: %g\n", hisMass->GetName(), dIntegral);

              if(dIntegral < 1)
              {
                printf("No entries, skipping\n");
              }
              else
              {
                printf("%s: BC Start\n", hisMass->GetName());
                bStatusBC = kTRUE;
                binCount = new BinCounterObject(hisMass, fMassV0[V0Type], fSigmaV0[V0Type]);
//                  binCount->SetVerbose(bVerboseBC);
                binCount->SetDegreePolMax(iDegreePolSideBands);
                bStatusBC &= binCount->SetXLimits(dBgOutMin, dBgOutMax);
                if(bFixedBC)
                  bStatusBC &= binCount->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax); // normal settings
                else
                  bStatusBC &= binCount->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
                binCount->SetFitOptionSideBands(sOptionFitSB.Data());
                sTitle = hisMass->GetName();
                bStatusBC &= binCount->EstimateParameters();
                bStatusBC &= binCount->Run();
                binCount->Plot("canBinCounter2", Form("%s, step %d", sTitle.Data(), 2));
                printf("%s: BC End\n", hisMass->GetName());
                if(!bStatusBC)
                {
                  printf("Something wrong with bin counting\n");
                  delete binCount;
                  continue;
                }
                dSignalErr = 0;
                dSignal = 0;
                dSignal = binCount->GetSignalAndError(&dSignalErr);
                dPuritySE = binCount->GetPurityAndError(&dPuritySEErr);
                delete binCount;
              }

              // ME
              hisMass = (TH1D*)spaCorrelME->Projection(0, "e");
              if(!hisMass)
                return;
              if(iNRebin > 1)
                hisMass->Rebin(iNRebin);
              dIntegral = hisMass->Integral();
              hisMass->SetName(Form("CorrelME%s-C%d-J%d-V%d-E%d", V0Name[V0Type].Data(), iCent, iJet, iV0, iEta));
              printf("%s: entries: %g\n", hisMass->GetName(), dIntegral);

              if(dIntegral < 1)
              {
                printf("No entries, skipping\n");
              }
              else
              {
                printf("%s: BC Start\n", hisMass->GetName());
                bStatusBC = kTRUE;
                binCount = new BinCounterObject(hisMass, fMassV0[V0Type], fSigmaV0[V0Type]);
//                  binCount->SetVerbose(bVerboseBC);
                binCount->SetDegreePolMax(iDegreePolSideBands);
                bStatusBC &= binCount->SetXLimits(dBgOutMin, dBgOutMax);
                if(bFixedBC)
                  bStatusBC &= binCount->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax); // normal settings
                else
                  bStatusBC &= binCount->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
                binCount->SetFitOptionSideBands(sOptionFitSB.Data());
                sTitle = hisMass->GetName();
                bStatusBC &= binCount->EstimateParameters();
                bStatusBC &= binCount->Run();
                binCount->Plot("canBinCounter2", Form("%s, step %d", sTitle.Data(), 2));
                printf("%s: BC End\n", hisMass->GetName());
                if(!bStatusBC)
                {
                  printf("Something wrong with bin counting\n");
                  delete binCount;
                  continue;
                }
                dSignalErr = 0;
                dSignal = 0;
                dSignal = binCount->GetSignalAndError(&dSignalErr);
                dPurityME = binCount->GetPurityAndError(&dPurityMEErr);
                delete binCount;
              }

              // project on delta-phi_V0-jet, axis 4
              Int_t iMassBinFirst = spaCorrelSE->GetAxis(0)->FindBin(dSigMin + dEpsilon);
              Int_t iMassBinLast = spaCorrelSE->GetAxis(0)->FindBin(dSigMax - dEpsilon);
              Double_t dMassMin = spaCorrelSE->GetAxis(0)->GetBinLowEdge(iMassBinFirst);
              Double_t dMassMax = spaCorrelSE->GetAxis(0)->GetBinUpEdge(iMassBinLast);
              printf("V0 mass range %g - %g, mass bins: %d-%d\n", dMassMin, dMassMax, iMassBinFirst, iMassBinLast);
              spaCorrelSE->GetAxis(0)->SetRange(iMassBinFirst, iMassBinLast);
              spaCorrelME->GetAxis(0)->SetRange(iMassBinFirst, iMassBinLast);
              TH1D* hisDPhiSETmp = (TH1D*)spaCorrelSE->Projection(4, "e");
              TH1D* hisDPhiMETmp = (TH1D*)spaCorrelME->Projection(4, "e");

              // scale with purity for given eta_V0
              hisDPhiSETmp = MultiplyHistogram(hisDPhiSETmp, dPuritySE, dPuritySEErr);
              hisDPhiMETmp = MultiplyHistogram(hisDPhiMETmp, dPurityME, dPurityMEErr);

              // correct for efficiency for given eta_V0
              if(iCorrEff)
              {
                Double_t dEfficiency = hisEffPtEtaInclBase[iCent]->GetBinContent(iV0 + 1, iEta + 1);
                Double_t dEfficiencyErr = hisEffPtEtaInclBase[iCent]->GetBinError(iV0 + 1, iEta + 1);
                hisDPhiSETmp = DivideHistogram(hisDPhiSETmp, dEfficiency, dEfficiencyErr);
                hisDPhiMETmp = DivideHistogram(hisDPhiMETmp, dEfficiency, dEfficiencyErr);
              }
              hisDPhiSE->Add(hisDPhiSETmp);
              hisDPhiME->Add(hisDPhiMETmp);
            }

            // integrated over eta
            dirOutSpectra->cd();
            hisDPhiSE->Write(Form("hisDPhiSE%s-C%d-J%d-V%d-Signal", V0Name[V0Type].Data(), iCent, iJet, iV0));
            hisDPhiME->Write(Form("hisDPhiME%s-C%d-J%d-V%d-Signal", V0Name[V0Type].Data(), iCent, iJet, iV0));

            // integrate over delta-phi and store the yield in the pt spectrum
            Double_t dNPairs = hisDPhiSE->Integral();
            Double_t dSignal = dNPairs;
            Double_t dSignalErr = GetHistError(hisDPhiSE);
            hisPtCorrelSE->SetBinContent(iV0 + 1, dSignal);
            hisPtCorrelSE->SetBinError(iV0 + 1, dSignalErr);
            TH1D* hisDPhiSENorm = DivideHistogram(hisDPhiSE, dSignal, dSignalErr);
            delete hisDPhiSE;
            if(!hisDPhiSENorm)
              continue;
            delete hisDPhiSETmp;
            dSignal = hisDPhiME->Integral();
            dSignalErr = GetHistError(hisDPhiME);
            hisPtCorrelME->SetBinContent(iV0 + 1, dSignal);
            hisPtCorrelME->SetBinError(iV0 + 1, dSignalErr);
            TH1D* hisDPhiMENorm = DivideHistogram(hisDPhiME, dSignal, dSignalErr);
            delete hisDPhiME;
            if(!hisDPhiMENorm)
              continue;
            delete hisDPhiMETmp;

            // get correlation function C
            TH1D* hisCorrelFun = DivideHistograms1D(hisDPhiSENorm, hisDPhiMENorm, Form("Correlation function C%d-J%d-V%d", iCent, iJet, iV0));
            // store integral of correlation function for normalisation
            Double_t dIntegralCorrelFun = hisCorrelFun->Integral();
            hisCorrelFun->Write(Form("hisCorrelFun_C%d-J%d-V%d", iCent, iJet, iV0));
            grDPhiME = MakeGraphErrors(hisDPhiMENorm, Form("#it{p}_{T}^{V^{0}}: %g-%g GeV/#it{c}", dPtV0Min, dPtV0Max), iMyColors[iV0 % iNMyColors], iMyMarkersEmpty[iV0 % iNMyMarkersEmpty]);
            mgrDPhiME->Add(grDPhiME);
            grCF = MakeGraphErrors(hisCorrelFun, Form("#it{p}_{T}^{V^{0}}: %g-%g GeV/#it{c}", dPtV0Min, dPtV0Max), iMyColors[iV0 % iNMyColors], iMyMarkersEmpty[iV0 % iNMyMarkersEmpty]);
            mgrCF->Add(grCF);

            // estimate ZYAM
            Double_t dZyam = hisCorrelFun->GetBinContent(hisCorrelFun->GetMinimumBin());

            // subtract ZYAM from C
            for(Int_t iBinPhi = 1; iBinPhi <= hisCorrelFun->GetNbinsX(); iBinPhi++)
              hisCorrelFun->AddBinContent(iBinPhi, -dZyam);
            hisCorrelFun->Write(Form("hisCorrelFun_C%d-J%d-V%d_Subtracted", iCent, iJet, iV0));

            // normalise
            hisCorrelFun->Scale(dNPairs / dIntegralCorrelFun);

            // integrate over delat-phi window
            Double_t dDPhiMin = -0.1;//-TMath::Pi()/2.; // minimum delta-phi window for integrating around the peak
            Double_t dDPhiMax = 0.1;//TMath::Pi()*3./2.; // maximum of delta-phi window for integrating around the peak
            Int_t iBinDPhiMin = hisCorrelFun->GetXaxis()->FindBin(dDPhiMin + dEpsilon);
            Int_t iBinDPhiMax = hisCorrelFun->GetXaxis()->FindBin(dDPhiMax - dEpsilon);
            printf("Integrating over dphi bins: %d-%d/%d\n", iBinDPhiMin, iBinDPhiMax, hisCorrelFun->GetNbinsX());
            Double_t dYieldInPeakErr;
            Double_t dYieldInPeak = hisCorrelFun->IntegralAndError(iBinDPhiMin, iBinDPhiMax, dYieldInPeakErr);

            // store yield in a pt spectrum
            hisPtCorrelYield->SetBinContent(iV0 + 1, dYieldInPeak);
            hisPtCorrelYield->SetBinError(iV0 + 1, dYieldInPeakErr);

          }
          dirOutSpectra->cd();
          hisPtCorrelSE->Write();
          hisPtCorrelME->Write();
          hisPtCorrelYield->Write();

          canDPhiME->cd();
          mgrDPhiME->SetTitle(Form("Azimuthal correlations jet-%s in mixed events;#Delta#it{#phi};arb. u.", V0Symbol[V0Type].Data()));
          mgrDPhiME->Draw("AP0");
          legend = canDPhiME->BuildLegend(0.5, 0.6, 0.85, 0.85);
          SetLegend(legend);
          canDPhiME->SaveAs(Form("canDPhiME%s-C%d-J%d.%s", V0Name[V0Type].Data(), iCent, iJet, sImageSuf.Data()));
          delete canDPhiME;
          delete mgrDPhiME;

          canCF->cd();
          mgrCF->SetTitle(Form("Correlation function for jet-%s pairs;#Delta#it{#phi};#it{C}", V0Symbol[V0Type].Data()));
          mgrCF->Draw("AP0");
          legend = canCF->BuildLegend(0.5, 0.6, 0.85, 0.85);
          SetLegend(legend);
          canCF->SaveAs(Form("canCF%s-C%d-J%d.%s", V0Name[V0Type].Data(), iCent, iJet, sImageSuf.Data()));
          delete canCF;
          delete mgrCF;
        }
      }
    }
  }

  if(iMC)
  {
    TH1D* hisMass = 0;

    // eff pT inclusive
    TCanvas* canEffPtIncl = new TCanvas("canEffPtIncl", "", iCanWidth, iCanHeight);
    TMultiGraph* mgrEffPtIncl = new TMultiGraph();
    TGraphErrors* grEffPtIncl;
    TString hisPtInclMCGenName = Form("fh1V0%sPtMCGen_%%d", V0Name[V0Type].Data());
    TString hisPtInclMCRec2DName = Form("fh2V0%sPtMassMCRec_%%d", V0Name[V0Type].Data());
    TH1D* hisPtInclMCGen;
    TH2D* hisPtInclMCRec2D;
    TH1D* hisPtInclMCRec1D;
    TH1D* hisEffPtIncl;

    // eff pT in jets
    TCanvas* canEffPtInJet[iNCentBins];
    TMultiGraph* mgrEffPtInJet[iNCentBins];
    TGraphErrors* grEffPtInJet;
    TString hisPtInJetsMCGen2DName = Form("fh2V0%sInJetPtMCGen_%%d", V0Name[V0Type].Data());
//      TString spaPtInJetsMCRec3DName = Form("fh3V0%sInJetPtMassMCRec_%%d",V0Name[V0Type].Data());
    TString spaPtInJetsMCRec3DName = Form("fh4V0%sInJetEtaPtMassMCRec_%%d", V0Name[V0Type].Data());
    TH2D* hisPtInJetsMCGen2D;
    THnSparseD* spaPtInJetsMCRec3D;
    TH1D* hisPtInJetsMCGen1D;
    TH1D* hisPtInJetsMCRec1D;
    TH1D* hisEffPtInJet;
    Double_t dArrayJetNormEff[iNBinsPtJet];

    // eff pT in jets/inclusive
    TCanvas* canEffPtRatio[iNCentBins];
    TMultiGraph* mgrEffPtRatio[iNCentBins];
    TGraphErrors* grEffPtRatio;
    TH1D* hisPtInclMCGenJ;
    TH1D* hisPtInclMCRec1DJ;
    TH1D* hisEffPtInclJ;
    TH1D* hisEffPtRatio;

    // inclusive eta-pT
    Double_t dBinsEta[] = { -fRangeEtaV0Max, -0.75, -0.7, -0.65, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, fRangeEtaV0Max};
    const Int_t iNBinsEta = sizeof(dBinsEta) / sizeof(dBinsEta[0]) - 1;
    Double_t dBinsPtELow[] =  {0.6, 1.,  2.,  3, 4, 6, 8};
    Double_t dBinsPtEHigh[] = {0.7, 1.1, 2.2, 4, 6, 8, 12};
    const Int_t iNBinsPtE = sizeof(dBinsPtELow) / sizeof(dBinsPtELow[0]);
    Double_t dBinsEtaTot[] = { -fRangeEtaV0Max, -0.75, 0.0, 0.75, fRangeEtaV0Max};
    const Int_t iNBinsEtaTot = sizeof(dBinsEtaTot) / sizeof(dBinsEtaTot[0]) - 1;
    // delta eta slices for eff pT-eta-deltaEta
    Double_t dBinsDeltaEtaInJet[] = { -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5};
    const Int_t iNBinsDeltaEtaInJet = sizeof(dBinsDeltaEtaInJet) / sizeof(dBinsDeltaEtaInJet[0]) - 1;

    // eff eta-pT inclusive
    TCanvas* canEffEtaPtIncl[iNCentBins];
    TMultiGraph* mgrEffEtaPtIncl[iNCentBins];
    TGraphErrors* grEffEtaPtIncl;
    TString hisEtaPtInclMCGen2DName = Form("fh2V0%sEtaPtMCGen_%%d", V0Name[V0Type].Data());
    TString spaEtaPtInclMCRec3DName = Form("fh3V0%sEtaPtMassMCRec_%%d", V0Name[V0Type].Data());
    TH2D* hisEtaPtInclMCGen2D;
    THnSparseD* spaEtaPtInclMCRec3D;
    TH1D* hisEtaPtInclMCGen1D;
    TH1D* hisEtaPtInclMCRec1D;
    TH1D* hisEffEtaPtIncl;
    TH2D* hisEffEtaPtIncl2D;

    // eff eta-pT injets
    TString spaEtaPtInJetsMCGen3DName = Form("fh3V0%sInJetEtaPtMCGen_%%d", V0Name[V0Type].Data());
    THnSparseD* spaEtaPtInJetsMCGen3D;

    // eff eta-pT daughters
    // inclusive
    TString spaEtaPtDInclMCRecName = Form("fhnV0%sInclDaughterEtaPtPtMCRec_%%d", V0Name[V0Type].Data());
    THnSparseD* spaEtaPtDInclMCRec;
    // in jets
    TString spaEtaPtDInJetsMCRecName = Form("fhnV0%sInJetsDaughterEtaPtPtMCRec_%%d", V0Name[V0Type].Data());
    THnSparseD* spaEtaPtDInJetsMCRec;

    // eff pT-eta inclusive
    TCanvas* canEffPtEtaIncl;
    TMultiGraph* mgrEffPtEtaIncl;
    TCanvas* canEffPtEtaInclPos;
    TMultiGraph* mgrEffPtEtaInclPos;
    TGraphErrors* grEffPtEtaIncl;
    TString hisPtEtaInclMCGen2DName = Form("fh2V0%sEtaPtMCGen_%%d", V0Name[V0Type].Data());
    TString spaPtEtaInclMCRec3DName = Form("fh3V0%sEtaPtMassMCRec_%%d", V0Name[V0Type].Data());
    TH2D* hisPtEtaInclMCGen2D;
    THnSparseD* spaPtEtaInclMCRec3D;
    TH1D* hisPtEtaInclMCGen1D;
    TH1D* hisPtEtaInclMCRec1D;
    TH1D* hisEffPtEtaIncl;

    // eff pT-eta in jets
    TCanvas* canEffPtEtaInJet;
    TMultiGraph* mgrEffPtEtaInJet;
    TCanvas* canEffPtEtaInJetPos;
    TMultiGraph* mgrEffPtEtaInJetPos;
    TGraphErrors* grEffPtEtaInJet;
    TString spaPtEtaInJetMCGen3DName = Form("fh3V0%sInJetEtaPtMCGen_%%d", V0Name[V0Type].Data());
    TString spaPtEtaInJetMCRec4DName = Form("fh4V0%sInJetEtaPtMassMCRec_%%d", V0Name[V0Type].Data());
    THnSparseD* spaPtEtaInJetMCGen3D;
    THnSparseD* spaPtEtaInJetMCRec4D;
    TH1D* hisPtEtaInJetMCGen1D;
    TH1D* hisPtEtaInJetMCRec1D;
    TH1D* hisEffPtEtaInJet;

    // eff pT-eta ratio in jets/inclusive
    TCanvas* canEffPtEtaRatio;
    TMultiGraph* mgrEffPtEtaRatio;
    TCanvas* canEffPtEtaRatioPos;
    TMultiGraph* mgrEffPtEtaRatioPos;
    TGraphErrors* grEffPtEtaRatio;
    TH1D* hisEffPtEtaRatio;
    // gen pT-eta ratio in jets/inclusive
    TCanvas* canGenPtEtaRatio;
    TMultiGraph* mgrGenPtEtaRatio;
    TCanvas* canGenPtEtaRatioPos;
    TMultiGraph* mgrGenPtEtaRatioPos;
    TGraphErrors* grGenPtEtaRatio;
    TH1D* hisGenPtEtaRatio;

    // eff eta/pT in jets/inclusive

    dirOutSpectra->cd();
    for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
    {
      if(!eventsCent[iCent])
        continue;

      if(iEffPtInJets)
      {
        if(!hisJetPt[iCent])
          return;
        Int_t iNBinsPtJetNormEff = hisJetPt[iCent]->GetXaxis()->GetNbins();
        for(Int_t i = 0; i < iNBinsPtJet; i++)
        {
          Int_t iBinFirst = hisJetPt[iCent]->GetXaxis()->FindBin(dBinsPtJet[i] + dEpsilon);
          Int_t iBinLast;
          if(bOpenPtBins)
            iBinLast = iNBinsPtJetNormEff + 1;
          else
            iBinLast = hisJetPt[iCent]->GetXaxis()->FindBin(dBinsPtJet[i + 1] - dEpsilon);
          dArrayJetNormEff[i] = hisJetPt[iCent]->Integral(iBinFirst, iBinLast);
        }
        printf("Statistics: Jets MC (Cent %d):", iCent);
        for(Int_t iJet = iJetMin; iJet <= iJetMax; iJet++)
          printf(" %g", dArrayJetNormEff[iJet]);
        printf("\n");
      }
      if(iEffPtInclusive)
      {
        hisPtInclMCGen = GetHistogram1D(listMC, Form(hisPtInclMCGenName.Data(), iCent));
        if(!hisPtInclMCGen)
          return;
        hisPtInclMCRec2D = GetHistogram2D(listMC, Form(hisPtInclMCRec2DName.Data(), iCent));
        if(!hisPtInclMCRec2D)
          return;
        hisPtInclMCRec1D = new TH1D(Form("hisCandRecAll%d", iCent), Form("MC %s signal candidates all, %s;#it{p}_{T} (GeV/c)", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()), iNBinsPtV0AllLF, dBinsPtV0AllLF);
        hisPtInclMCRec1DJ = new TH1D(Form("hisCandRecJ%d", iCent), Form("MC %s signal candidates J, %s;#it{p}_{T} (GeV/c)", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()), iNBinsPtV0InJet, dBinsPtV0InJet);

        hisPtInclMCGenJ = (TH1D*)hisPtInclMCGen->Clone(Form("%sJ", hisPtInclMCGen->GetName()));
        hisPtInclMCGenJ = (TH1D*)hisPtInclMCGenJ->Rebin(iNBinsPtV0InJet, "", dBinsPtV0InJet);
        hisPtInclMCGen = (TH1D*)hisPtInclMCGen->Rebin(iNBinsPtV0AllLF, "", dBinsPtV0AllLF);
        if(bZeroGenErr)
        {
          for(Int_t iB = 1; iB <= hisPtInclMCGen->GetNbinsX(); iB++)
            hisPtInclMCGen->SetBinError(iB, 0);
          for(Int_t iB = 1; iB <= hisPtInclMCGenJ->GetNbinsX(); iB++)
            hisPtInclMCGenJ->SetBinError(iB, 0);
        }
        CheckHistogram(hisPtInclMCGen);
//              dirOutSpectra->cd();
//              hisPtInclMCGen->Write(Form("fh1PtGenInclusive%s_C%d",V0Name[V0Type].Data(),iCent));
//              continue;

        // eff inclusive
        for(Int_t iPt = iBinPtInclFirst; iPt <= iBinPtInclLast; iPt++)
        {
          Int_t iBinPtFirst = hisPtInclMCRec2D->GetXaxis()->FindBin(dBinsPtV0AllLF[iPt - 1] + dEpsilon);
          Int_t iBinPtLast = hisPtInclMCRec2D->GetXaxis()->FindBin(dBinsPtV0AllLF[iPt] - dEpsilon);
          hisMass = (TH1D*)hisPtInclMCRec2D->ProjectionY(Form("MCPtInclusive%s_C%d-Pt%d", V0Name[V0Type].Data(), iCent, iPt), iBinPtFirst, iBinPtLast, "e");
          if(!hisMass)
            return;
          if(iNRebin > 1)
            hisMass->Rebin(iNRebin);
          printf("iEffPtInclusive: %s: %g - %g : %g, under: %g, over: %g\n", V0Name[V0Type].Data(), dBinsPtV0AllLF[iPt - 1], dBinsPtV0AllLF[iPt], hisMass->Integral(), hisMass->GetBinContent(0), hisMass->GetBinContent(hisMass->GetNbinsX()));

          Double_t dSigmaOld = MassPeakSigmaOld((dBinsPtV0AllLF[iPt - 1] + dBinsPtV0AllLF[iPt]) / 2, V0Type);
          Double_t dSigMinOld = fMassV0[V0Type] - 3 * dSigmaOld;
          Double_t dSigMaxOld = fMassV0[V0Type] + 3 * dSigmaOld;

          Bool_t bStatusBC = kTRUE;
          BinCounterObject* binCount = new BinCounterObject(hisMass, fMassV0[V0Type], fSigmaV0[V0Type]);
          binCount->SetVerbose(bVerboseBC);
          binCount->SetDegreePolMax(0);
          bStatusBC &= binCount->SetXLimits(dBgOutMin, dBgOutMax);
          if(bFixedBC)
            bStatusBC &= binCount->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax);
//                        bStatusBC &= binCount->FixRegions(dSigMinOld,dSigMaxOld,dBgInMin,dBgInMax,dBgOutMin,dBgOutMax);
          else
            bStatusBC &= binCount->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
//                  binCount->SetConsiderBackground(0);
          binCount->SetBackgroundSubtraction(0);
//                  binCount->Plot("canBinCounter0");
          bStatusBC &= binCount->EstimateParameters();
//                  binCount->Plot("canBinCounter1");
          bStatusBC &= binCount->Run();
          binCount->Plot("canBinCounter2");
          if(!bStatusBC)
          {
            printf("Something wrong with bin counting\n");
            delete binCount;
            continue;
          }
          Double_t dSignalErr = 0;
          Double_t dSignal = binCount->GetSignalAndError(&dSignalErr);
          delete binCount;

          hisPtInclMCRec1D->SetBinContent(iPt, dSignal);
          hisPtInclMCRec1D->SetBinError(iPt, dSignalErr);
        }
        // eff inclusive, rebin for jets and bulk

        for(Int_t iPt = iBinPtInJetsFirst; iPt <= iBinPtInJetsLast; iPt++)
        {
          Int_t iBinPtFirst = hisPtInclMCRec2D->GetXaxis()->FindBin(dBinsPtV0InJet[iPt - 1] + dEpsilon);
          Int_t iBinPtLast = hisPtInclMCRec2D->GetXaxis()->FindBin(dBinsPtV0InJet[iPt] - dEpsilon);
          hisMass = (TH1D*)hisPtInclMCRec2D->ProjectionY(Form("MCPtInclusiveRebin%s_C%d-Pt%d", V0Name[V0Type].Data(), iCent, iPt), iBinPtFirst, iBinPtLast, "e");
          if(!hisMass)
            return;
          if(iNRebin > 1)
            hisMass->Rebin(iNRebin);

          Double_t dSigmaOld = MassPeakSigmaOld((dBinsPtV0InJet[iPt - 1] + dBinsPtV0InJet[iPt]) / 2, V0Type);
          Double_t dSigMinOld = fMassV0[V0Type] - 3 * dSigmaOld;
          Double_t dSigMaxOld = fMassV0[V0Type] + 3 * dSigmaOld;

          Bool_t bStatusBC = kTRUE;
          BinCounterObject* binCount = new BinCounterObject(hisMass, fMassV0[V0Type], fSigmaV0[V0Type]);
          binCount->SetVerbose(bVerboseBC);
          binCount->SetDegreePolMax(0);
          bStatusBC &= binCount->SetXLimits(dBgOutMin, dBgOutMax);
          if(bFixedBC)
            bStatusBC &= binCount->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax);
//                        bStatusBC &= binCount->FixRegions(dSigMinOld,dSigMaxOld,dBgInMin,dBgInMax,dBgOutMin,dBgOutMax);
          else
            bStatusBC &= binCount->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
//                  binCount->SetConsiderBackground(0);
          binCount->SetBackgroundSubtraction(0);
//                  binCount->Plot("canBinCounter0");
          bStatusBC &= binCount->EstimateParameters();
//                  binCount->Plot("canBinCounter1");
          bStatusBC &= binCount->Run();
          binCount->Plot("canBinCounter2");
          if(!bStatusBC)
          {
            printf("Something wrong with bin counting\n");
            delete binCount;
            continue;
          }
          Double_t dSignalErr = 0;
          Double_t dSignal = binCount->GetSignalAndError(&dSignalErr);
          delete binCount;

          hisPtInclMCRec1DJ->SetBinContent(iPt, dSignal);
          hisPtInclMCRec1DJ->SetBinError(iPt, dSignalErr);
        }

        hisEffPtIncl = DivideHistograms1D(hisPtInclMCRec1D, hisPtInclMCGen, Form("hisEffPtIncl-%d", iCent));
        if(!hisEffPtIncl)
          return;
        hisEffPtInclJ = DivideHistograms1D(hisPtInclMCRec1DJ, hisPtInclMCGenJ, Form("hisEffJ-%d", iCent));
        if(!hisEffPtInclJ)
          return;
        hisPtInclMCRec1D->Write(Form("fh1RecPtIncl%s_%d", V0Name[V0Type].Data(), iCent));
        hisPtInclMCGen->Write(Form("fh1GenPtIncl%s_%d", V0Name[V0Type].Data(), iCent));
        hisEffPtIncl->Write(Form(sNameHisEffPtIncl.Data(), V0Name[V0Type].Data(), iCent));
        hisEffPtInclJ->Write(Form(sNameHisEffPtInJets.Data(), V0Name[V0Type].Data(), iCent));
//              hisEffPtInclJ->Write(Form(sNameHisEffPtInPC.Data(),V0Name[V0Type].Data(),iCent));
        grEffPtIncl = MakeGraphErrors(hisEffPtIncl, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersFull[iCent]);
        mgrEffPtIncl->Add(grEffPtIncl);
      }

      if(iEffPtInJets)
      {
        // eff in jets
        canEffPtInJet[iCent] = new TCanvas(Form("canEffPtInJet%d", iCent), "", iCanWidth, iCanHeight);
        mgrEffPtInJet[iCent] = new TMultiGraph();
        canEffPtRatio[iCent] = new TCanvas(Form("canEffPtRatio%d", iCent), "", iCanWidth, iCanHeight);
        mgrEffPtRatio[iCent] = new TMultiGraph();
        TLine* lineOne = new TLine(fPtInJetsXMin, 1, fPtInJetsXMax, 1);
        hisPtInJetsMCGen2D = GetHistogram2D(listMC, Form(hisPtInJetsMCGen2DName.Data(), iCent));
        if(!hisPtInJetsMCGen2D)
          return;
        spaPtInJetsMCRec3D = GetSparseD(listMC, Form(spaPtInJetsMCRec3DName.Data(), iCent));
        if(!spaPtInJetsMCRec3D)
          return;
        for(Int_t iJet = iJetMin; iJet <= iJetMax; iJet++)
        {
          if(!dArrayJetNormEff[iJet])
            continue;
//                  printf("Start: %d-%d\n",iCent,iJet);
          Int_t iJetPtFirstBin = hisPtInJetsMCGen2D->GetYaxis()->FindBin(dBinsPtJet[iJet] + dEpsilon);
          Int_t iJetPtLastBin;
          if(bOpenPtBins)
            iJetPtLastBin = hisPtInJetsMCGen2D->GetYaxis()->GetNbins() + 1;
          else
            iJetPtLastBin = hisPtInJetsMCGen2D->GetYaxis()->FindBin(dBinsPtJet[iJet + 1] - dEpsilon);
          spaPtInJetsMCRec3D->GetAxis(3)->SetRange(iJetPtFirstBin, iJetPtLastBin);
          hisPtInJetsMCGen1D = (TH1D*)hisPtInJetsMCGen2D->ProjectionX(Form("hisInJetGen_%d_%d", iCent, iJet), iJetPtFirstBin, iJetPtLastBin, "e");
          hisPtInJetsMCGen1D = (TH1D*)hisPtInJetsMCGen1D->Rebin(iNBinsPtV0InJet, "", dBinsPtV0InJet);
          if(bZeroGenErr)
          {
            for(Int_t iB = 1; iB <= hisPtInJetsMCGen1D->GetNbinsX(); iB++)
              hisPtInJetsMCGen1D->SetBinError(iB, 0);
          }
          hisPtInJetsMCRec1D = new TH1D("hisPtInJetsMCRec1D", "", iNBinsPtV0InJet, dBinsPtV0InJet);
          for(Int_t iPt = iBinPtInJetsFirst; iPt <= iBinPtInJetsLast; iPt++)
          {
            Int_t iBinPtFirst = spaPtInJetsMCRec3D->GetAxis(1)->FindBin(dBinsPtV0InJet[iPt - 1] + dEpsilon);
            Int_t iBinPtLast = spaPtInJetsMCRec3D->GetAxis(1)->FindBin(dBinsPtV0InJet[iPt] - dEpsilon);
            spaPtInJetsMCRec3D->GetAxis(1)->SetRange(iBinPtFirst, iBinPtLast);
            hisMass = (TH1D*)spaPtInJetsMCRec3D->Projection(0, "E");
            if(!hisMass)
              return;
            if(iNRebin > 1)
              hisMass->Rebin(iNRebin);
            hisMass->SetName(Form("MCPtInJets%s_C%d_J%d-Pt%d", V0Name[V0Type].Data(), iCent, iJet, iPt));

//                      Double_t dSigmaOld = MassPeakSigmaOld((dBinsPtV0AllLF[iPt-1]+dBinsPtV0AllLF[iPt])/2,V0Type);
//                      Double_t dSigMinOld = fMassV0[V0Type]-3*dSigmaOld;
//                      Double_t dSigMaxOld = fMassV0[V0Type]+3*dSigmaOld;

            Bool_t bStatusBC = kTRUE;
            BinCounterObject* binCount = new BinCounterObject(hisMass, fMassV0[V0Type], fSigmaV0[V0Type]);
            binCount->SetVerbose(bVerboseBC);
            binCount->SetDegreePolMax(0);
            bStatusBC &= binCount->SetXLimits(dBgOutMin, dBgOutMax);
            if(bFixedBC)
              bStatusBC &= binCount->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax);
//                        bStatusBC &= binCount->FixRegions(dSigMinOld,dSigMaxOld,dBgInMin,dBgInMax,dBgOutMin,dBgOutMax);
            else
              bStatusBC &= binCount->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
//                      binCount->SetConsiderBackground(0);
            binCount->SetBackgroundSubtraction(0);
//                      binCount->Plot("canBinCounter0");
            bStatusBC &= binCount->EstimateParameters();
//                      binCount->Plot("canBinCounter1");
            bStatusBC &= binCount->Run();
            binCount->Plot("canBinCounter2");
            if(!bStatusBC)
            {
              printf("Something wrong with bin counting\n");
              delete binCount;
              continue;
            }
            Double_t dSignalErr = 0;
            Double_t dSignal = binCount->GetSignalAndError(&dSignalErr);
            delete binCount;

            hisPtInJetsMCRec1D->SetBinContent(iPt, dSignal);
            hisPtInJetsMCRec1D->SetBinError(iPt, dSignalErr);
          }
//                  TH1D* hisPtInJetsMCRec1DOut = DivideHistogram(hisPtInJetsMCRec1D,dArrayJetNormEff[iJet],0,1);
//                  hisPtInJetsMCRec1DOut->Write(Form("hisPtInJetsMCRec1D_C%d_J%d",iCent,iJet));
//                  TH1D* hisPtInJetsMCGen1DOut = DivideHistogram(hisPtInJetsMCGen1D,dArrayJetNormEff[iJet],0,1);
//                  hisPtInJetsMCGen1DOut->Write(Form("hisPtInJetsMCGen1D_C%d_J%d",iCent,iJet));
//                  hisEffPtInJet = DivideHistograms1D(hisPtInJetsMCRec1DOut,hisPtInJetsMCGen1DOut,Form("hisEffPtInJet-%d-%d",iCent,iJet));
          hisEffPtInJet = DivideHistograms1D(hisPtInJetsMCRec1D, hisPtInJetsMCGen1D, Form("hisEffPtInJet-%d-%d", iCent, iJet));
          if(!hisEffPtInJet)
            return;
          hisEffPtInJet->Write(Form(sNameHisEffPtInJetsTrue.Data(), V0Name[V0Type].Data(), iCent, iJet));
          if(bOpenPtBins)
            grEffPtInJet = MakeGraphErrors(hisEffPtInJet, Form("#it{p}_{T}^{jet} > %1.f GeV/#it{c}", dBinsPtJet[iJet]), iMyColors[iJet], iMyMarkersFull[iJet]);
          else
            grEffPtInJet = MakeGraphErrors(hisEffPtInJet, Form("#it{p}_{T}^{jet}: %1.f-%1.f GeV/#it{c}", dBinsPtJet[iJet], dBinsPtJet[iJet + 1]), iMyColors[iJet], iMyMarkersFull[iJet]);
          mgrEffPtInJet[iCent]->Add(grEffPtInJet);
          hisEffPtRatio = DivideHistograms1D(hisEffPtInJet, hisEffPtInclJ, Form("hisEffPtRatio-%d-%d", iCent, iJet));
          if(!hisEffPtRatio)
            return;
          if(bOpenPtBins)
            grEffPtRatio = MakeGraphErrors(hisEffPtRatio, Form("#it{p}_{T}^{jet} > %1.f GeV/#it{c}", dBinsPtJet[iJet]), iMyColors[iJet], iMyMarkersFull[iJet], 0.5);
          else
            grEffPtRatio = MakeGraphErrors(hisEffPtRatio, Form("#it{p}_{T}^{jet}: %1.f-%1.f GeV/#it{c}", dBinsPtJet[iJet], dBinsPtJet[iJet + 1]), iMyColors[iJet], iMyMarkersFull[iJet], 0.5);
          grEffPtRatio->SetMarkerSize(1);
//                  mgrEffPtRatio[iCent]->Add(grEffPtRatio,"[]");
          mgrEffPtRatio[iCent]->Add(grEffPtRatio, "");
//                  mgrEffPtRatio[iCent]->Add(grEffPtRatio,"l");
//                  printf("End: %d-%d\n",iCent,iJet);
        }
        canEffPtInJet[iCent]->cd();
        canEffPtInJet[iCent]->SetLeftMargin(0.15);
        mgrEffPtInJet[iCent]->SetTitle(Form("%s reconstruction efficiency, in JC, c. %s;#it{p}_{T}^{gen.} (GeV/#it{c});efficiency", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrEffPtInJet[iCent]->SetMinimum(0);
        mgrEffPtInJet[iCent]->SetMaximum(0.6);
        mgrEffPtInJet[iCent]->Draw("AP0");
        mgrEffPtInJet[iCent]->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
        mgrEffPtInJet[iCent]->GetYaxis()->SetTitleOffset(2);
        legend = canEffPtInJet[iCent]->BuildLegend(0.2, 0.6, 0.55, 0.85);
        SetLegend(legend);
        canEffPtInJet[iCent]->SaveAs(Form("canEffPtInJets%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canEffPtInJet[iCent];
        delete mgrEffPtInJet[iCent];

        canEffPtRatio[iCent]->cd();
        canEffPtRatio[iCent]->SetLeftMargin(0.15);
        mgrEffPtRatio[iCent]->SetTitle(Form("%s reconstruction efficiency, in JC/inclusive, c. %s;#it{p}_{T}^{gen.} (GeV/#it{c});ratio", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrEffPtRatio[iCent]->SetMinimum(0.8);
        mgrEffPtRatio[iCent]->SetMaximum(1.3);
        mgrEffPtRatio[iCent]->Draw("AP0");
        mgrEffPtRatio[iCent]->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
//              mgrEffPtRatio[iCent]->GetXaxis()->SetLabelSize(0.1);
//              mgrEffPtRatio[iCent]->GetYaxis()->SetLabelSize(0.1);
//              mgrEffPtRatio[iCent]->GetXaxis()->SetTitleSize(0.1);
//              mgrEffPtRatio[iCent]->GetYaxis()->SetTitleSize(0.1);
        mgrEffPtRatio[iCent]->GetYaxis()->SetTitleOffset(2);
        legend = canEffPtRatio[iCent]->BuildLegend(0.5, 0.6, 0.85, 0.85);
        SetLegend(legend);
        lineOne->Draw();
        canEffPtRatio[iCent]->SaveAs(Form("canEffPtRatio%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canEffPtRatio[iCent];
        delete mgrEffPtRatio[iCent];
      }

      if(iEffEtaPtInclusive)
      {
        Int_t iNBinsPtV0InclTmp, iBinPtInclFirstTmp, iBinPtInclLastTmp;
        Double_t* dBinsPtV0InclTmp;
        TString sAppendix = "";

        if(iSwitchPtV0 == 0)
        {
          iBinPtInclFirstTmp = iBinPtInclFirst;
          iBinPtInclLastTmp = iBinPtInclLast;
          iNBinsPtV0InclTmp = iNBinsPtV0AllLF;
          dBinsPtV0InclTmp = dBinsPtV0AllLF;
          sAppendix = "-InclBins";
        }
        else if(iSwitchPtV0 == 1)
        {
          iBinPtInclFirstTmp = iBinPtInJetsFirst;
          iBinPtInclLastTmp = iBinPtInJetsLast;
          iNBinsPtV0InclTmp = iNBinsPtV0InJet;
          dBinsPtV0InclTmp = dBinsPtV0InJet;
          sAppendix = "-JetBins";
        }

        canEffEtaPtIncl[iCent] = new TCanvas(Form("canEffEtaPtIncl%d", iCent), "", iCanWidth, iCanHeight);
        mgrEffEtaPtIncl[iCent] = new TMultiGraph();
        TLine* lineCutEtaLeft = new TLine(-0.7, 0, -0.7, 0.4);
        TLine* lineCutEtaRight = new TLine(0.7, 0, 0.7, 0.4);
        // eff in eta
        hisEtaPtInclMCGen2D = GetHistogram2D(listMC, Form(hisEtaPtInclMCGen2DName.Data(), iCent));
        if(!hisEtaPtInclMCGen2D)
          return;
        spaEtaPtInclMCRec3D = GetSparseD(listMC, Form(spaEtaPtInclMCRec3DName.Data(), iCent));
        if(!spaEtaPtInclMCRec3D)
          return;
        hisEffEtaPtIncl2D = new TH2D(Form("EffEtaPt_%d", iCent), "Reconstruction efficiency", iNBinsPtV0InclTmp, dBinsPtV0InclTmp, iNBinsEtaInJet, dBinsEtaInJet);
        for(Int_t iPt = iBinPtInclFirstTmp; iPt <= iBinPtInclLastTmp; iPt++)
        {
          printf("Start Eta: %d-%d\n", iCent, iPt);
          Int_t iBinPtFirst = spaEtaPtInclMCRec3D->GetAxis(1)->FindBin(dBinsPtV0InclTmp[iPt - 1] + dEpsilon);
          Int_t iBinPtLast = spaEtaPtInclMCRec3D->GetAxis(1)->FindBin(dBinsPtV0InclTmp[iPt] - dEpsilon);
          spaEtaPtInclMCRec3D->GetAxis(1)->SetRange(iBinPtFirst, iBinPtLast);
          hisEtaPtInclMCGen1D = (TH1D*)hisEtaPtInclMCGen2D->ProjectionY(Form("hisEtaPtInclMCGen_%d_%d", iCent, iPt), iBinPtFirst, iBinPtLast, "e");
          hisEtaPtInclMCGen1D = (TH1D*)hisEtaPtInclMCGen1D->Rebin(iNBinsEtaInJet, "", dBinsEtaInJet);
          if(bZeroGenErr)
          {
            for(Int_t iB = 1; iB <= hisEtaPtInclMCGen1D->GetNbinsX(); iB++)
              hisEtaPtInclMCGen1D->SetBinError(iB, 0);
          }
          hisEtaPtInclMCRec1D = new TH1D(Form("hisEtaPtInclMCRec1D_C%d_Pt%d", iCent, iPt), "", iNBinsEtaInJet, dBinsEtaInJet);
          for(Int_t iEta = 0; iEta < iNBinsEtaInJet; iEta++)
          {
            Int_t iEtaPtFirstBin = spaEtaPtInclMCRec3D->GetAxis(2)->FindBin(dBinsEtaInJet[iEta] + dEpsilon);
            Int_t iEtaPtLastBin = spaEtaPtInclMCRec3D->GetAxis(2)->FindBin(dBinsEtaInJet[iEta + 1] - dEpsilon);
            spaEtaPtInclMCRec3D->GetAxis(2)->SetRange(iEtaPtFirstBin, iEtaPtLastBin);

            hisMass = (TH1D*)spaEtaPtInclMCRec3D->Projection(0, "E");
            if(!hisMass)
              return;
            if(iNRebin > 1)
              hisMass->Rebin(iNRebin);
            hisMass->SetName(Form("MCEtaPtInclusive%s_C%d_Pt%d_E%d", V0Name[V0Type].Data(), iCent, iPt, iEta));

            Bool_t bStatusBC = kTRUE;
            BinCounterObject* binCount = new BinCounterObject(hisMass, fMassV0[V0Type], fSigmaV0[V0Type]);
            binCount->SetVerbose(bVerboseBC);
            binCount->SetDegreePolMax(0);
            bStatusBC &= binCount->SetXLimits(dBgOutMin, dBgOutMax);
            if(bFixedBC)
              bStatusBC &= binCount->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax);
            else
              bStatusBC &= binCount->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
//                      binCount->SetConsiderBackground(0);
            binCount->SetBackgroundSubtraction(0);
//                      binCount->Plot("canBinCounter0");
            bStatusBC &= binCount->EstimateParameters();
//                      binCount->Plot("canBinCounter1");
            bStatusBC &= binCount->Run();
            binCount->Plot("canBinCounter2");
            if(!bStatusBC)
            {
              printf("Something wrong with bin counting\n");
              delete binCount;
              continue;
            }
            Double_t dSignalErr = 0;
            Double_t dSignal = binCount->GetSignalAndError(&dSignalErr);
            delete binCount;

            hisEtaPtInclMCRec1D->SetBinContent(iEta + 1, dSignal);
            hisEtaPtInclMCRec1D->SetBinError(iEta + 1, dSignalErr);
          }
          hisEffEtaPtIncl = DivideHistograms1D(hisEtaPtInclMCRec1D, hisEtaPtInclMCGen1D, Form("hisEffEtaPt-%d-%d", iCent, iPt));
          if(!hisEffEtaPtIncl)
            return;
          for(Int_t iEta = 1; iEta <= iNBinsEtaInJet; iEta++)
          {
            hisEffEtaPtIncl2D->SetBinContent(iPt, iEta, hisEffEtaPtIncl->GetBinContent(iEta));
            hisEffEtaPtIncl2D->SetBinError(iPt, iEta, hisEffEtaPtIncl->GetBinError(iEta));
          }

          grEffEtaPtIncl = MakeGraphErrors(hisEffEtaPtIncl, Form("#it{p}_{T}^{gen.}: %g-%g GeV/#it{c}", dBinsPtV0InclTmp[iPt - 1], dBinsPtV0InclTmp[iPt]), iMyColors[(iPt - iBinPtInclFirstTmp) % iNMyColors], iMyMarkersFull[(iPt - iBinPtInclFirstTmp) % iNMyMarkersFull]);
          mgrEffEtaPtIncl[iCent]->Add(grEffEtaPtIncl);
          printf("End: %d-%d\n", iCent, iPt);
        }
        canEffEtaPtIncl[iCent]->cd();
        canEffEtaPtIncl[iCent]->SetLeftMargin(0.15);
        mgrEffEtaPtIncl[iCent]->SetTitle(Form("%s reconstruction efficiency, c. %s;#it{#eta}^{gen.};efficiency", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrEffEtaPtIncl[iCent]->SetMinimum(0);
        mgrEffEtaPtIncl[iCent]->SetMaximum(0.6);
        mgrEffEtaPtIncl[iCent]->Draw("AP0");
        mgrEffEtaPtIncl[iCent]->GetXaxis()->SetLimits(-fRangeEtaV0Max, fRangeEtaV0Max);
        mgrEffEtaPtIncl[iCent]->GetYaxis()->SetTitleOffset(2);
        legend = canEffEtaPtIncl[iCent]->BuildLegend(0.2, 0.6, 0.55, 0.85);
        SetLegend(legend);
        lineCutEtaLeft->Draw();
        lineCutEtaRight->Draw();
        canEffEtaPtIncl[iCent]->SaveAs(Form("canEffEtaPtInclusive%s_C%d%s.%s", V0Name[V0Type].Data(), iCent, sAppendix.Data(), sImageSuf.Data()));
        delete canEffEtaPtIncl[iCent];
        delete mgrEffEtaPtIncl[iCent];

        dirOutSpectra->cd();
        hisEffEtaPtIncl2D->Write(Form("%s%s", Form(sNameHisEffPtEtaIncl2D.Data(), V0Name[V0Type].Data(), iCent), sAppendix.Data()));
      }

      if(iEffEtaPtDaughters)
      {
        TCanvas* canEtaV0InclGen = new TCanvas(Form("canEtaV0InclGen%d", iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrEtaV0InclGen = new TMultiGraph();
        TCanvas* canEtaV0InclRec = new TCanvas(Form("canEtaV0InclRec%d", iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrEtaV0InclRec = new TMultiGraph();
        TCanvas* canEtaDInclRec = new TCanvas(Form("canEtaDInclRec%d", iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrEtaDInclRec = new TMultiGraph();
        TCanvas* canEtaV0InclEff = new TCanvas(Form("canEtaV0InclEff%d", iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrEtaV0InclEff = new TMultiGraph();

        TCanvas* canEtaV0InJetsGen = new TCanvas(Form("canEtaV0InJetsGen%d", iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrEtaV0InJetsGen = new TMultiGraph();
        TCanvas* canEtaV0InJetsRec = new TCanvas(Form("canEtaV0InJetsRec%d", iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrEtaV0InJetsRec = new TMultiGraph();
        TCanvas* canEtaDInJetsRec = new TCanvas(Form("canEtaDInJetsRec%d", iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrEtaDInJetsRec = new TMultiGraph();
        TCanvas* canEtaV0InJetsEff = new TCanvas(Form("canEtaV0InJetsEff%d", iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrEtaV0InJetsEff = new TMultiGraph();

        TCanvas* canEtaV0RatioEff = new TCanvas(Form("canEtaV0RatioEff%d", iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrEtaV0RatioEff = new TMultiGraph();

        TCanvas* canMeanEta = new TCanvas(Form("canMeanEta%d", iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrMeanEta = new TMultiGraph();
        TGraphErrors* grMeanEta;
        TH1D* hisMeanEtaIncl;
        TH1D* hisMeanEtaInJC;
        TCanvas* canMeanEtaEff = new TCanvas(Form("canMeanEtaEff%d", iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrMeanEtaEff = new TMultiGraph();
        TGraphErrors* grMeanEtaEff;
        TH1D* hisMeanEtaEffIncl;
        TH1D* hisMeanEtaEffInJC;
        TCanvas* canMeanEtaEffRatio = new TCanvas(Form("canMeanEtaEffRatio%d", iCent), "", iCanWidth, iCanHeight);
        TGraphErrors* grMeanEtaEffRatio;
        TH1D* hisMeanEtaEffRatio;

        hisEtaPtInclMCGen2D = GetHistogram2D(listMC, Form(hisEtaPtInclMCGen2DName.Data(), iCent));
        if(!hisEtaPtInclMCGen2D)
          return;
        spaEtaPtDInclMCRec = GetSparseD(listMC, Form(spaEtaPtDInclMCRecName.Data(), iCent));
        if(!spaEtaPtDInclMCRec)
          return;
        spaEtaPtInJetsMCGen3D = GetSparseD(listMC, Form(spaEtaPtInJetsMCGen3DName.Data(), iCent));
        if(!spaEtaPtInJetsMCGen3D)
          return;
        spaEtaPtDInJetsMCRec = GetSparseD(listMC, Form(spaEtaPtDInJetsMCRecName.Data(), iCent));
        if(!spaEtaPtDInJetsMCRec)
          return;

        // restrict range of jet pT
        Int_t iBinPtJetsFirst = spaEtaPtDInJetsMCRec->GetAxis(5)->FindBin(dBinsPtJet[1] + dEpsilon);
        Int_t iBinPtJetsLast;
        if(bOpenPtBins)
          iBinPtJetsLast = spaEtaPtDInJetsMCRec->GetAxis(5)->GetNbins() + 1;
        else
          iBinPtJetsLast = spaEtaPtDInJetsMCRec->GetAxis(5)->FindBin(dBinsPtJet[1 + 1] - dEpsilon);
        spaEtaPtDInJetsMCRec->GetAxis(5)->SetRange(iBinPtJetsFirst, iBinPtJetsLast);
        spaEtaPtInJetsMCGen3D->GetAxis(2)->SetRange(iBinPtJetsFirst, iBinPtJetsLast);

        // subset corresponding to negative daughters
        spaEtaPtDInclMCRec->GetAxis(0)->SetRange(1, 1);
        spaEtaPtDInJetsMCRec->GetAxis(0)->SetRange(1, 1);

        TLine* lineOneEta = new TLine(-fRangeEtaV0Max, 1, fRangeEtaV0Max, 1);
        TLine* lineOnePt = new TLine(fPtInJetsXMin, 1, fPtInJetsXMax, 1);

        TH1D* hisPtInclGen = (TH1D*)hisEtaPtInclMCGen2D->ProjectionX(Form("hisEtaPtInclMCGen_%d", iCent), 0, -1, "e");
        if(bZeroGenErr)
        {
          for(Int_t iB = 1; iB <= hisPtInclGen->GetNbinsX(); iB++)
            hisPtInclGen->SetBinError(iB, 0);
        }
        TH1D* hisPtInclRec = (TH1D*)spaEtaPtDInclMCRec->Projection(4, "E");
        TH1D* hisPtInclEff = DivideHistograms1D(hisPtInclRec, hisPtInclGen);
        TCanvas* canEffPtIn = new TCanvas(Form("canEffPtIn_C%d", iCent), "", iCanWidth, iCanHeight);
        canEffPtIn->cd();
        TGraphErrors* grEffPtIn = MakeGraphErrors(hisPtInclEff, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersEmpty[iCent]);
        grEffPtIn->SetTitle(Form("%s: efficiency inclusive, %s;#it{p}_{T}^{gen.};efficiency", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        grEffPtIn->SetMinimum(0);
        grEffPtIn->SetMaximum(0.6);
        grEffPtIn->Draw("AP0");
        grEffPtIn->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
        canEffPtIn->SaveAs(Form("canEffDPtIn%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canEffPtIn;

        TH1D* hisPtInJetsGen = (TH1D*)spaEtaPtInJetsMCGen3D->Projection(0, "E");
        if(bZeroGenErr)
        {
          for(Int_t iB = 1; iB <= hisPtInJetsGen->GetNbinsX(); iB++)
            hisPtInJetsGen->SetBinError(iB, 0);
        }
        TH1D* hisPtInJetsRec = (TH1D*)spaEtaPtDInJetsMCRec->Projection(4, "E");
        TH1D* hisPtInJetsEff = DivideHistograms1D(hisPtInJetsRec, hisPtInJetsGen);
        TCanvas* canEffPtJC = new TCanvas(Form("canEffPtJC_C%d", iCent), "", iCanWidth, iCanHeight);
        canEffPtJC->cd();
        TGraphErrors* grEffPtJC = MakeGraphErrors(hisPtInJetsEff, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersEmpty[iCent]);
        if(bOpenPtBins)
          grEffPtJC->SetTitle(Form("%s: efficiency in JC, #it{p}_{T}^{jet} > %.0f GeV/c, %s;#it{p}_{T}^{gen.};efficiency", V0Symbol[V0Type].Data(), dBinsPtJet[1], GetCentBinLabel(iCent).Data()));
        else
          grEffPtJC->SetTitle(Form("%s: efficiency in JC, #it{p}_{T}^{jet}: %.0f-%.0f GeV/c, %s;#it{p}_{T}^{gen.};efficiency", V0Symbol[V0Type].Data(), dBinsPtJet[1], dBinsPtJet[1 + 1], GetCentBinLabel(iCent).Data()));
        grEffPtJC->SetMinimum(0);
        grEffPtJC->SetMaximum(0.6);
        grEffPtJC->Draw("AP0");
        grEffPtJC->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
        canEffPtJC->SaveAs(Form("canEffDPtJC%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canEffPtJC;

        TH1D* hisRatioEff = DivideHistograms1D(hisPtInJetsEff, hisPtInclEff);
        TCanvas* canEffDPtRatio = new TCanvas(Form("canEffPtRatio_C%d", iCent), "", iCanWidth, iCanHeight);
        canEffDPtRatio->cd();
        TGraphErrors* grEffDPtRatio = MakeGraphErrors(hisRatioEff, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersEmpty[iCent]);
        grEffDPtRatio->SetTitle(Form("%s: efficiency ratio in JC/inclusive, %s;#it{p}_{T}^{gen.};ratio", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        grEffDPtRatio->SetMinimum(0.8);
        grEffDPtRatio->SetMaximum(1.3);
        grEffDPtRatio->Draw("AP0");
        grEffDPtRatio->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
        lineOnePt->Draw();
        canEffDPtRatio->SaveAs(Form("canEffDPtRatio%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canEffDPtRatio;

//              hisPtInclGen->Scale(1./hisPtInclGen->Integral(),"width");
//              hisPtInJetsGen->Scale(1./hisPtInJetsGen->Integral(),"width");
        TH1D* hisPtRatioGen = DivideHistograms1D(hisPtInJetsGen, hisPtInclGen);

        TCanvas* canPtGen = new TCanvas(Form("canPtGen_C%d", iCent), "", iCanWidth, iCanHeight);
        TMultiGraph* mgrPtGen = new TMultiGraph();
        TGraphErrors* grPtGen;

        grPtGen = MakeGraphErrors(hisPtInclGen, Form("gen. incl. V0"), iMyColors[0], iMyMarkersEmpty[0]);
        mgrPtGen->Add(grPtGen);
        grPtGen = MakeGraphErrors(hisPtInJetsGen, Form("gen. in JC V0"), iMyColors[1], iMyMarkersEmpty[0]);
        mgrPtGen->Add(grPtGen);
        grPtGen = MakeGraphErrors(hisPtInclRec, Form("rec. incl. V0"), iMyColors[0], iMyMarkersEmpty[1]);
        mgrPtGen->Add(grPtGen);
        grPtGen = MakeGraphErrors(hisPtInJetsRec, Form("rec. in JC V0"), iMyColors[1], iMyMarkersEmpty[1]);
        mgrPtGen->Add(grPtGen);

        // negative daughters
        spaEtaPtDInclMCRec->GetAxis(0)->SetRange(1, 1);
        spaEtaPtDInJetsMCRec->GetAxis(0)->SetRange(1, 1);
        TH1D* hisPtDNegInclRec = spaEtaPtDInclMCRec->Projection(2, "E");
        TH1D* hisPtDNegInJetsRec = spaEtaPtDInJetsMCRec->Projection(2, "E");
        // positive daughters
        spaEtaPtDInclMCRec->GetAxis(0)->SetRange(2, 2);
        spaEtaPtDInJetsMCRec->GetAxis(0)->SetRange(2, 2);
        TH1D* hisPtDPosInclRec = spaEtaPtDInclMCRec->Projection(2, "E");
        TH1D* hisPtDPosInJetsRec = spaEtaPtDInJetsMCRec->Projection(2, "E");

        grPtGen = MakeGraphErrors(hisPtDNegInclRec, Form("rec. incl. d-"), iMyColors[0], iMyMarkersEmpty[2]);
        mgrPtGen->Add(grPtGen);
        grPtGen = MakeGraphErrors(hisPtDNegInJetsRec, Form("rec. in JC d-"), iMyColors[1], iMyMarkersEmpty[2]);
        mgrPtGen->Add(grPtGen);
        grPtGen = MakeGraphErrors(hisPtDPosInclRec, Form("rec. incl. d+"), iMyColors[0], iMyMarkersEmpty[3]);
        mgrPtGen->Add(grPtGen);
        grPtGen = MakeGraphErrors(hisPtDPosInJetsRec, Form("rec. in JC d+"), iMyColors[1], iMyMarkersEmpty[3]);
        mgrPtGen->Add(grPtGen);

        canPtGen->cd();
        mgrPtGen->SetTitle(Form("%s: MC #it{p}_{T}^{gen.} spectra, %s;#it{p}_{T}^{gen.};counts", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrPtGen->SetMinimum(1e1);
        mgrPtGen->SetMaximum(1e7);
        mgrPtGen->Draw("AP0");
        mgrPtGen->GetXaxis()->SetLimits(fPtInJetsXMin, 5);
        legend = canPtGen->BuildLegend(0.55, 0.65, 0.85, 0.85);
        SetLegend(legend);
        canPtGen->SetLogy();
        canPtGen->SaveAs(Form("canDPt%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));

        TCanvas* canPtRatioGen = new TCanvas(Form("canPtRatioGen_C%d", iCent), "", iCanWidth, iCanHeight);
        canPtRatioGen->cd();
        grPtGen = MakeGraphErrors(hisPtRatioGen, Form("in JC/incl."), iMyColors[iCent], iMyMarkersEmpty[iCent]);
        grPtGen->Draw("AP0");
        grPtGen->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
        legend = canPtRatioGen->BuildLegend(0.35, 0.15, 0.65, 0.35);
        SetLegend(legend);
        canPtRatioGen->SaveAs(Form("canPtRatio%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));

        Double_t* dArray = dBinsPtV0InJet;
        Int_t iNSize = iNBinsPtV0InJet;
        Int_t iCounterStyle = -1;
//              dArray = dBinsPtV0AllLF;
//              iNSize = iNBinsPtV0AllLF;
        printf("size %d\n", iNSize);

        hisMeanEtaIncl = new TH1D(Form("hisMeanEtaIncl%d", iCent), "Mean eta", iNSize, dArray);
        hisMeanEtaInJC = new TH1D(Form("hisMeanEtaInJC%d", iCent), "Mean eta", iNSize, dArray);
        hisMeanEtaEffIncl = new TH1D(Form("hisMeanEtaEffIncl%d", iCent), "Mean eta efficiency", iNSize, dArray);
        hisMeanEtaEffInJC = new TH1D(Form("hisMeanEtaEffInJC%d", iCent), "Mean eta efficiency", iNSize, dArray);
//              for (Int_t iPt=iBinPtInJetsFirst; iPt<=iBinPtInJetsLast; iPt++)
//                {
//                  if (iPt!=2 && iPt<iBinPtInJetsLast-1)
//                    continue;
        for(Int_t iPt = 2; iPt <= iNSize - 2; iPt++)
        {
          printf("%f-%f\n", dArray[iPt - 1], dArray[iPt]);
//                  if ( !( (dArray[iPt-1]>=0.6 && dArray[iPt]<=1) || (dArray[iPt-1]>=4 && dArray[iPt]<=5) ) )
//                    continue;
//                  if (iPt==3 || iPt==5 || iPt==6 || iPt>7)
//                    continue;
          iCounterStyle++;
          Int_t iPtBinFirst = hisEtaPtInclMCGen2D->GetXaxis()->FindBin(dArray[iPt - 1] + dEpsilon);
          Int_t iPtBinLast = hisEtaPtInclMCGen2D->GetXaxis()->FindBin(dArray[iPt] - dEpsilon);
          printf("%d-0: %d-%d\n", iPt, iPtBinFirst, iPtBinLast);
          /*
                            iPtBinFirst = spaEtaPtDInclMCRec->GetAxis(4)->FindBin(dArray[iPt-1]+dEpsilon);
                            iPtBinLast = spaEtaPtDInclMCRec->GetAxis(4)->FindBin(dArray[iPt]-dEpsilon);
                            printf("%d-1: %d-%d\n",iPt,iPtBinFirst,iPtBinLast);
                            iPtBinFirst = spaEtaPtInJetsMCGen3D->GetAxis(0)->FindBin(dArray[iPt-1]+dEpsilon);
                            iPtBinLast = spaEtaPtInJetsMCGen3D->GetAxis(0)->FindBin(dArray[iPt]-dEpsilon);
                            printf("%d-2: %d-%d\n",iPt,iPtBinFirst,iPtBinLast);
                            iPtBinFirst = spaEtaPtDInJetsMCRec->GetAxis(4)->FindBin(dArray[iPt-1]+dEpsilon);
                            iPtBinLast = spaEtaPtDInJetsMCRec->GetAxis(4)->FindBin(dArray[iPt]-dEpsilon);
                            printf("%d-3: %d-%d\n",iPt,iPtBinFirst,iPtBinLast);
          */
          Int_t iSwitchBinning = 1; // 0 -fine, 1 - rough, 2 - total

          spaEtaPtDInclMCRec->GetAxis(0)->SetRange(1, 1);
          spaEtaPtDInJetsMCRec->GetAxis(0)->SetRange(1, 1);

          spaEtaPtDInclMCRec->GetAxis(4)->SetRange(iPtBinFirst, iPtBinLast);
          spaEtaPtInJetsMCGen3D->GetAxis(0)->SetRange(iPtBinFirst, iPtBinLast);
          spaEtaPtDInJetsMCRec->GetAxis(4)->SetRange(iPtBinFirst, iPtBinLast);

          TH1D* hisEtaInclGen = (TH1D*)hisEtaPtInclMCGen2D->ProjectionY(Form("hisEtaPtInclMCGen_%d_%d", iCent, iPt), iPtBinFirst, iPtBinLast, "e");
          // mean eta
          Double_t dMeanEta = GetMeanValue(hisEtaInclGen);
          hisMeanEtaIncl->SetBinContent(iPt, dMeanEta);
          hisMeanEtaIncl->SetBinError(iPt, 0);

          TH1D* hisEtaInclRec = (TH1D*)spaEtaPtDInclMCRec->Projection(3, "E");
          spaEtaPtDInclMCRec->GetAxis(0)->SetRange(1, 1); // neg
          TH1D* hisEtaDNegInclRec = (TH1D*)spaEtaPtDInclMCRec->Projection(1, "E");
          spaEtaPtDInclMCRec->GetAxis(0)->SetRange(2, 2); // pos
          TH1D* hisEtaDPosInclRec = (TH1D*)spaEtaPtDInclMCRec->Projection(1, "E");
          if(iSwitchBinning == 1)
          {
            hisEtaInclGen = (TH1D*)hisEtaInclGen->Rebin(iNBinsEtaInJet, "", dBinsEtaInJet);
            hisEtaInclRec = (TH1D*)hisEtaInclRec->Rebin(iNBinsEtaInJet, "", dBinsEtaInJet);
          }
          if(iSwitchBinning == 2)
          {
            hisEtaInclGen = (TH1D*)hisEtaInclGen->Rebin(iNBinsEtaTot, "", dBinsEtaTot);
            hisEtaInclRec = (TH1D*)hisEtaInclRec->Rebin(iNBinsEtaTot, "", dBinsEtaTot);
          }
          if(bZeroGenErr)
          {
            for(Int_t iB = 1; iB <= hisEtaInclGen->GetNbinsX(); iB++)
              hisEtaInclGen->SetBinError(iB, 0);
          }

          TH1D* hisEtaInJetsGen = (TH1D*)spaEtaPtInJetsMCGen3D->Projection(1, "E");
          // mean eta
          dMeanEta = GetMeanValue(hisEtaInJetsGen);
          hisMeanEtaInJC->SetBinContent(iPt, dMeanEta);
          hisMeanEtaInJC->SetBinError(iPt, 0);
          TH1D* hisEtaInJetsRec = (TH1D*)spaEtaPtDInJetsMCRec->Projection(3, "E");
          spaEtaPtDInJetsMCRec->GetAxis(0)->SetRange(1, 1); // neg
          TH1D* hisEtaDNegInJetsRec = (TH1D*)spaEtaPtDInJetsMCRec->Projection(1, "E");
          spaEtaPtDInJetsMCRec->GetAxis(0)->SetRange(2, 2); // pos
          TH1D* hisEtaDPosInJetsRec = (TH1D*)spaEtaPtDInJetsMCRec->Projection(1, "E");
          if(iSwitchBinning == 1)
          {
            hisEtaInJetsGen = (TH1D*)hisEtaInJetsGen->Rebin(iNBinsEtaInJet, "", dBinsEtaInJet);
            hisEtaInJetsRec = (TH1D*)hisEtaInJetsRec->Rebin(iNBinsEtaInJet, "", dBinsEtaInJet);
          }
          if(iSwitchBinning == 2)
          {
            hisEtaInJetsGen = (TH1D*)hisEtaInJetsGen->Rebin(iNBinsEtaTot, "", dBinsEtaTot);
            hisEtaInJetsRec = (TH1D*)hisEtaInJetsRec->Rebin(iNBinsEtaTot, "", dBinsEtaTot);
          }
          if(bZeroGenErr)
          {
            for(Int_t iB = 1; iB <= hisEtaInJetsGen->GetNbinsX(); iB++)
              hisEtaInJetsGen->SetBinError(iB, 0);
          }

          TH1D* hisEtaInclEff = DivideHistograms1D(hisEtaInclRec, hisEtaInclGen);
          TH1D* hisEtaInJetsEff = DivideHistograms1D(hisEtaInJetsRec, hisEtaInJetsGen);
          TH1D* hisEtaRatioEff = DivideHistograms1D(hisEtaInJetsEff, hisEtaInclEff);

          // efficiency at mean eta
          Int_t iBinEffPlus;
          Int_t iBinEffMinus;
          iBinEffPlus = hisEtaInclEff->FindBin(hisMeanEtaIncl->GetBinContent(iPt));
          iBinEffMinus = hisEtaInclEff->FindBin(-hisMeanEtaIncl->GetBinContent(iPt));
          Double_t dMeanEtaEff = (hisEtaInclEff->GetBinContent(iBinEffPlus) + hisEtaInclEff->GetBinContent(iBinEffMinus)) / 2.;
          Double_t dMeanEtaEffErr = (SqrtPlus(hisEtaInclEff->GetBinError(iBinEffPlus), hisEtaInclEff->GetBinError(iBinEffMinus))) / 2.;
          hisMeanEtaEffIncl->SetBinContent(iPt, dMeanEtaEff);
          hisMeanEtaEffIncl->SetBinError(iPt, dMeanEtaEffErr);
          iBinEffPlus = hisEtaInJetsEff->FindBin(hisMeanEtaInJC->GetBinContent(iPt));
          iBinEffMinus = hisEtaInJetsEff->FindBin(-hisMeanEtaInJC->GetBinContent(iPt));
          dMeanEtaEff = (hisEtaInJetsEff->GetBinContent(iBinEffPlus) + hisEtaInJetsEff->GetBinContent(iBinEffMinus)) / 2.;
          dMeanEtaEffErr = (SqrtPlus(hisEtaInJetsEff->GetBinError(iBinEffPlus), hisEtaInJetsEff->GetBinError(iBinEffMinus))) / 2.;
          hisMeanEtaEffInJC->SetBinContent(iPt, dMeanEtaEff);
          hisMeanEtaEffInJC->SetBinError(iPt, dMeanEtaEffErr);

          hisEtaInclGen->Scale(1. / (hisEtaInclGen->Integral() == 0. ? 1 : hisEtaInclGen->Integral()), "width");
          Double_t dIntIncl = hisEtaInclRec->Integral();
          dIntIncl = (dIntIncl == 0. ? 1 : dIntIncl);
          hisEtaInclRec->Scale(1. / dIntIncl, "width");
          hisEtaDNegInclRec->Scale(1. / dIntIncl, "width");
          hisEtaDPosInclRec->Scale(1. / dIntIncl, "width");

          hisEtaInJetsGen->Scale(1. / (hisEtaInJetsGen->Integral() == 0 ? 1 : hisEtaInJetsGen->Integral()), "width");
          Double_t dIntInJets = hisEtaInJetsRec->Integral();
          dIntInJets = (dIntInJets == 0. ? 1 : dIntInJets);
          hisEtaInJetsRec->Scale(1. / dIntInJets, "width");
          hisEtaDNegInJetsRec->Scale(1. / dIntInJets, "width");
          hisEtaDPosInJetsRec->Scale(1. / dIntInJets, "width");
//                  printf("%f %f %f\n",hisEtaInclRec->Integral(0,-1),hisEtaDNegInclRec->Integral(0,-1),hisEtaDPosInclRec->Integral(0,-1));

          TGraphErrors* grEtaV0InclGen = MakeGraphErrors(hisEtaInclGen, Form("%.1f-%.1f GeV/c", dArray[iPt - 1], dArray[iPt]), iMyColors[iCounterStyle % iNMyColors], iMyMarkersEmpty[iCounterStyle % iNMyMarkersEmpty]);
          mgrEtaV0InclGen->Add(grEtaV0InclGen);
          TGraphErrors* grEtaV0InclRec = MakeGraphErrors(hisEtaInclRec, Form("%.1f-%.1f GeV/c", dArray[iPt - 1], dArray[iPt]), iMyColors[iCounterStyle % iNMyColors], iMyMarkersEmpty[iCounterStyle % iNMyMarkersEmpty]);
          mgrEtaV0InclRec->Add(grEtaV0InclRec);
          TGraphErrors* grEtaV0InclEff = MakeGraphErrors(hisEtaInclEff, Form("%.1f-%.1f GeV/c", dArray[iPt - 1], dArray[iPt]), iMyColors[iCounterStyle % iNMyColors], iMyMarkersEmpty[iCounterStyle % iNMyMarkersEmpty]);
          mgrEtaV0InclEff->Add(grEtaV0InclEff);

          TGraphErrors* grEtaDInclRec = MakeGraphErrors(hisEtaDNegInclRec, Form("d- %.1f-%.1f GeV/c", dArray[iPt - 1], dArray[iPt]), iMyColors[iCounterStyle % iNMyColors], iMyMarkersEmpty[1]);
          mgrEtaDInclRec->Add(grEtaDInclRec);
          grEtaDInclRec = MakeGraphErrors(hisEtaDPosInclRec, Form("d+ %.1f-%.1f GeV/c", dArray[iPt - 1], dArray[iPt]), iMyColors[iCounterStyle % iNMyColors], iMyMarkersEmpty[2]);
          mgrEtaDInclRec->Add(grEtaDInclRec);

          TGraphErrors* grEtaV0InJetsGen = MakeGraphErrors(hisEtaInJetsGen, Form("%.1f-%.1f GeV/c", dArray[iPt - 1], dArray[iPt]), iMyColors[iCounterStyle % iNMyColors], iMyMarkersEmpty[iCounterStyle % iNMyMarkersEmpty]);
          mgrEtaV0InJetsGen->Add(grEtaV0InJetsGen);
          TGraphErrors* grEtaV0InJetsRec = MakeGraphErrors(hisEtaInJetsRec, Form("%.1f-%.1f GeV/c", dArray[iPt - 1], dArray[iPt]), iMyColors[iCounterStyle % iNMyColors], iMyMarkersEmpty[iCounterStyle % iNMyMarkersEmpty]);
          mgrEtaV0InJetsRec->Add(grEtaV0InJetsRec);
          TGraphErrors* grEtaV0InJetsEff = MakeGraphErrors(hisEtaInJetsEff, Form("%.1f-%.1f GeV/c", dArray[iPt - 1], dArray[iPt]), iMyColors[iCounterStyle % iNMyColors], iMyMarkersEmpty[iCounterStyle % iNMyMarkersEmpty]);
          mgrEtaV0InJetsEff->Add(grEtaV0InJetsEff);

          TGraphErrors* grEtaDInJetsRec = MakeGraphErrors(hisEtaDNegInJetsRec, Form("d- %.1f-%.1f GeV/c", dArray[iPt - 1], dArray[iPt]), iMyColors[iCounterStyle % iNMyColors], iMyMarkersEmpty[1]);
          mgrEtaDInJetsRec->Add(grEtaDInJetsRec);
          grEtaDInJetsRec = MakeGraphErrors(hisEtaDPosInJetsRec, Form("d+ %.1f-%.1f GeV/c", dArray[iPt - 1], dArray[iPt]), iMyColors[iCounterStyle % iNMyColors], iMyMarkersEmpty[2]);
          mgrEtaDInJetsRec->Add(grEtaDInJetsRec);

          TGraphErrors* grEtaV0RatioEff = MakeGraphErrors(hisEtaRatioEff, Form("%.1f-%.1f GeV/c", dArray[iPt - 1], dArray[iPt]), iMyColors[iCounterStyle % iNMyColors], iMyMarkersEmpty[iCounterStyle % iNMyMarkersEmpty]);
//                  if (iPt==iBinPtInJetsFirst)
          mgrEtaV0RatioEff->Add(grEtaV0RatioEff);

//                      hisEtaInclGen = (TH1D*)hisEtaInclGen->Rebin(iNBinsEta,"",dBinsEta);
        }
        Float_t fMaxGenRecIncl = 0.025;
        Float_t fMaxGenRecInJets = 0.04;

        mgrEtaV0InclGen->SetTitle(Form("%s: generated inclusive, %s;#it{#eta}^{gen.};prob. dens.", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrEtaV0InclGen->SetMinimum(0);
//              mgrEtaV0InclGen->SetMaximum(fMaxGenRecIncl);
        canEtaV0InclGen->cd();
        mgrEtaV0InclGen->Draw("AP0");
        legend = canEtaV0InclGen->BuildLegend(0.35, 0.15, 0.65, 0.35);
        SetLegend(legend);
        canEtaV0InclGen->SaveAs(Form("canEtaV0Incl%sGen_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canEtaV0InclGen;
        delete mgrEtaV0InclGen;
        mgrEtaV0InclRec->SetTitle(Form("%s: reconstructed inclusive, %s;#it{#eta}^{gen.};prob. dens.", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrEtaV0InclRec->SetMinimum(0);
//              mgrEtaV0InclRec->SetMaximum(fMaxGenRecIncl);
        canEtaV0InclRec->cd();
        mgrEtaV0InclRec->Draw("AP0");
//              legend = canEtaV0InclRec->BuildLegend(0.35,0.15,0.65,0.35);
//              SetLegend(legend);
        canEtaV0InclRec->SaveAs(Form("canEtaV0Incl%sRec_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canEtaV0InclRec;
        delete mgrEtaV0InclRec;
        mgrEtaDInclRec->SetTitle(Form("%s: reconstructed inclusive, daughters, %s;#it{#eta}^{gen.};prob. dens.", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrEtaDInclRec->SetMinimum(0);
//              mgrEtaDInclRec->SetMaximum(fMaxGenRecIncl);
        canEtaDInclRec->cd();
        mgrEtaDInclRec->Draw("AP0");
        legend = canEtaDInclRec->BuildLegend(0.35, 0.15, 0.65, 0.35);
        SetLegend(legend);
        canEtaDInclRec->SaveAs(Form("canEtaDIncl%sRec_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canEtaDInclRec;
        delete mgrEtaDInclRec;
        mgrEtaV0InclEff->SetTitle(Form("%s: efficiency, inclusive, %s;#it{#eta}^{gen.};efficiency", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrEtaV0InclEff->SetMinimum(0);
        mgrEtaV0InclEff->SetMaximum(0.4);
        canEtaV0InclEff->cd();
        mgrEtaV0InclEff->Draw("AP0");
        canEtaV0InclEff->SaveAs(Form("canEtaV0Incl%sEff_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canEtaV0InclEff;
        delete mgrEtaV0InclEff;

        if(bOpenPtBins)
          mgrEtaV0InJetsGen->SetTitle(Form("%s: generated in JC, #it{p}_{T}^{jet} > %.0f GeV/c, %s;#it{#eta}^{gen.};prob. dens.", V0Symbol[V0Type].Data(), dBinsPtJet[1], GetCentBinLabel(iCent).Data()));
        else
          mgrEtaV0InJetsGen->SetTitle(Form("%s: generated in JC, #it{p}_{T}^{jet}: %.0f-%.0f GeV/c, %s;#it{#eta}^{gen.};prob. dens.", V0Symbol[V0Type].Data(), dBinsPtJet[1], dBinsPtJet[1 + 1], GetCentBinLabel(iCent).Data()));
        mgrEtaV0InJetsGen->SetMinimum(0);
//              mgrEtaV0InJetsGen->SetMaximum(fMaxGenRecInJets);
        canEtaV0InJetsGen->cd();
        mgrEtaV0InJetsGen->Draw("AP0");
        legend = canEtaV0InJetsGen->BuildLegend(0.35, 0.15, 0.65, 0.35);
        SetLegend(legend);
        canEtaV0InJetsGen->SaveAs(Form("canEtaV0InJets%sGen_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canEtaV0InJetsGen;
        delete mgrEtaV0InJetsGen;
        if(bOpenPtBins)
          mgrEtaV0InJetsRec->SetTitle(Form("%s: reconstructed in JC, #it{p}_{T}^{jet} > %.0f GeV/c, %s;#it{#eta}^{gen.};prob. dens.", V0Symbol[V0Type].Data(), dBinsPtJet[1], GetCentBinLabel(iCent).Data()));
        else
          mgrEtaV0InJetsRec->SetTitle(Form("%s: reconstructed in JC, #it{p}_{T}^{jet}: %.0f-%.0f GeV/c, %s;#it{#eta}^{gen.};prob. dens.", V0Symbol[V0Type].Data(), dBinsPtJet[1], dBinsPtJet[1 + 1], GetCentBinLabel(iCent).Data()));
        mgrEtaV0InJetsRec->SetMinimum(0);
//              mgrEtaV0InJetsRec->SetMaximum(fMaxGenRecInJets);
        canEtaV0InJetsRec->cd();
        mgrEtaV0InJetsRec->Draw("AP0");
        canEtaV0InJetsRec->SaveAs(Form("canEtaV0InJets%sRec_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canEtaV0InJetsRec;
        delete mgrEtaV0InJetsRec;
        if(bOpenPtBins)
          mgrEtaDInJetsRec->SetTitle(Form("%s: reconstructed in JC, daughters, #it{p}_{T}^{jet} > %.0f GeV/c, %s;#it{#eta}^{gen.};prob. dens.", V0Symbol[V0Type].Data(), dBinsPtJet[1], GetCentBinLabel(iCent).Data()));
        else
          mgrEtaDInJetsRec->SetTitle(Form("%s: reconstructed in JC, daughters, #it{p}_{T}^{jet}: %.0f-%.0f GeV/c, %s;#it{#eta}^{gen.};prob. dens.", V0Symbol[V0Type].Data(), dBinsPtJet[1], dBinsPtJet[1 + 1], GetCentBinLabel(iCent).Data()));
        mgrEtaDInJetsRec->SetMinimum(0);
//              mgrEtaDInJetsRec->SetMaximum(fMaxGenRecInJets);
        canEtaDInJetsRec->cd();
        mgrEtaDInJetsRec->Draw("AP0");
        legend = canEtaDInJetsRec->BuildLegend(0.35, 0.15, 0.65, 0.35);
        SetLegend(legend);
        canEtaDInJetsRec->SaveAs(Form("canEtaDInJets%sRec_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canEtaDInJetsRec;
        delete mgrEtaDInJetsRec;
        if(bOpenPtBins)
          mgrEtaV0InJetsEff->SetTitle(Form("%s: efficiency, in JC, #it{p}_{T}^{jet} > %.0f GeV/c, %s;#it{#eta}^{gen.};efficiency", V0Symbol[V0Type].Data(), dBinsPtJet[1], GetCentBinLabel(iCent).Data()));
        else
          mgrEtaV0InJetsEff->SetTitle(Form("%s: efficiency, in JC, #it{p}_{T}^{jet}: %.0f-%.0f GeV/c, %s;#it{#eta}^{gen.};efficiency", V0Symbol[V0Type].Data(), dBinsPtJet[1], dBinsPtJet[1 + 1], GetCentBinLabel(iCent).Data()));
        mgrEtaV0InJetsEff->SetMinimum(0);
        mgrEtaV0InJetsEff->SetMaximum(0.4);
        canEtaV0InJetsEff->cd();
        mgrEtaV0InJetsEff->Draw("AP0");
        canEtaV0InJetsEff->SaveAs(Form("canEtaV0InJets%sEff_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canEtaV0InJetsEff;
        delete mgrEtaV0InJetsEff;

        mgrEtaV0RatioEff->SetTitle(Form("%s: efficiency ratio, in JC/inclusive, %s;#it{#eta}^{gen.};ratio", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrEtaV0RatioEff->SetMinimum(0.0);
        mgrEtaV0RatioEff->SetMaximum(2);
        canEtaV0RatioEff->cd();
        mgrEtaV0RatioEff->Draw("AP0");
        legend = canEtaV0RatioEff->BuildLegend(0.35, 0.65, 0.65, 0.85);
        SetLegend(legend);
        lineOneEta->Draw();
        canEtaV0RatioEff->SaveAs(Form("canEtaV0Ratio%sEff_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canEtaV0RatioEff;
        delete mgrEtaV0RatioEff;

        grMeanEta = MakeGraphErrors(hisMeanEtaIncl, "inclusive", iMyColors[0], iMyMarkersFull[0]);
        mgrMeanEta->Add(grMeanEta);
        grMeanEta = MakeGraphErrors(hisMeanEtaInJC, "in JC", iMyColors[1], iMyMarkersFull[2]);
        mgrMeanEta->Add(grMeanEta);
        mgrMeanEta->SetTitle(Form("%s: mean |#it{#eta}^{gen.}|, %s;#it{p}_{T}^{gen.} (GeV/#it{c});mean |#it{#eta}^{gen.}|", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrMeanEta->SetMinimum(0.0);
        mgrMeanEta->SetMaximum(0.4);
        canMeanEta->cd();
        canMeanEta->SetLeftMargin(0.15);
        mgrMeanEta->Draw("AP0");
        mgrMeanEta->GetYaxis()->SetTitleOffset(2);
        mgrMeanEta->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
        legend = canMeanEta->BuildLegend(0.55, 0.15, 0.85, 0.35);
        SetLegend(legend);
        canMeanEta->SaveAs(Form("canMeanEta%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canMeanEta;
        delete mgrMeanEta;
        grMeanEtaEff = MakeGraphErrors(hisMeanEtaEffIncl, "inclusive", iMyColors[0], iMyMarkersFull[0]);
        mgrMeanEtaEff->Add(grMeanEtaEff);
        grMeanEtaEff = MakeGraphErrors(hisMeanEtaEffInJC, "in JC", iMyColors[1], iMyMarkersFull[2]);
        mgrMeanEtaEff->Add(grMeanEtaEff);
        mgrMeanEtaEff->SetTitle(Form("%s: efficiency at mean |#it{#eta}^{gen.}|, %s;#it{p}_{T}^{gen.} (GeV/#it{c});efficiency at mean |#it{#eta}^{gen.}|", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrMeanEtaEff->SetMinimum(0.0);
        mgrMeanEtaEff->SetMaximum(0.4);
        canMeanEtaEff->cd();
        canMeanEtaEff->SetLeftMargin(0.15);
        mgrMeanEtaEff->Draw("AP0");
        mgrMeanEtaEff->GetYaxis()->SetTitleOffset(2);
        mgrMeanEtaEff->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
        legend = canMeanEtaEff->BuildLegend(0.55, 0.15, 0.85, 0.35);
        SetLegend(legend);
        canMeanEtaEff->SaveAs(Form("canMeanEtaEff%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canMeanEtaEff;
        delete mgrMeanEtaEff;
        hisMeanEtaEffRatio = DivideHistograms1D(hisMeanEtaEffInJC, hisMeanEtaEffIncl);
        grMeanEtaEffRatio = MakeGraphErrors(hisMeanEtaEffRatio, "in JC/inclusive", iMyColors[0], iMyMarkersFull[0]);
        grMeanEtaEffRatio->SetTitle(Form("%s: ratio of efficiencies at mean |#it{#eta}^{gen.}|, in JC/inclusive, %s;#it{p}_{T}^{gen.} (GeV/#it{c});ratio", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        grMeanEtaEffRatio->SetMinimum(0.8);
        grMeanEtaEffRatio->SetMaximum(1.3);
        canMeanEtaEffRatio->cd();
        canMeanEtaEffRatio->SetLeftMargin(0.15);
        grMeanEtaEffRatio->Draw("AP0");
        grMeanEtaEffRatio->GetYaxis()->SetTitleOffset(2);
        grMeanEtaEffRatio->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
        TLine* lineOneEffRatio = new TLine(fPtInJetsXMin, 1, fPtInJetsXMax, 1);
        lineOneEffRatio->Draw();
//              legend = canMeanEtaEff->BuildLegend(0.55,0.15,0.85,0.35);
//              SetLegend(legend);
        canMeanEtaEffRatio->SaveAs(Form("canMeanEtaEffRatio%s_C%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canMeanEtaEffRatio;
        delete grMeanEtaEffRatio;
      }

      if(iEffPtEtaInclusive)
      {
        hisPtEtaInclMCGen2D = GetHistogram2D(listMC, Form(hisPtEtaInclMCGen2DName.Data(), iCent));
        if(!hisPtEtaInclMCGen2D)
          return;
        spaPtEtaInclMCRec3D = GetSparseD(listMC, Form(spaPtEtaInclMCRec3DName.Data(), iCent));
        if(!spaPtEtaInclMCRec3D)
          return;
        printf("iEffPtEtaIncl: C%d: Start\n", iCent);
        canEffPtEtaIncl = new TCanvas(Form("canEffPtEtaIncl_%d", iCent), "", iCanWidth, iCanHeight);
        mgrEffPtEtaIncl = new TMultiGraph();
        canEffPtEtaInclPos = new TCanvas(Form("canEffPtEtaInclPos_%d", iCent), "", iCanWidth, iCanHeight);
        mgrEffPtEtaInclPos = new TMultiGraph();
        Int_t iBinEtaMax = iNBinsEtaInJet - 1;
        if(bEffEtaWindow)
          iBinEtaMax = iNBinsEtaInJet / 2;
        for(Int_t iEta = 1; iEta < iBinEtaMax; iEta++)
        {
          Int_t iEtaPtFirstBin = spaPtEtaInclMCRec3D->GetAxis(2)->FindBin(dBinsEtaInJet[iEta] + dEpsilon);
          Int_t iEtaPtLastBin = spaPtEtaInclMCRec3D->GetAxis(2)->FindBin(dBinsEtaInJet[iEta + 1] - dEpsilon);
          if(bEffEtaWindow)
            iEtaPtLastBin = spaPtEtaInclMCRec3D->GetAxis(2)->FindBin(-dBinsEtaInJet[iEta] - dEpsilon);
          spaPtEtaInclMCRec3D->GetAxis(2)->SetRange(iEtaPtFirstBin, iEtaPtLastBin);
          hisPtEtaInclMCGen1D = (TH1D*)hisPtEtaInclMCGen2D->ProjectionX(Form("hisPtEtaInclMCGen_%d_%d", iCent, iEta), iEtaPtFirstBin, iEtaPtLastBin, "e");
          hisPtEtaInclMCGen1D = (TH1D*)hisPtEtaInclMCGen1D->Rebin(iNBinsPtV0InJet, "", dBinsPtV0InJet);
          if(bZeroGenErr)
          {
            for(Int_t iB = 1; iB <= hisPtEtaInclMCGen1D->GetNbinsX(); iB++)
              hisPtEtaInclMCGen1D->SetBinError(iB, 0);
          }
//                  hisPtEtaInclMCGen1D->Write(Form(sNameHisGenPtEtaIncl.Data(),V0Name[V0Type].Data(),iCent,iEta));
          hisPtEtaInclMCRec1D = new TH1D(Form("hisPtEtaInclMCRec1D_C%d_E%d", iCent, iEta), "", iNBinsPtV0InJet, dBinsPtV0InJet);
          for(Int_t iPt = iBinPtInJetsFirst - 1; iPt <= iBinPtInJetsLast - 1; iPt++)
          {
            Int_t iPtBinFirst = spaPtEtaInclMCRec3D->GetAxis(1)->FindBin(dBinsPtV0InJet[iPt] + dEpsilon);
            Int_t iPtBinLast = spaPtEtaInclMCRec3D->GetAxis(1)->FindBin(dBinsPtV0InJet[iPt + 1] - dEpsilon);
            spaPtEtaInclMCRec3D->GetAxis(1)->SetRange(iPtBinFirst, iPtBinLast);
            hisMass = (TH1D*)spaPtEtaInclMCRec3D->Projection(0, "E");
            if(!hisMass)
              return;
            if(iNRebin > 1)
              hisMass->Rebin(iNRebin);
            hisMass->SetName(Form("MCPtEtaIncl%s_C%d_E%d-Pt%d", V0Name[V0Type].Data(), iCent, iEta, iPt));

            Bool_t bStatusBC = kTRUE;
            BinCounterObject* binCount = new BinCounterObject(hisMass, fMassV0[V0Type], fSigmaV0[V0Type]);
            binCount->SetVerbose(bVerboseBC);
            binCount->SetDegreePolMax(0);
            bStatusBC &= binCount->SetXLimits(dBgOutMin, dBgOutMax);
            if(bFixedBC)
              bStatusBC &= binCount->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax);
            else
              bStatusBC &= binCount->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
//                      binCount->SetConsiderBackground(0);
            binCount->SetBackgroundSubtraction(0);
//                      binCount->Plot("canBinCounter0");
            bStatusBC &= binCount->EstimateParameters();
//                      binCount->Plot("canBinCounter1");
            bStatusBC &= binCount->Run();
            binCount->Plot("canBinCounter2");
            if(!bStatusBC)
            {
              printf("Something wrong with bin counting\n");
              delete binCount;
              continue;
            }
            Double_t dSignalErr = 0;
            Double_t dSignal = binCount->GetSignalAndError(&dSignalErr);
            delete binCount;

            hisPtEtaInclMCRec1D->SetBinContent(iPt + 1, dSignal);
            hisPtEtaInclMCRec1D->SetBinError(iPt + 1, dSignalErr);
          }
          hisEffPtEtaIncl = DivideHistograms1D(hisPtEtaInclMCRec1D, hisPtEtaInclMCGen1D, Form("hisEffPtEtaIncl_C%d_E%d", iCent, iEta));
          if(!hisEffPtEtaIncl)
            return;
          hisEffPtEtaIncl->Write(Form(sNameHisEffPtEtaIncl.Data(), V0Name[V0Type].Data(), iCent, iEta));
          if(iEta < iNBinsEtaInJet / 2)
          {
            if(!bEffEtaWindow)
              grEffPtEtaIncl = MakeGraphErrors(hisEffPtEtaIncl, Form("#it{#eta}: %.2f-%.2f", dBinsEtaInJet[iEta], dBinsEtaInJet[iEta + 1]), iMyColors[(iEta - 1) % (iNBinsEtaInJet / 2) % iNMyColors], iMyMarkersFull[(iEta - 1) % (iNBinsEtaInJet / 2) % iNMyMarkersFull]);
            else
              grEffPtEtaIncl = MakeGraphErrors(hisEffPtEtaIncl, Form("|#it{#eta}| < %.2f", -dBinsEtaInJet[iEta]), iMyColors[(iEta - 1) % iNMyColors], iMyMarkersFull[(iEta - 1) % iNMyMarkersFull]);
            mgrEffPtEtaIncl->Add(grEffPtEtaIncl);
          }
          else
          {
            if(!bEffEtaWindow)
              grEffPtEtaIncl = MakeGraphErrors(hisEffPtEtaIncl, Form("#it{#eta}: %.2f-%.2f", dBinsEtaInJet[iEta], dBinsEtaInJet[iEta + 1]), iMyColors[(iNBinsEtaInJet - iEta - 1 - 1) % (iNBinsEtaInJet / 2) % iNMyColors], iMyMarkersFull[(iNBinsEtaInJet - iEta - 1) % (iNBinsEtaInJet / 2) % iNMyMarkersFull]);
            else
              grEffPtEtaIncl = MakeGraphErrors(hisEffPtEtaIncl, Form("|#it{#eta}| < %.2f", -dBinsEtaInJet[iEta]), iMyColors[(iEta - 1) % iNMyColors], iMyMarkersFull[(iEta - 1) % iNMyMarkersFull]);
            mgrEffPtEtaInclPos->Add(grEffPtEtaIncl);
          }
        }
        canEffPtEtaIncl->cd();
        canEffPtEtaIncl->SetLeftMargin(0.15);
        mgrEffPtEtaIncl->SetTitle(Form("%s reconstruction efficiency, inclusive, c. %s;#it{p}_{T}^{gen.} (GeV/#it{c});efficiency", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
        mgrEffPtEtaIncl->SetMinimum(0);
        mgrEffPtEtaIncl->SetMaximum(0.6);
        mgrEffPtEtaIncl->Draw("AP0");
        mgrEffPtEtaIncl->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
        mgrEffPtEtaIncl->GetYaxis()->SetTitleOffset(2);
        legend = canEffPtEtaIncl->BuildLegend(0.2, 0.6, 0.55, 0.85);
        SetLegend(legend);
        canEffPtEtaIncl->SaveAs(Form("canEffPtEtaIncl%s_C%d_Neg.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canEffPtEtaIncl;
        delete mgrEffPtEtaIncl;

        if(!bEffEtaWindow)
        {
          canEffPtEtaInclPos->cd();
          canEffPtEtaInclPos->SetLeftMargin(0.15);
          mgrEffPtEtaInclPos->SetTitle(Form("%s reconstruction efficiency, inclusive, c. %s;#it{p}_{T}^{gen.} (GeV/#it{c});efficiency", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data()));
          mgrEffPtEtaInclPos->SetMinimum(0);
          mgrEffPtEtaInclPos->SetMaximum(0.6);
          mgrEffPtEtaInclPos->Draw("AP0");
          mgrEffPtEtaInclPos->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
          mgrEffPtEtaInclPos->GetYaxis()->SetTitleOffset(2);
          legend = canEffPtEtaInclPos->BuildLegend(0.2, 0.6, 0.55, 0.85);
          SetLegend(legend);
          canEffPtEtaInclPos->SaveAs(Form("canEffPtEtaIncl%s_C%d_Pos.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        }
        delete canEffPtEtaInclPos;
        delete mgrEffPtEtaInclPos;

        printf("iEffPtEtaIncl: C%d: End\n", iCent);
      }

      if(iEffPtEtaInJets)
      {
        spaPtEtaInJetMCGen3D = GetSparseD(listMC, Form(spaPtEtaInJetMCGen3DName.Data(), iCent));
        if(!spaPtEtaInJetMCGen3D)
          return;
        spaPtEtaInJetMCRec4D = GetSparseD(listMC, Form(spaPtEtaInJetMCRec4DName.Data(), iCent));
        if(!spaPtEtaInJetMCRec4D)
          return;
        for(Int_t iJet = iJetMin; iJet <= iJetMax; iJet++)
        {
          printf("iEffPtEtaInJets: C%d J%d: Start\n", iCent, iJet);
          canEffPtEtaInJet = new TCanvas(Form("canEffPtEtaInJet_%d_%d", iCent, iJet), "", iCanWidth, iCanHeight);
          mgrEffPtEtaInJet = new TMultiGraph();
          canEffPtEtaInJetPos = new TCanvas(Form("canEffPtEtaInJetPos_%d_%d", iCent, iJet), "", iCanWidth, iCanHeight);
          mgrEffPtEtaInJetPos = new TMultiGraph();
          canEffPtEtaRatio = new TCanvas(Form("canEffPtEtaRatio_%d_%d", iCent, iJet), "", iCanWidth, iCanHeight);
          mgrEffPtEtaRatio = new TMultiGraph();
          canEffPtEtaRatioPos = new TCanvas(Form("canEffPtEtaRatioPos_%d_%d", iCent, iJet), "", iCanWidth, iCanHeight);
          mgrEffPtEtaRatioPos = new TMultiGraph();
          canGenPtEtaRatio = new TCanvas(Form("canGenPtEtaRatio_%d_%d", iCent, iJet), "", iCanWidth, iCanHeight);
          mgrGenPtEtaRatio = new TMultiGraph();
          canGenPtEtaRatioPos = new TCanvas(Form("canGenPtEtaRatioPos_%d_%d", iCent, iJet), "", iCanWidth, iCanHeight);
          mgrGenPtEtaRatioPos = new TMultiGraph();
          Int_t iJetPtFirstBin = spaPtEtaInJetMCRec4D->GetAxis(3)->FindBin(dBinsPtJet[iJet] + dEpsilon);
          Int_t iJetPtLastBin;
          if(bOpenPtBins)
            iJetPtLastBin = spaPtEtaInJetMCRec4D->GetAxis(3)->GetNbins() + 1;
          else
            iJetPtLastBin = spaPtEtaInJetMCRec4D->GetAxis(3)->FindBin(dBinsPtJet[iJet + 1] - dEpsilon);
          spaPtEtaInJetMCRec4D->GetAxis(3)->SetRange(iJetPtFirstBin, iJetPtLastBin);
          spaPtEtaInJetMCGen3D->GetAxis(2)->SetRange(iJetPtFirstBin, iJetPtLastBin);
          Int_t iBinEtaMax = iNBinsEtaInJet - 1;
          if(bEffEtaWindow)
            iBinEtaMax = iNBinsEtaInJet / 2;
          for(Int_t iEta = 1; iEta < iBinEtaMax; iEta++)
          {
            Int_t iEtaPtFirstBin = spaPtEtaInJetMCRec4D->GetAxis(2)->FindBin(dBinsEtaInJet[iEta] + dEpsilon);
            Int_t iEtaPtLastBin = spaPtEtaInJetMCRec4D->GetAxis(2)->FindBin(dBinsEtaInJet[iEta + 1] - dEpsilon);
            if(bEffEtaWindow)
              iEtaPtLastBin = spaPtEtaInJetMCRec4D->GetAxis(2)->FindBin(-dBinsEtaInJet[iEta] - dEpsilon);
            spaPtEtaInJetMCRec4D->GetAxis(2)->SetRange(iEtaPtFirstBin, iEtaPtLastBin);
            spaPtEtaInJetMCGen3D->GetAxis(1)->SetRange(iEtaPtFirstBin, iEtaPtLastBin);
            /*
                                  // delta eta
                                  Int_t iDEtaBinSel = 4;
                                  Bool_t bEffDEtaWindow = 0;
                                  Int_t iDEtaPtFirstBin = spaPtEtaInJetMCRec4D->GetAxis(4)->FindBin(dBinsDeltaEtaInJet[iDEtaBinSel]+dEpsilon);
                                  Int_t iDEtaPtLastBin = spaPtEtaInJetMCRec4D->GetAxis(4)->FindBin(dBinsDeltaEtaInJet[iDEtaBinSel+1]-dEpsilon);
                                  if (bEffDEtaWindow)
                                    iDEtaPtLastBin = spaPtEtaInJetMCRec4D->GetAxis(4)->FindBin(-dBinsDeltaEtaInJet[iDEtaBinSel]-dEpsilon);
                                  printf("Delta eta restricted to: %f - %f\n",dBinsDeltaEtaInJet[iDEtaBinSel],(bEffDEtaWindow ? -dBinsDeltaEtaInJet[iDEtaBinSel] : dBinsDeltaEtaInJet[iDEtaBinSel+1]));
                                  spaPtEtaInJetMCRec4D->GetAxis(4)->SetRange(iDEtaPtFirstBin,iDEtaPtLastBin);
                                  spaPtEtaInJetMCGen3D->GetAxis(3)->SetRange(iDEtaPtFirstBin,iDEtaPtLastBin);
            */
            hisPtEtaInJetMCGen1D = (TH1D*)spaPtEtaInJetMCGen3D->Projection(0, "E");
            hisPtEtaInJetMCGen1D = (TH1D*)hisPtEtaInJetMCGen1D->Rebin(iNBinsPtV0InJet, "", dBinsPtV0InJet);
            if(bZeroGenErr)
            {
              for(Int_t iB = 1; iB <= hisPtEtaInJetMCGen1D->GetNbinsX(); iB++)
                hisPtEtaInJetMCGen1D->SetBinError(iB, 0);
            }
            hisPtEtaInJetMCRec1D = new TH1D(Form("hisPtEtaInJetMCRec1D_C%d_J%d_E%d", iCent, iJet, iEta), "", iNBinsPtV0InJet, dBinsPtV0InJet);
            for(Int_t iPt = iBinPtInJetsFirst - 1; iPt <= iBinPtInJetsLast - 1; iPt++)
            {
              Int_t iPtBinFirst = spaPtEtaInJetMCRec4D->GetAxis(1)->FindBin(dBinsPtV0InJet[iPt] + dEpsilon);
              Int_t iPtBinLast = spaPtEtaInJetMCRec4D->GetAxis(1)->FindBin(dBinsPtV0InJet[iPt + 1] - dEpsilon);
              spaPtEtaInJetMCRec4D->GetAxis(1)->SetRange(iPtBinFirst, iPtBinLast);
              hisMass = (TH1D*)spaPtEtaInJetMCRec4D->Projection(0, "E");
              if(!hisMass)
                return;
              if(iNRebin > 1)
                hisMass->Rebin(iNRebin);
              hisMass->SetName(Form("MCPtEtaInJet%s_C%d_J%d_E%d-Pt%d", V0Name[V0Type].Data(), iCent, iJet, iEta, iPt));

              Bool_t bStatusBC = kTRUE;
              BinCounterObject* binCount = new BinCounterObject(hisMass, fMassV0[V0Type], fSigmaV0[V0Type]);
              binCount->SetVerbose(bVerboseBC);
              binCount->SetDegreePolMax(0);
              bStatusBC &= binCount->SetXLimits(dBgOutMin, dBgOutMax);
              if(bFixedBC)
                bStatusBC &= binCount->FixRegions(dSigMin, dSigMax, dBgInMin, dBgInMax, dBgOutMin, dBgOutMax);
              else
                bStatusBC &= binCount->SetNSigmas(dSigmaSig, dSigmaBgIn, dSigmaBgOut);
//                          binCount->SetConsiderBackground(0);
              binCount->SetBackgroundSubtraction(0);
//                          binCount->Plot("canBinCounter0");
              bStatusBC &= binCount->EstimateParameters();
//                          binCount->Plot("canBinCounter1");
              bStatusBC &= binCount->Run();
              binCount->Plot("canBinCounter2");
              if(!bStatusBC)
              {
                printf("Something wrong with bin counting\n");
                delete binCount;
                continue;
              }
              Double_t dSignalErr = 0;
              Double_t dSignal = binCount->GetSignalAndError(&dSignalErr);
              delete binCount;

              hisPtEtaInJetMCRec1D->SetBinContent(iPt + 1, dSignal);
              hisPtEtaInJetMCRec1D->SetBinError(iPt + 1, dSignalErr);
            }
            hisEffPtEtaInJet = DivideHistograms1D(hisPtEtaInJetMCRec1D, hisPtEtaInJetMCGen1D, Form("hisEffPtEtaInJet_C%d_J%d_E%d", iCent, iJet, iEta));
            if(!hisEffPtEtaInJet)
              return;
            hisEffPtEtaInJet->Write(Form(sNameHisEffPtEtaInJets.Data(), V0Name[V0Type].Data(), iCent, iJet, iEta));
            hisEffPtEtaIncl = (TH1D*)dirOutSpectra->Get(Form(sNameHisEffPtEtaIncl.Data(), V0Name[V0Type].Data(), iCent, iEta));
            if(!hisEffPtEtaIncl)
            {
              printf("Loading inclusive pt-eta eff failed. (Error)\n");
              return;
            }
            hisEffPtEtaRatio = DivideHistograms1D(hisEffPtEtaInJet, hisEffPtEtaIncl, Form("hisEffPtRatio_C%d_J%d_E%d", iCent, iJet, iEta));
            if(!hisEffPtEtaRatio)
              return;
            hisPtEtaInclMCGen1D = (TH1D*)dirOutSpectra->Get(Form(sNameHisGenPtEtaIncl.Data(), V0Name[V0Type].Data(), iCent, iEta));
            if(!hisPtEtaInclMCGen1D)
            {
              printf("Loading inclusive pt-eta gen failed. (Error)\n");
              return;
            }
            hisGenPtEtaRatio = DivideHistograms1D(hisPtEtaInJetMCGen1D, hisPtEtaInclMCGen1D, Form("hisGenPtRatio_C%d_J%d_E%d", iCent, iJet, iEta));
            if(!hisGenPtEtaRatio)
              return;
//            if (iEta==1)
//              hisGenPtEtaRatio->Write(Form("hisGenPtRatio_C%d_J%d_E%d",iCent,iJet,iEta));
//            TH1D* hisGenPtEtaRatioRef = (TH1D*)dirOutSpectra->Get(Form("hisGenPtRatio_C%d_J%d_E%d",iCent,iJet,1));
//            if (!hisGenPtEtaRatioRef)
//              return;
//            hisGenPtEtaRatio = DivideHistograms1D(hisGenPtEtaRatio,hisGenPtEtaRatioRef);
//            if (!hisGenPtEtaRatio)
//              return;
//            hisGenPtEtaRatio->Scale(1./hisGenPtEtaRatio->GetBinContent(9));
            if(iEta < iNBinsEtaInJet / 2)
            {
              if(!bEffEtaWindow)
              {
                grEffPtEtaInJet = MakeGraphErrors(hisEffPtEtaInJet, Form("#it{#eta}^{gen.}: %.2f-%.2f", dBinsEtaInJet[iEta], dBinsEtaInJet[iEta + 1]), iMyColors[(iEta - 1) % (iNBinsEtaInJet / 2) % iNMyColors], iMyMarkersFull[(iEta - 1) % (iNBinsEtaInJet / 2) % iNMyMarkersFull]);
                grEffPtEtaRatio = MakeGraphErrors(hisEffPtEtaRatio, Form("#it{#eta}^{gen.}: %.2f-%.2f", dBinsEtaInJet[iEta], dBinsEtaInJet[iEta + 1]), iMyColors[(iEta - 1) % (iNBinsEtaInJet / 2) % iNMyColors], iMyMarkersFull[(iEta - 1) % (iNBinsEtaInJet / 2) % iNMyMarkersFull]);
                grGenPtEtaRatio = MakeGraphErrors(hisGenPtEtaRatio, Form("#it{#eta}^{gen.}: %.2f-%.2f", dBinsEtaInJet[iEta], dBinsEtaInJet[iEta + 1]), iMyColors[(iEta - 1) % (iNBinsEtaInJet / 2) % iNMyColors], iMyMarkersFull[(iEta - 1) % (iNBinsEtaInJet / 2) % iNMyMarkersFull]);
              }
              else
              {
                grEffPtEtaInJet = MakeGraphErrors(hisEffPtEtaInJet, Form("|#it{#eta}^{gen.}| < %.2f", -dBinsEtaInJet[iEta]), iMyColors[(iEta - 1) % iNMyColors], iMyMarkersFull[(iEta - 1) % iNMyMarkersFull]);
                grEffPtEtaRatio = MakeGraphErrors(hisEffPtEtaRatio, Form("|#it{#eta}^{gen.}| < %.2f", -dBinsEtaInJet[iEta]), iMyColors[(iEta - 1) % iNMyColors], iMyMarkersFull[(iEta - 1) % iNMyMarkersFull]);
                grGenPtEtaRatio = MakeGraphErrors(hisGenPtEtaRatio, Form("|#it{#eta}^{gen.}| < %.2f", -dBinsEtaInJet[iEta]), iMyColors[(iEta - 1) % iNMyColors], iMyMarkersFull[(iEta - 1) % iNMyMarkersFull]);
              }
              mgrEffPtEtaInJet->Add(grEffPtEtaInJet);
              mgrEffPtEtaRatio->Add(grEffPtEtaRatio);
              mgrGenPtEtaRatio->Add(grGenPtEtaRatio);
            }
            else
            {
              if(!bEffEtaWindow)
              {
                grEffPtEtaInJet = MakeGraphErrors(hisEffPtEtaInJet, Form("#it{#eta}^{gen.}: %.2f-%.2f", dBinsEtaInJet[iEta], dBinsEtaInJet[iEta + 1]), iMyColors[(iNBinsEtaInJet - iEta - 1 - 1) % (iNBinsEtaInJet / 2) % iNMyColors], iMyMarkersFull[(iNBinsEtaInJet - iEta - 1) % (iNBinsEtaInJet / 2) % iNMyMarkersFull]);
                grEffPtEtaRatio = MakeGraphErrors(hisEffPtEtaRatio, Form("#it{#eta}^{gen.}: %.2f-%.2f", dBinsEtaInJet[iEta], dBinsEtaInJet[iEta + 1]), iMyColors[(iNBinsEtaInJet - iEta - 1 - 1) % (iNBinsEtaInJet / 2) % iNMyColors], iMyMarkersFull[(iNBinsEtaInJet - iEta - 1) % (iNBinsEtaInJet / 2) % iNMyMarkersFull]);
                grGenPtEtaRatio = MakeGraphErrors(hisGenPtEtaRatio, Form("#it{#eta}^{gen.}: %.2f-%.2f", dBinsEtaInJet[iEta], dBinsEtaInJet[iEta + 1]), iMyColors[(iNBinsEtaInJet - iEta - 1 - 1) % (iNBinsEtaInJet / 2) % iNMyColors], iMyMarkersFull[(iNBinsEtaInJet - iEta - 1) % (iNBinsEtaInJet / 2) % iNMyMarkersFull]);
              }
              else
              {
                grEffPtEtaInJet = MakeGraphErrors(hisEffPtEtaInJet, Form("|#it{#eta}^{gen.}| < %.2f", -dBinsEtaInJet[iEta]), iMyColors[(iEta - 1) % iNMyColors], iMyMarkersFull[(iEta - 1) % iNMyMarkersFull]);
                grEffPtEtaRatio = MakeGraphErrors(hisEffPtEtaRatio, Form("|#it{#eta}^{gen.}| < %.2f", -dBinsEtaInJet[iEta]), iMyColors[(iEta - 1) % iNMyColors], iMyMarkersFull[(iEta - 1) % iNMyMarkersFull]);
                grGenPtEtaRatio = MakeGraphErrors(hisGenPtEtaRatio, Form("|#it{#eta}^{gen.}| < %.2f", -dBinsEtaInJet[iEta]), iMyColors[(iEta - 1) % iNMyColors], iMyMarkersFull[(iEta - 1) % iNMyMarkersFull]);
              }
              mgrEffPtEtaInJetPos->Add(grEffPtEtaInJet);
              mgrEffPtEtaRatioPos->Add(grEffPtEtaRatio);
              mgrGenPtEtaRatioPos->Add(grGenPtEtaRatio);
            }
          }
          canEffPtEtaInJet->cd();
          canEffPtEtaInJet->SetLeftMargin(0.15);
          if(bOpenPtBins)
            mgrEffPtEtaInJet->SetTitle(Form("%s reconstruction efficiency, in JC, c. %s, pT jet > %.0f GeV/c;#it{p}_{T}^{gen.} (GeV/#it{c});efficiency", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet]));
          else
            mgrEffPtEtaInJet->SetTitle(Form("%s reconstruction efficiency, in JC, c. %s, pT jet: %.0f-%.0f GeV/c;#it{p}_{T}^{gen.} (GeV/#it{c});efficiency", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet], dBinsPtJet[iJet + 1]));
          mgrEffPtEtaInJet->SetMinimum(0);
          mgrEffPtEtaInJet->SetMaximum(0.6);
          mgrEffPtEtaInJet->Draw("AP0");
          mgrEffPtEtaInJet->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
          mgrEffPtEtaInJet->GetYaxis()->SetTitleOffset(2);
          legend = canEffPtEtaInJet->BuildLegend(0.2, 0.6, 0.55, 0.85);
          SetLegend(legend);
          canEffPtEtaInJet->SaveAs(Form("canEffPtEtaInJets%s_C%d_J%d_Neg.%s", V0Name[V0Type].Data(), iCent, iJet, sImageSuf.Data()));
          delete canEffPtEtaInJet;
          delete mgrEffPtEtaInJet;

          if(!bEffEtaWindow)
          {
            canEffPtEtaInJetPos->cd();
            canEffPtEtaInJetPos->SetLeftMargin(0.15);
            if(bOpenPtBins)
              mgrEffPtEtaInJetPos->SetTitle(Form("%s reconstruction efficiency, in JC, c. %s, pT jet > %.0f GeV/c;#it{p}_{T}^{gen.} (GeV/#it{c});efficiency", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet]));
            else
              mgrEffPtEtaInJetPos->SetTitle(Form("%s reconstruction efficiency, in JC, c. %s, pT jet: %.0f-%.0f GeV/c;#it{p}_{T}^{gen.} (GeV/#it{c});efficiency", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet], dBinsPtJet[iJet + 1]));
            mgrEffPtEtaInJetPos->SetMinimum(0);
            mgrEffPtEtaInJetPos->SetMaximum(0.6);
            mgrEffPtEtaInJetPos->Draw("AP0");
            mgrEffPtEtaInJetPos->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
            mgrEffPtEtaInJetPos->GetYaxis()->SetTitleOffset(2);
            legend = canEffPtEtaInJetPos->BuildLegend(0.2, 0.6, 0.55, 0.85);
            SetLegend(legend);
            canEffPtEtaInJetPos->SaveAs(Form("canEffPtEtaInJets%s_C%d_J%d_Pos.%s", V0Name[V0Type].Data(), iCent, iJet, sImageSuf.Data()));
          }
          delete canEffPtEtaInJetPos;
          delete mgrEffPtEtaInJetPos;

          canEffPtEtaRatio->cd();
          canEffPtEtaRatio->SetLeftMargin(0.15);
          if(bOpenPtBins)
            mgrEffPtEtaRatio->SetTitle(Form("%s reconstruction efficiency, in JC/inclusive, c. %s, pT jet > %.0f GeV/c;#it{p}_{T}^{gen.} (GeV/#it{c});ratio", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet]));
          else
            mgrEffPtEtaRatio->SetTitle(Form("%s reconstruction efficiency, in JC/inclusive, c. %s, pT jet: %.0f-%.0f GeV/c;#it{p}_{T}^{gen.} (GeV/#it{c});ratio", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet], dBinsPtJet[iJet + 1]));
          mgrEffPtEtaRatio->SetMinimum(0.8);
          mgrEffPtEtaRatio->SetMaximum(1.3);
          mgrEffPtEtaRatio->Draw("AP0");
          mgrEffPtEtaRatio->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
          mgrEffPtEtaRatio->GetYaxis()->SetTitleOffset(2);
          legend = canEffPtEtaRatio->BuildLegend(0.5, 0.6, 0.85, 0.85);
          SetLegend(legend);
          canEffPtEtaRatio->SaveAs(Form("canEffPtEtaRatio%s_C%d_J%d_Neg.%s", V0Name[V0Type].Data(), iCent, iJet, sImageSuf.Data()));
          delete canEffPtEtaRatio;
          delete mgrEffPtEtaRatio;

          if(!bEffEtaWindow)
          {
            canEffPtEtaRatioPos->cd();
            canEffPtEtaRatioPos->SetLeftMargin(0.15);
            if(bOpenPtBins)
              mgrEffPtEtaRatioPos->SetTitle(Form("%s reconstruction efficiency, in JC/inclusive, c. %s, pT jet > %.0f GeV/c;#it{p}_{T}^{gen.} (GeV/#it{c});ratio", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet]));
            else
              mgrEffPtEtaRatioPos->SetTitle(Form("%s reconstruction efficiency, in JC/inclusive, c. %s, pT jet: %.0f-%.0f GeV/c;#it{p}_{T}^{gen.} (GeV/#it{c});ratio", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet], dBinsPtJet[iJet + 1]));
            mgrEffPtEtaRatioPos->SetMinimum(0.8);
            mgrEffPtEtaRatioPos->SetMaximum(1.3);
            mgrEffPtEtaRatioPos->Draw("AP0");
            mgrEffPtEtaRatioPos->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
            mgrEffPtEtaRatioPos->GetYaxis()->SetTitleOffset(2);
            legend = canEffPtEtaRatioPos->BuildLegend(0.5, 0.6, 0.85, 0.85);
            SetLegend(legend);
            canEffPtEtaRatioPos->SaveAs(Form("canEffPtEtaRatio%s_C%d_J%d_Pos.%s", V0Name[V0Type].Data(), iCent, iJet, sImageSuf.Data()));
          }
          delete canEffPtEtaRatioPos;
          delete mgrEffPtEtaRatioPos;

          canGenPtEtaRatio->cd();
          canGenPtEtaRatio->SetLeftMargin(0.15);
          if(bOpenPtBins)
            mgrGenPtEtaRatio->SetTitle(Form("%s generated, in JC/inclusive, c. %s, pT jet > %.0f GeV/c;#it{p}_{T}^{gen.} (GeV/#it{c});ratio", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet]));
          else
            mgrGenPtEtaRatio->SetTitle(Form("%s generated, in JC/inclusive, c. %s, pT jet: %.0f-%.0f GeV/c;#it{p}_{T}^{gen.} (GeV/#it{c});ratio", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet], dBinsPtJet[iJet + 1]));
//                  mgrGenPtEtaRatio->SetMinimum(0.8);
//                  mgrGenPtEtaRatio->SetMaximum(1.3);
          mgrGenPtEtaRatio->Draw("AP0");
          mgrGenPtEtaRatio->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
          mgrGenPtEtaRatio->GetYaxis()->SetTitleOffset(2);
          legend = canGenPtEtaRatio->BuildLegend(0.5, 0.6, 0.85, 0.85);
          SetLegend(legend);
          canGenPtEtaRatio->SaveAs(Form("canGenPtEtaRatio%s_C%d_J%d_Neg.%s", V0Name[V0Type].Data(), iCent, iJet, sImageSuf.Data()));
          delete canGenPtEtaRatio;
          delete mgrGenPtEtaRatio;

          if(!bEffEtaWindow)
          {
            canGenPtEtaRatioPos->cd();
            canGenPtEtaRatioPos->SetLeftMargin(0.15);
            if(bOpenPtBins)
              mgrGenPtEtaRatioPos->SetTitle(Form("%s generated, in JC/inclusive, c. %s, pT jet > %.0f GeV/c;#it{p}_{T}^{gen.} (GeV/#it{c});ratio", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet]));
            else
              mgrGenPtEtaRatioPos->SetTitle(Form("%s generated, in JC/inclusive, c. %s, pT jet: %.0f-%.0f GeV/c;#it{p}_{T}^{gen.} (GeV/#it{c});ratio", V0Symbol[V0Type].Data(), GetCentBinLabel(iCent).Data(), dBinsPtJet[iJet], dBinsPtJet[iJet + 1]));
//                      mgrGenPtEtaRatioPos->SetMinimum(0.8);
//                      mgrGenPtEtaRatioPos->SetMaximum(1.3);
            mgrGenPtEtaRatioPos->Draw("AP0");
            mgrGenPtEtaRatioPos->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
            mgrGenPtEtaRatioPos->GetYaxis()->SetTitleOffset(2);
            legend = canGenPtEtaRatioPos->BuildLegend(0.5, 0.6, 0.85, 0.85);
            SetLegend(legend);
            canGenPtEtaRatioPos->SaveAs(Form("canGenPtEtaRatio%s_C%d_J%d_Pos.%s", V0Name[V0Type].Data(), iCent, iJet, sImageSuf.Data()));
          }
          delete canGenPtEtaRatioPos;
          delete mgrGenPtEtaRatioPos;
          printf("iEffPtEtaInJets: C%d J%d: End\n", iCent, iJet);
        }
      }
    }
    if(iEffPtInclusive)
    {
      canEffPtIncl->cd();
      canEffPtIncl->SetLeftMargin(0.15);
      mgrEffPtIncl->SetTitle(Form("%s reconstruction efficiency, inclusive;#it{p}_{T}^{gen.} (GeV/#it{c});efficiency", V0Symbol[V0Type].Data()));
      mgrEffPtIncl->SetMinimum(0);
      mgrEffPtIncl->SetMaximum(0.6);
      mgrEffPtIncl->Draw("AP0");
      mgrEffPtIncl->GetXaxis()->SetLimits(fPtAllXMin, fPtAllXMax);
      mgrEffPtIncl->GetYaxis()->SetTitleOffset(2);
      legend = canEffPtIncl->BuildLegend(0.2, 0.6, 0.35, 0.85);
      SetLegend(legend);
//          labelSystem=labelCollision->DrawLatex(fPtAllXMax-(fPtAllXMax-fPtAllXMin)/5.,0.5,sLabelCollisionText.Data());
      canEffPtIncl->SaveAs(Form("canEffPtInclusive%s.%s", V0Name[V0Type].Data(), sImageSuf.Data()));
    }
    delete canEffPtIncl;
    delete mgrEffPtIncl;

    if(iMassResol)
    {
      TCanvas* canMassResolPt = new TCanvas("canResol", "", iCanWidth, iCanHeight);
      TMultiGraph* mgrMassResolPt = new TMultiGraph();
      TGraphErrors* grMassResolPt;
      TCanvas* canMassResolPtPeaks[iNCentBins];
      TMultiGraph* mgrMassResolPtPeaks;
      TGraphErrors* grMassResolPtPeaks;
      TString hisMassResolName = Form("fh2V0%sMCResolMPt_%%d", V0Name[V0Type].Data());
      TH1D* hisMassResol1D;
      TH1D* hisMassResolPt;
      for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
      {
        if(!eventsCent[iCent])
          continue;
        TH2D* hisMassResol2D = GetHistogram2D(listMC, Form(hisMassResolName.Data(), iCent));
        if(!hisMassResol2D)
          return;
        canMassResolPtPeaks[iCent] = new TCanvas(Form("canMRPeaks%d", iCent), "", iCanWidth, iCanHeight);
        hisMassResolPt = new TH1D(Form("hisResolPt_%d", iCent), Form("resolution;pt"), 100, 0, 10);
        hisMassResolPt = (TH1D*)hisMassResolPt->Rebin(iNBinsPtV0AllLF, "", dBinsPtV0AllLF);
        mgrMassResolPtPeaks = new TMultiGraph();

        for(Int_t iPt = iBinPtInclFirst - 1; iPt <= iBinPtInclLast - 1; iPt++)
        {
          Int_t iPtFirstBin = hisMassResol2D->GetYaxis()->FindBin(dBinsPtV0AllLF[iPt] + dEpsilon);
          Int_t iPtLastBin = hisMassResol2D->GetYaxis()->FindBin(dBinsPtV0AllLF[iPt + 1] - dEpsilon);
          hisMassResol1D = (TH1D*)hisMassResol2D->ProjectionX(Form("%d_%d", iCent, iPt), iPtFirstBin, iPtLastBin, "e");
          grMassResolPtPeaks = new TGraphErrors(hisMassResol1D);
          grMassResolPtPeaks->SetTitle(Form("%.1f-%.1f", dBinsPtV0AllLF[iPt], dBinsPtV0AllLF[iPt + 1]));
          grMassResolPtPeaks->SetLineColor(iMyColors[iPt % iNMyColors]);
//                  grMassResolPtPeaks->SetMarkerColor(iMyColors[iPt%7]);
//                  grMassResolPtPeaks->SetMarkerStyle(iMyMarkersFull[iPt/7]);
          mgrMassResolPtPeaks->Add(grMassResolPtPeaks);
          Float_t fEnt = hisMassResol1D->GetEntries();
          if(fEnt > 10)
          {
            hisMassResol1D->Fit("gaus");
            TF1* funFitResultResol = hisMassResol1D->GetFunction("gaus");
            if(funFitResultResol)
            {
              hisMassResolPt->SetBinContent(iPt + 1, funFitResultResol->GetParameter(2));
              hisMassResolPt->SetBinError(iPt + 1, funFitResultResol->GetParError(2));
            }
          }
          else
          {
            printf("Too few entries to fit: %d\n", int(fEnt));
          }
        }
        grMassResolPt = new TGraphErrors(hisMassResolPt);
        grMassResolPt->SetTitle(GetCentBinLabel(iCent).Data());
        grMassResolPt->SetLineColor(iMyColors[iCent]);
        grMassResolPt->SetMarkerColor(iMyColors[iCent]);
        grMassResolPt->SetMarkerStyle(iMyMarkersFull[iCent]);
        mgrMassResolPt->Add(grMassResolPt);
        canMassResolPtPeaks[iCent]->cd();
        mgrMassResolPtPeaks->SetTitle(Form("Mass resolution, c. %s", GetCentBinLabel(iCent).Data()));
        mgrMassResolPtPeaks->Draw("AP0");
        legend = canMassResolPtPeaks[iCent]->BuildLegend(0.15, 0.2, 0.3, 0.85);
        legend->SetBorderSize(0);
        legend->SetFillColor(kWhite);
        canMassResolPtPeaks[iCent]->SetLogy();
        canMassResolPtPeaks[iCent]->SaveAs(Form("canMassResolPt%s_%d.%s", V0Name[V0Type].Data(), iCent, sImageSuf.Data()));
        delete canMassResolPtPeaks[iCent];
      }
      canMassResolPt->cd();
      mgrMassResolPt->SetTitle(Form("Mass resolution;pT;sigma"));
      if(V0Type == 0)
      {
        mgrMassResolPt->SetMinimum(0.004);
        mgrMassResolPt->SetMaximum(0.009);
//                      mgrMassResolPt->SetMaximum(0.019);
      }
      else
      {
        mgrMassResolPt->SetMinimum(0.0013);
        mgrMassResolPt->SetMaximum(0.0033);
      }
      mgrMassResolPt->Draw("AP0");
      legend = canMassResolPt->BuildLegend(0.15, 0.6, 0.3, 0.85);
      legend->SetBorderSize(0);
      legend->SetFillColor(kWhite);
      canMassResolPt->SaveAs(Form("canMassResolPt%s.%s", V0Name[V0Type].Data(), sImageSuf.Data()));
      delete canMassResolPt;
    }
  }
  /*
    if (iCuts)
      {
        if (iTestCut)
          {
            TCanvas* canTestCanvas[2][iNPtBins];
            TGraphErrors* grTest;
            TMultiGraph* mgrTest[2][iNPtBins];
            TString hisTestHistName = Form("fh2CutTest%s_%%d_Pt%%d",V0Name[V0Type].Data());
            TString testVariable = "#tau";
            TH2D* hisTestHist;
            for (Int_t i=0; i<1; i++)
              for (Int_t iPt=0; iPt<iNPtBins; iPt++)
                {
                  canTestCanvas[i][iPt] = new TCanvas(Form("canTestCanvas_%d_%d",i,iPt),"",iCanWidth,iCanHeight);
                  mgrTest[i][iPt] = new TMultiGraph();
                  hisTestHist = GetHistogram2D(listCuts,Form(hisTestHistName.Data(),i,iPt));
                  if (!hisTestHist)
                    return;
  //          CheckHistogram(hisTestHist);
  //                    hisTestHist->SetTitle(Form("%s: %s distribution vs mass, cent: %s",V0Name[V0Type].Data(),testVariable.Data(),GetCentBinLabel(iPt).Data()));
                  hisTestHist->GetYaxis()->SetTitle(testVariable.Data());
  //          hisTestHist->GetYaxis()->SetTitleOffset(1.7);
                  hisTestHist->GetXaxis()->SetTitle(V0LabelM[V0Type].Data());
                  for (Int_t iCut = 1; iCut<=10; iCut++)
                    {
                      TH1D* hisProjCut = (TH1D*)hisTestHist->ProjectionX(Form("p%d%d",i,iPt),1,iCut*10,"e");
                      grTest = new TGraphErrors(hisProjCut);
                      grTest->SetLineColor(iMyColors[(iCut-1)%iNMyColors]);
                      grTest->SetTitle(Form("%d #tau",iCut));
                      mgrTest[i][iPt]->Add(grTest);
                    }
                  canTestCanvas[i][iPt]->cd();
                  mgrTest[i][iPt]->Draw("AP0");
                  canTestCanvas[i][iPt]->SetLogz(1);
                  canTestCanvas[i][iPt]->SetRightMargin(0.15);
                  legend = canTestCanvas[i][iPt]->BuildLegend(0.7,0.6,0.85,0.85);
                  SetLegend(legend);
                }
            for (Int_t i=0; i<1; i++)
              for (Int_t iPt=0; iPt<iNPtBins; iPt++)
                {
                  canTestCanvas[i][iPt]->SaveAs(Form("canTestCut%s_%d-Pt%d.%s",V0Name[V0Type].Data(),i,iPt,sImageSuf.Data()));
                }
            //====================================================

            TCanvas* canTuneCutAll[2][iNPtBins];
            TCanvas* canTuneCutBg[2][iNPtBins];
            TCanvas* canTuneCutSignal[2][iNPtBins];
            TCanvas* canTuneCutSB[2][iNPtBins];
            TCanvas* canTuneCutCand[2][iNPtBins];

            TGraphErrors* grTuneCutAll;
            TGraphErrors* grTuneCutBg;
            TGraphErrors* grTuneCutSignal;

            for (Int_t i=0; i<2; i++)
              for (Int_t iPt=0; iPt<iNPtBins; iPt++)
                {
                  canTuneCutAll[i][iPt] = new TCanvas(Form("can%d%d",i,iPt),"",iCanWidth,iCanHeight);
                  canTuneCutBg[i][iPt] = new TCanvas(Form("canBg%d%d",i,iPt),"",iCanWidth,iCanHeight);
                  canTuneCutSignal[i][iPt] = new TCanvas(Form("canSig%d%d",i,iPt),"",iCanWidth,iCanHeight);
                  canTuneCutSB[i][iPt] = new TCanvas(Form("canSB%d%d",i,iPt),"",iCanWidth,iCanHeight);
                  canTuneCutCand[i][iPt] = new TCanvas(Form("canCand%d%d",i,iPt),"",iCanWidth,iCanHeight);

                  TH1D* hisTuneCut;
                  TH1D* hisTuneCutSB = new TH1D(Form("hisSB%d%d",i,iPt),Form("%s signal purity, %s",V0Symbol[V0Type].Data(),GetPtBinLabel(iPt).Data()),10,0,10);
                  TH1D* hisTuneCutCand = new TH1D(Form("hisCand%d%d",i,iPt),Form("%s signal candidates, %s",V0Symbol[V0Type].Data(),GetPtBinLabel(iPt).Data()),10,0,10);
                  TMultiGraph* mgrTuneCutAll = new TMultiGraph();
                  TMultiGraph* mgrTuneCutBg = new TMultiGraph();
                  TMultiGraph* mgrTuneCutSignal = new TMultiGraph();
                  TString hisTuneCutName = Form("fh2CutTest%s_%%d_Pt%%d",V0Name[V0Type].Data());

                  hisTestHist = GetHistogram2D(listCuts,Form(hisTestHistName.Data(),i,iPt));
                  if (!hisTestHist)
                    return;
                  hisTestHist->GetYaxis()->SetTitle(testVariable.Data());
                  hisTestHist->GetXaxis()->SetTitle(V0LabelM[V0Type].Data());

                  for (Int_t iCut=1; iCut<=10; iCut++)
                    {
                      hisTuneCut = (TH1D*)hisTestHist->ProjectionX(Form("p%d%d",i,iPt),1,iCut*10,"e");
                      if (!hisTuneCut)
                        return;
                      grTuneCutAll = new TGraphErrors(hisTuneCut);
                      grTuneCutAll->SetLineColor(iMyColors[(iCut-1)%iNMyColors]);
                      grTuneCutAll->SetMarkerColor(iMyColors[(iCut-1)%iNMyColors]);
                      grTuneCutAll->SetTitle(Form("%d#tau",iCut));
                      mgrTuneCutAll->Add(grTuneCutAll);

                      // copy of the histogram with excluded peak
                      TH1D* hisTuneCutBg = new TH1D("hisTuneCutBg","",hisTuneCut->GetNbinsX(),hisTuneCut->GetBinLowEdge(1),hisTuneCut->GetBinLowEdge(hisTuneCut->GetNbinsX()+1));
                      TH1D* hisTuneCutSignal = new TH1D("hisTuneCutSignal","",hisTuneCut->GetNbinsX(),hisTuneCut->GetBinLowEdge(1),hisTuneCut->GetBinLowEdge(hisTuneCut->GetNbinsX()+1));
                      for (Int_t j = 1; j<=hisTuneCutBg->GetNbinsX(); j++)
                        {
                          if ((hisTuneCutBg->GetBinLowEdge(j) > fRangeSignal[V0Type][0]) && (hisTuneCutBg->GetBinLowEdge(j) < fRangeSignal[V0Type][1]))
                            continue;
                          hisTuneCutBg->SetBinContent(j,hisTuneCut->GetBinContent(j));
                          hisTuneCutBg->SetBinError(j,hisTuneCut->GetBinError(j));
                        }
                      // fit background
                      hisTuneCutBg->Fit(funBg,"R");
                      grTuneCutBg = new TGraphErrors(hisTuneCutBg);
                      grTuneCutBg->SetLineColor(iMyColors[(iCut-1)%iNMyColors]);
                      grTuneCutBg->SetMarkerColor(iMyColors[(iCut-1)%iNMyColors]);
                      grTuneCutBg->SetTitle(Form("%d#tau",iCut));
                      mgrTuneCutBg->Add(grTuneCutBg);
                      // get and plot the fit function
  //                        TF1* funFitResultBg = hisTuneCutBg->GetFunction(funBg);
                      if (funBg)
                        {
                          TGraph* grFitResultBg = new TGraph(funBg);
                          grFitResultBg->SetLineColor(iMyColors[(iCut-1)%7]);
                          grFitResultBg->SetLineStyle(2);
                          mgrTuneCutBg->Add(grFitResultBg);
                          mgrTuneCutAll->Add(grFitResultBg);
                        }

                      Double_t fIntegralAllBC = 0; // sum of bin content under the peak
                      Double_t fIntegralSignalBC = 0; // sum of bin content in the peak only
                      Double_t fBinDiff; // signal height in a bin
                      Double_t fIntegralAllBCError2 = 0; // square error of the sum under peak
                      for (Int_t j = 1; j<=hisTuneCut->GetNbinsX(); j++)
                        {
                          if ((hisTuneCut->GetBinLowEdge(j) < fRangeSignal[V0Type][0]) || (hisTuneCut->GetBinLowEdge(j) > fRangeSignal[V0Type][1]))
                            continue;
                          fIntegralAllBC+=hisTuneCut->GetBinContent(j);
                          fIntegralAllBCError2+=(hisTuneCut->GetBinError(j))*(hisTuneCut->GetBinError(j));
                          if (funBg)
                            fBinDiff=hisTuneCut->GetBinContent(j)-funBg->Eval(hisTuneCut->GetBinCenter(j));
                          else
                            fBinDiff=hisTuneCut->GetBinContent(j);
                          fIntegralSignalBC+=fBinDiff;
                          hisTuneCutSignal->SetBinContent(j,fBinDiff);
                          hisTuneCutSignal->SetBinError(j,hisTuneCut->GetBinError(j));
                        }
                      Double_t fSB = fIntegralSignalBC/fIntegralAllBC;
                      Double_t fSBError = (fIntegralAllBC-fIntegralSignalBC)/(fIntegralAllBC*fIntegralAllBC)*sqrt(fIntegralAllBCError2);
                      printf("Signal purity for %s: %f +- %f\n",Form("%d#tau",iCut),fSB,fSBError);
                      hisTuneCutSB->SetBinContent(iCut,fSB);
                      hisTuneCutSB->SetBinError(iCut,fSBError);
                      hisTuneCutSB->GetXaxis()->SetBinLabel(iCut,Form("%d#tau",iCut));
                      hisTuneCutCand->SetBinContent(iCut,fIntegralSignalBC);
                      hisTuneCutCand->SetBinError(iCut,sqrt(fIntegralAllBCError2));
                      hisTuneCutCand->GetXaxis()->SetBinLabel(iCut,Form("%d#tau",iCut));
                      // plot the extracted signal
                      grTuneCutSignal = new TGraphErrors(hisTuneCutSignal);
                      grTuneCutSignal->SetTitle(Form("%d#tau",iCut));
                      grTuneCutSignal->SetLineColor(iMyColors[(iCut-1)%iNMyColors]);
                      grTuneCutSignal->SetMarkerColor(iMyColors[(iCut-1)%iNMyColors]);
                      mgrTuneCutSignal->Add(grTuneCutSignal);
                    }
                  canTuneCutAll[i][iPt]->cd();
  //                    canTuneCutAll[i][iPt]->SetLogy();
                  mgrTuneCutAll->SetTitle(Form("%s signal + background, %s",V0Symbol[V0Type].Data(),GetPtBinLabel(iPt).Data()));
                  mgrTuneCutAll->Draw("apl");
                  legend = canTuneCutAll[i][iPt]->BuildLegend(0.7,0.2,0.85,0.85);
                  legend->SetBorderSize(0);
                  legend->SetFillColor(kWhite);
                  canTuneCutAll[i][iPt]->SaveAs(Form("canTuneCutAll%s-%d-%d.%s",V0Name[V0Type].Data(),i,iPt,sImageSuf.Data()));
                  canTuneCutBg[i][iPt]->cd();
  //                    canTuneCutBg[i][iPt]->SetLogy();
                  mgrTuneCutBg->SetTitle(Form("%s background, %s",V0Symbol[V0Type].Data(),GetPtBinLabel(iPt).Data()));
                  mgrTuneCutBg->Draw("apl");
                  legend = canTuneCutBg[i][iPt]->BuildLegend(0.7,0.2,0.85,0.85);
                  legend->SetBorderSize(0);
                  legend->SetFillColor(kWhite);
                  canTuneCutBg[i][iPt]->SaveAs(Form("canTuneCutBg%s-%d-%d.%s",V0Name[V0Type].Data(),i,iPt,sImageSuf.Data()));
                  canTuneCutSignal[i][iPt]->cd();
  //                    canTuneCutSignal[i][iPt]->SetLogy();
                  mgrTuneCutSignal->SetMinimum(1);
                  mgrTuneCutSignal->SetTitle(Form("%s extracted signal, %s",V0Symbol[V0Type].Data(),GetPtBinLabel(iPt).Data()));
                  mgrTuneCutSignal->Draw("apl");
                  legend = canTuneCutSignal[i][iPt]->BuildLegend(0.7,0.2,0.85,0.85);
                  legend->SetBorderSize(0);
                  legend->SetFillColor(kWhite);
                  canTuneCutSignal[i][iPt]->SaveAs(Form("canTuneCutSignal%s-%d-%d.%s",V0Name[V0Type].Data(),i,iPt,sImageSuf.Data()));
                  canTuneCutSB[i][iPt]->cd();
                  hisTuneCutSB->Draw();
                  canTuneCutSB[i][iPt]->SaveAs(Form("canTuneCutSB%s-%d-%d.%s",V0Name[V0Type].Data(),i,iPt,sImageSuf.Data()));
                  canTuneCutCand[i][iPt]->cd();
                  hisTuneCutCand->Draw();
                  canTuneCutCand[i][iPt]->SaveAs(Form("canTuneCutCand%s-%d-%d.%s",V0Name[V0Type].Data(),i,iPt,sImageSuf.Data()));
                }
          }
      }
  */
  if(iJetSpectrum)
  {
    Double_t dEtaJetMax = 0.8;

    TCanvas* canPtJet = new TCanvas("canPtJet", "", iCanWidth, iCanHeight);
    TMultiGraph* mgrPtJet = new TMultiGraph();
    TGraphErrors* grPtJet;
    TString hisPtJetName = "fh1PtJet_%d";
    TH1D* hisPtJet;

    TCanvas* canPhiJet = new TCanvas("canPhiJet", "", iCanWidth, iCanHeight);
    TMultiGraph* mgrPhiJet = new TMultiGraph();
    TGraphErrors* grPhiJet;
    TString hisPhiJetName = "fh1PhiJet_%d";
    TH1D* hisPhiJet;

    TCanvas* canEtaJet = new TCanvas("canEtaJet", "", iCanWidth, iCanHeight);
    TMultiGraph* mgrEtaJet = new TMultiGraph();
    TGraphErrors* grEtaJet;
    TString hisEtaJetName = "fh1EtaJet_%d";
    TH1D* hisEtaJet;

    TCanvas* canEtaPtJet;
    TMultiGraph* mgrEtaPtJet;
    TGraphErrors* grEtaPtJet;
    TString hisEtaPtJetName = "fh2EtaPtJet_%d";
    TH2D* hisEtaPtJet;

    TCanvas* canNJet = new TCanvas("canNJet", "", iCanWidth, iCanHeight);
    TMultiGraph* mgrNJet = new TMultiGraph();
    TGraphErrors* grNJet;
    TString hisNJetName = "fh1NJetPerEvent_%d";
    TH1D* hisNJet;

    TCanvas* canPtJetPtTrack;
    TCanvas* canPtTrack;
    TMultiGraph* mgrPtTrack;
    TGraphErrors* grPtTrack;
    TString hisPtJetPtTrackName = "fh2PtJetPtTrackLeading_%d";
    TH2D* hisPtJetPtTrack;
    TH1D* hisPtTrack;

    TCanvas* canPtJetPtTrigger;
    TCanvas* canPtTrigger;
    TMultiGraph* mgrPtTrigger;
    TGraphErrors* grPtTrigger;
    TString hisPtJetPtTriggerName = "fh2PtJetPtTrigger_%d";
    TH2D* hisPtJetPtTrigger;
    TH1D* hisPtTrigger;
    TString hisPtTriggerName = "fh1PtTrigger_%d";

    printf("Statistics: Jets (pt bins):");
    for(Int_t iJet = iJetMin; iJet <= iJetMax; iJet++)
      printf(" %.1f", dBinsPtJet[iJet]);
    printf(", open bins: %d\n", bOpenPtBins);

    for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
    {
      if(!eventsCent[iCent])
        continue;
      hisPtJet = GetHistogram1D(listStd, Form(hisPtJetName.Data(), iCent));
      hisPhiJet = GetHistogram1D(listStd, Form(hisPhiJetName.Data(), iCent));
      hisEtaJet = GetHistogram1D(listStd, Form(hisEtaJetName.Data(), iCent));
      hisEtaPtJet = GetHistogram2D(listStd, Form(hisEtaPtJetName.Data(), iCent));
      hisNJet = GetHistogram1D(listStd, Form(hisNJetName.Data(), iCent));
      if(!hisPtJet || !hisPhiJet || !hisEtaJet || !hisEtaPtJet || !hisNJet)
        return;
      if(!hisEtaJet->GetEntries())
      {
        printf("Histograms empty\n");
        continue;
      }

      if(iPtCorrel)
      {
        hisPtJetPtTrack = GetHistogram2D(listStd, Form(hisPtJetPtTrackName.Data(), iCent));
        hisPtJetPtTrigger = GetHistogram2D(listStd, Form(hisPtJetPtTriggerName.Data(), iCent));
        hisPtTrigger = GetHistogram1D(listStd, Form(hisPtTriggerName.Data(), iCent));
        if(!hisPtJetPtTrack || !hisPtJetPtTrigger || !hisPtTrigger)
          return;
      }

      CheckHistogram(hisPtJet);
      if(!bYieldsOnly)
        hisPtJet->Scale(1. / eventsCent[iCent]);

      Double_t dArrayJetNorm[iNBinsPtJet];
      for(Int_t i = 0; i < iNBinsPtJet; i++)
      {
        Int_t iBinFirst = hisPtJet->GetXaxis()->FindBin(dBinsPtJet[i] + dEpsilon);
        Int_t iBinLast;
        if(bOpenPtBins)
          iBinLast = hisPtJet->GetNbinsX() + 1;
        else
          iBinLast = hisPtJet->GetXaxis()->FindBin(dBinsPtJet[i + 1] - dEpsilon);
        dArrayJetNorm[i] = hisPtJet->Integral(iBinFirst, iBinLast);
      }
      printf("Statistics: Jets per event (Cent %d):", iCent);
      for(Int_t iJet = iJetMin; iJet <= iJetMax; iJet++)
        printf(" J%d %g,", iJet, dArrayJetNorm[iJet]);
      printf("\n");

      if(!bYieldsOnly)
        hisPtJet->Scale(1., "width");

      dirOutSpectra->cd();
      hisPtJet->Write();
      grPtJet = MakeGraphErrors(hisPtJet, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersFull[iCent]);
      mgrPtJet->Add(grPtJet);

      hisPhiJet->Scale(1. / (hisPhiJet->Integral() == 0 ? 1 : hisPhiJet->Integral()));
//              hisPhiJet->Scale(1./eventsCent[iCent]);
      grPhiJet = MakeGraphErrors(hisPhiJet, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersFull[iCent]);
      mgrPhiJet->Add(grPhiJet);
//          hisEtaJet->Scale(1./hisEtaJet->Integral());
      hisEtaJet->Scale(1. / eventsCent[iCent], "width");
      grEtaJet = MakeGraphErrors(hisEtaJet, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersFull[iCent]);
      mgrEtaJet->Add(grEtaJet);
      CheckHistogram(hisNJet);
      hisNJet->Scale(1. / eventsCent[iCent]);
      grNJet = MakeGraphErrors(hisNJet, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersFull[iCent]);
      mgrNJet->Add(grNJet);
      printf("Statistics: Jets per event (Cent %d):", iCent);
      Double_t dJetEvents = 1 - hisNJet->GetBinContent(1); // fraction of jet events = 1 - fraction of no-jet events
      Double_t dMultiJetEvents = 0;
      for(Int_t iNJet = 0; iNJet <= 10; iNJet++)
      {
        printf(" %d jets %g,", iNJet, hisNJet->GetBinContent(iNJet + 1));
        if(iNJet > 1)
          dMultiJetEvents += hisNJet->GetBinContent(iNJet + 1);
      }
      printf("\nStatistics: Jets per event (Cent %d): No jet events: %g, jet events %g, 1-jet events: %g, multi-jet events %g (%g %% of jet events), total %g\n", iCent, hisNJet->GetBinContent(1), dJetEvents, hisNJet->GetBinContent(2), dMultiJetEvents, 100 * dMultiJetEvents / dJetEvents, hisNJet->Integral());

      canEtaPtJet = new TCanvas(Form("canEtaPtJet_C%d", iCent), "", iCanWidth, iCanHeight);
      mgrEtaPtJet = new TMultiGraph();
      for(Int_t iJet = iJetMin - 1; iJet <= iJetMax; iJet++)
      {
        Int_t iBinJetFirst = hisEtaPtJet->GetYaxis()->FindBin(dBinsPtJet[iJet] + dEpsilon);
        Int_t iBinJetLast;
        if(bOpenPtBins)
          iBinJetLast = hisEtaPtJet->GetYaxis()->GetNbins() + 1;
        else
          iBinJetLast = hisEtaPtJet->GetYaxis()->FindBin(dBinsPtJet[iJet + 1] - dEpsilon);
        TH1D* hisEtaJetProj = (TH1D*)hisEtaPtJet->ProjectionX(Form("%s_px-%d", hisEtaPtJet->GetName(), iJet), iBinJetFirst, iBinJetLast, "e");
        if(!hisEtaJetProj)
          return;
        hisEtaJetProj->Scale(1. / (hisEtaJetProj->Integral() == 0 ? 1 : hisEtaJetProj->Integral()), "width");
        TString sLabelPtJet;
        if(bOpenPtBins)
          sLabelPtJet = Form("#it{p}_{T}^{jet} > %g GeV/#it{c}", dBinsPtJet[iJet]);
        else
          sLabelPtJet = Form("#it{p}_{T}^{jet}: %g-%g GeV/#it{c}", dBinsPtJet[iJet], dBinsPtJet[iJet + 1]);
        TGraphErrors* grEtaJetProj = MakeGraphErrors(hisEtaJetProj, sLabelPtJet.Data(), iMyColors[iJet % iNMyColors], iMyMarkersFull[iJet % iNMyMarkersFull]);
        mgrEtaPtJet->Add(grEtaJetProj);
      }
      canEtaPtJet->cd();
      canEtaPtJet->SetLeftMargin(0.15);
      //      mgrEtaPtJet->SetTitle("Charged-jet #it{#eta} spectrum;#it{#eta}^{jet,ch};arb. unit");
      mgrEtaPtJet->SetTitle(Form("Charged-jet #it{#eta} spectrum, c. %s;#it{#eta}^{jet,ch};arb. u.", GetCentBinLabel(iCent).Data()));
      mgrEtaPtJet->SetMinimum(0.);
      //      mgrEtaPtJet->SetMaximum(0.05);
      //      mgrEtaPtJet->SetMaximum(0.3);
      mgrEtaPtJet->Draw("AP0");
      mgrEtaPtJet->GetYaxis()->SetTitleOffset(2);
      mgrEtaPtJet->GetXaxis()->SetLimits(-dEtaJetMax, dEtaJetMax);
      legend = canEtaPtJet->BuildLegend(0.6, 0.15, 0.75, 0.35);
      SetLegend(legend);
      labelSystem = labelCollision->DrawLatex(-0.2, 0.4, sLabelCollisionText.Data());
      canEtaPtJet->SaveAs(Form("canEtaPtJet_C%d.%s", iCent, sImageSuf.Data()));
      delete canEtaPtJet;
      delete mgrEtaPtJet;

      if(iPtCorrel)
      {
        canPtJetPtTrack = new TCanvas(Form("canPtJetPtTrack_C%d", iCent), "", iCanWidth, iCanHeight);
        canPtJetPtTrack->cd();
        hisPtJetPtTrack->GetZaxis()->SetRangeUser(1, 1e5);
        hisPtJetPtTrack->Draw("colz");
        canPtJetPtTrack->SetLogz();
        canPtJetPtTrack->SetRightMargin(0.15);
        canPtJetPtTrack->SaveAs(Form("canPtJetPtTrack_C%d.%s", iCent, sImageSuf.Data()));
        delete canPtJetPtTrack;

        canPtTrack = new TCanvas(Form("canPtTrack_C%d", iCent), "", iCanWidth, iCanHeight);
        mgrPtTrack = new TMultiGraph();
        canPtTrigger = new TCanvas(Form("canPtTrigger_C%d", iCent), "", iCanWidth, iCanHeight);
        mgrPtTrigger = new TMultiGraph();
        for(Int_t iJet = iJetMin; iJet <= iJetMax; iJet++)
        {
          Int_t iBinJetFirst = hisPtJetPtTrack->GetXaxis()->FindBin(dBinsPtJet[iJet] + dEpsilon);
          Int_t iBinJetLast;
          if(bOpenPtBins)
            iBinJetLast = hisPtJetPtTrack->GetXaxis()->GetNbins() + 1;
          else
            iBinJetLast = hisPtJetPtTrack->GetXaxis()->FindBin(dBinsPtJet[iJet + 1] - dEpsilon);
          TH1D* hisPtTrackProj = (TH1D*)hisPtJetPtTrack->ProjectionY(Form("%s_py-%d", hisPtJetPtTrack->GetName(), iJet), iBinJetFirst, iBinJetLast, "e");
          if(!hisPtTrackProj)
            return;
          TString sLabelPtJet;
          if(bOpenPtBins)
            sLabelPtJet = Form("#it{p}_{T}^{jet} > %g GeV/#it{c}", dBinsPtJet[iJet]);
          else
            sLabelPtJet = Form("#it{p}_{T}^{jet}: %g-%g GeV/#it{c}", dBinsPtJet[iJet], dBinsPtJet[iJet + 1]);
          TGraphErrors* grPtTrackProj = MakeGraphErrors(hisPtTrackProj, sLabelPtJet.Data(), iMyColors[iJet % iNMyColors], iMyMarkersEmpty[iJet % iNMyMarkersEmpty]);
          TGraphErrors* grPtTrackProjTrig = MakeGraphErrors(hisPtTrackProj, Form("LT, %s", sLabelPtJet.Data()), iMyColors[iJet % iNMyColors], iMyMarkersEmpty[iJet % iNMyMarkersEmpty]);
          mgrPtTrack->Add(grPtTrackProj);
          mgrPtTrigger->Add(grPtTrackProjTrig);
        }
        canPtTrack->cd();
        mgrPtTrack->SetTitle("Leading tracks in jets;#it{p}_{T}^{leading track} (GeV/#it{c});counts");
        mgrPtTrack->Draw("AP0");
        canPtTrack->SetLogy();
        legend = canPtTrack->BuildLegend(0.6, 0.6, 0.75, 0.85);
        SetLegend(legend);
        labelSystem = labelCollision->DrawLatex(-0.2, 0.4, sLabelCollisionText.Data());
        canPtTrack->SaveAs(Form("canPtTrack_C%d.%s", iCent, sImageSuf.Data()));
        delete canPtTrack;

        canPtJetPtTrigger = new TCanvas(Form("canPtJetPtTrigger_C%d", iCent), "", iCanWidth, iCanHeight);
        canPtJetPtTrigger->cd();
        hisPtJetPtTrigger->GetZaxis()->SetRangeUser(1, 1e5);
        hisPtJetPtTrigger->Draw("colz");
        canPtJetPtTrigger->SetLogz();
        canPtJetPtTrigger->SetRightMargin(0.15);
        canPtJetPtTrigger->SaveAs(Form("canPtJetPtTrigger_C%d.%s", iCent, sImageSuf.Data()));
        delete canPtJetPtTrigger;

        grPtTrigger = MakeGraphErrors(hisPtTrigger, "TT in all events", iMyColors[0], iMyMarkersFull[0]);
        mgrPtTrigger->Add(grPtTrigger);
        for(Int_t iJet = iJetMin; iJet <= iJetMax; iJet++)
        {
          Int_t iBinJetFirst = hisPtJetPtTrigger->GetXaxis()->FindBin(dBinsPtJet[iJet] + dEpsilon);
          Int_t iBinJetLast;
          if(bOpenPtBins)
            iBinJetLast = hisPtJetPtTrigger->GetXaxis()->GetNbins() + 1;
          else
            iBinJetLast = hisPtJetPtTrigger->GetXaxis()->FindBin(dBinsPtJet[iJet + 1] - dEpsilon);
          TH1D* hisPtTriggerProj = (TH1D*)hisPtJetPtTrigger->ProjectionY(Form("%s_py-%d", hisPtJetPtTrigger->GetName(), iJet), iBinJetFirst, iBinJetLast, "e");
          if(!hisPtTriggerProj)
            return;
          TString sLabelPtJet;
          if(bOpenPtBins)
            sLabelPtJet = Form("#it{p}_{T}^{jet} > %g GeV/#it{c}", dBinsPtJet[iJet]);
          else
            sLabelPtJet = Form("#it{p}_{T}^{jet}: %g-%g GeV/#it{c}", dBinsPtJet[iJet], dBinsPtJet[iJet + 1]);
          TGraphErrors* grPtTriggerProj = MakeGraphErrors(hisPtTriggerProj, Form("TT, %s", sLabelPtJet.Data()), iMyColors[iJet % iNMyColors], iMyMarkersFull[iJet % iNMyMarkersFull]);
          mgrPtTrigger->Add(grPtTriggerProj);
        }
        hisPtTrigger = (TH1D*)hisPtJetPtTrigger->ProjectionY(Form("%s_py1", hisPtJetPtTrigger->GetName()), 1, 1, "e"); // bin 1 - triggers in no jet events
        Double_t dNEventNoJetNoTrigger = hisPtTrigger->GetBinContent(1); // number of events where neither trigger nor jet was selected
        hisPtTrigger->SetBinContent(1, 0);
        hisPtTrigger->SetBinError(1, 0);
        grPtTrigger = MakeGraphErrors(hisPtTrigger, "TT in no-jet events", iMyColors[(iJetMax + 1) % iNMyColors], iMyMarkersFull[(iJetMax + 1) % iNMyMarkersFull]);
        mgrPtTrigger->Add(grPtTrigger);

        canPtTrigger->cd();
        mgrPtTrigger->SetTitle("Leading tracks in jets and trigger tracks;#it{p}_{T}^{track} (GeV/#it{c});counts");
        mgrPtTrigger->SetMinimum(1);
        mgrPtTrigger->Draw("AP0");
        canPtTrigger->SetLogy();
        legend = canPtTrigger->BuildLegend(0.6, 0.6, 0.75, 0.85);
        SetLegend(legend);
//      legend->AddEntry(Form("soft events: %g",dNEventNoJetNoTrigger));
        labelSystem = labelCollision->DrawLatex(5, 10, sLabelCollisionText.Data());
        labelSystem->DrawLatex(5, 5e2, Form("soft events: %.0f %%", 100 * dNEventNoJetNoTrigger / eventsCent[iCent]));
        canPtTrigger->SaveAs(Form("canPtTrigger_C%d.%s", iCent, sImageSuf.Data()));
        delete canPtTrigger;
        delete mgrPtTrigger;
        delete mgrPtTrack;
      }
    }
    canPtJet->cd();
    canPtJet->SetLeftMargin(0.2);
    mgrPtJet->SetTitle("Charged-jet #it{p}_{T} spectrum;#it{p}_{T}^{jet,ch} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{jet,ch}}{d#it{p}_{T}^{jet,ch}} (#it{c} GeV^{#minus1})");
    mgrPtJet->SetMaximum(fPtJetYMax);
    mgrPtJet->SetMinimum(fPtJetYMin);
    mgrPtJet->Draw("AP0");
    mgrPtJet->GetXaxis()->SetLimits(fPtJetXMin, fPtJetXMax);
    mgrPtJet->GetYaxis()->SetTitleOffset(2);
    canPtJet->SetLogy();
    legend = canPtJet->BuildLegend(0.7, 0.6, 0.85, 0.85);
    SetLegend(legend);
    labelSystem = labelCollision->DrawLatex(20, 1e-6, sLabelCollisionText.Data());
    canPtJet->SaveAs(Form("canPtJet.%s", sImageSuf.Data()));
    delete canPtJet;
    delete mgrPtJet;

    canPhiJet->cd();
    canPhiJet->SetLeftMargin(0.15);
    mgrPhiJet->SetTitle("Charged-jet #it{#phi} spectrum;#it{#phi}^{jet,ch} (rad);probability");
    mgrPhiJet->SetMinimum(0.);
    mgrPhiJet->Draw("AP0");
    mgrPhiJet->GetYaxis()->SetTitleOffset(2);
    legend = canPhiJet->BuildLegend(0.7, 0.15, 0.85, 0.35);
    SetLegend(legend);
    labelSystem = labelCollision->DrawLatex(3, 0.005, sLabelCollisionText.Data());
    canPhiJet->SaveAs(Form("canPhiJet.%s", sImageSuf.Data()));
    delete canPhiJet;
    delete mgrPhiJet;

    canEtaJet->cd();
    canEtaJet->SetLeftMargin(0.15);
//      mgrEtaJet->SetTitle("Charged-jet #it{#eta} spectrum;#it{#eta}^{jet,ch};arb. unit");
    mgrEtaJet->SetTitle("Charged-jet #it{#eta} spectrum;#it{#eta}^{jet,ch};#frac{1}{#it{N}_{ev.}} #frac{d#it{N}}{d#it{#eta}}");
    mgrEtaJet->SetMinimum(0.);
//      mgrEtaJet->SetMaximum(0.05);
//      mgrEtaJet->SetMaximum(0.3);
    mgrEtaJet->Draw("AP0");
    mgrEtaJet->GetYaxis()->SetTitleOffset(2);
    mgrEtaJet->GetXaxis()->SetLimits(-dEtaJetMax, dEtaJetMax);
    legend = canEtaJet->BuildLegend(0.7, 0.15, 0.85, 0.35);
    SetLegend(legend);
    labelSystem = labelCollision->DrawLatex(0, 0.07, sLabelCollisionText.Data());
    canEtaJet->SaveAs(Form("canEtaJet.%s", sImageSuf.Data()));
    delete canEtaJet;
    delete mgrEtaJet;

    canNJet->cd();
    canNJet->SetLeftMargin(0.15);
    mgrNJet->SetTitle("# selected charged jets per event;# selected ch. jets/event;probability");
    mgrNJet->SetMinimum(1e-7);
    mgrNJet->SetMaximum(2);
    mgrNJet->Draw("AP0");
    mgrNJet->GetYaxis()->SetTitleOffset(2);
    mgrNJet->GetXaxis()->SetLimits(0, 6);
    canNJet->SetLogy();
    legend = canNJet->BuildLegend(0.7, 0.65, 0.85, 0.85);
    SetLegend(legend);
    labelSystem = labelCollision->DrawLatex(2, 1e-6, sLabelCollisionText.Data());
    canNJet->SaveAs(Form("canNJet.%s", sImageSuf.Data()));
    delete canNJet;
    delete mgrNJet;
  }

  printf("DrawResults: Closing files\n");
  fileHisto->Close();
  if(fileEff)
    fileEff->Close();
  fileOutput->Close();
  printf("DrawResults: End\n");

  gSystem->cd(sPwd.Data());

  // stop counting and print out
  timer.Stop();
  timer.Print();
}

Double_t GetMeanFeedDownLambdaLF(Double_t dPtMin, Double_t dPtMax)
{
  if(dPtMin >= dPtMax)
    return 0;
  TF1* funFD = new TF1("funFD", "pol4(0)", 0., 100.);
  funFD->SetParameter(0, 0.22);
  funFD->SetParameter(1, 0.0323);
  funFD->SetParameter(2, -0.0174);
  funFD->SetParameter(3, 0.00179);
  funFD->SetParameter(4, -5.68e-5);
  Double_t dPtLimit = 10;
  Double_t dResult = 0;
  Double_t dConst = 0.;
  if(dPtMin > dPtLimit)
    dResult = dConst;
  if(dPtMax <= dPtLimit)
    dResult = funFD->Integral(dPtMin, dPtMax) / (dPtMax - dPtMin);
  if(dPtMax > dPtLimit)
    dResult = (funFD->Integral(dPtMin, dPtLimit) + dConst * (dPtMax - dPtLimit)) / (dPtMax - dPtMin);
  delete funFD;
  return dResult;
}

Double_t GetMeanFeedDownLambdaPYTHIA(Double_t dPtMin, Double_t dPtMax)
{
  Double_t dFDConst = 0.142;
  return dFDConst;
}

Double_t GetMeanValue(TH1* his)
{
  if(!his)
    return -1;
  Double_t dSumVal = 0;
  Double_t dSumWeight = 0;
  for(Int_t iB = 1; iB <= his->GetNbinsX(); iB++)
  {
    dSumVal += TMath::Abs(his->GetBinCenter(iB)) * his->GetBinContent(iB);
//    dSumVal += his->GetBinCenter(iB)*his->GetBinContent(iB);
    dSumWeight += his->GetBinContent(iB);
  }
  if(dSumWeight == 0)
    return -1;
  return dSumVal / dSumWeight;
}

TH1D* GetScaledEfficiency(TH2* hisPtEtaMeasured, TH2* hisPtEtaEffBase, Int_t iSwitchPt)
{
  // Expected axes:
  // x - pT
  // y - eta
  // Notation:
  // i - eta bin (code: iEta)
  // m(i) -bin counts of measured uncorrected spectrum of particles of interest (code: hisPtEtaMeasured)
  // a(i) - bin counts of spectrum of associated inclusive particles (not in code)
  // g(i) - bin counts of spectrum of generated inclusive particles (not in code)
  // eff(i) - bin counts of efficiency of inclusive particles (hisPtEtaEffBase) = a(i)/g(i)
  // a_s(i) - bin counts of scaled spectrum of associated particles of interest (code: hisPtEtaRecScaled)
  // g_s(i) - bin counts of scaled spectrum of generated particles of interest (code: hisPtEtaGenScaled)
  // Procedure:
  // - for each i:
  //   - Calculate a_s(i) = m(i)
  //   - Calculate g_s(i) = a_s(i)/a(i)*g(i) = m(i)/eff(i)

  printf("GetScaledEfficiency: Start\n");
  if(!CompareAxes2D(hisPtEtaMeasured, hisPtEtaEffBase))
    return NULL;

  // corrected MC histograms
  TH2D* hisPtEtaRecScaled = (TH2D*)hisPtEtaMeasured->Clone(Form("%s-RecScaled", hisPtEtaMeasured->GetName()));
  hisPtEtaRecScaled->Reset();
  hisPtEtaRecScaled->Sumw2();
  TH2D* hisPtEtaGenScaled = (TH2D*)hisPtEtaMeasured->Clone(Form("%s-GenScaled", hisPtEtaMeasured->GetName()));
  hisPtEtaGenScaled->Reset();
  hisPtEtaGenScaled->Sumw2();

  Int_t iNBinsPtV0Tmp, iBinPtFirstTmp, iBinPtLastTmp;
  if(iSwitchPt == 0) // inclusive pt bins
  {
    iNBinsPtV0Tmp = iNBinsPtV0AllLF;
    iBinPtFirstTmp = iBinPtInclFirst;
    iBinPtLastTmp = iBinPtInclLast;
  }
  if(iSwitchPt == 1) // in-jet pt bins
  {
    iNBinsPtV0Tmp = iNBinsPtV0InJet;
    iBinPtFirstTmp = iBinPtInJetsFirst;
    iBinPtLastTmp = iBinPtInJetsLast;
  }

  for(Int_t iPt = iBinPtFirstTmp; iPt <= iBinPtLastTmp; iPt++)
  {
    for(Int_t iEta = 1; iEta <= hisPtEtaMeasured->GetNbinsY(); iEta++)
    {
      // calculate correct yield of associated: a_s(i) = m(i)
      Double_t dRecNewErr = 0;
      Double_t dRecNew = hisPtEtaMeasured->GetBinContent(iPt, iEta);
      if(dRecNew <= dEpsilon)
        continue;
      hisPtEtaRecScaled->SetBinContent(iPt, iEta, dRecNew);
      hisPtEtaRecScaled->SetBinError(iPt, iEta, dRecNewErr);
      // calculate correct yield of generated: g_s(i) = a_s(i)/eff(i)
      Double_t dGenNewErr = 0;
      Double_t dGenNew = DivideNumbersError(hisPtEtaRecScaled->GetBinContent(iPt, iEta), hisPtEtaRecScaled->GetBinError(iPt, iEta), hisPtEtaEffBase->GetBinContent(iPt, iEta), hisPtEtaEffBase->GetBinError(iPt, iEta), &dGenNewErr, 0);
      if(dGenNew == 0)
        dGenNewErr = 0;
      printf("Pt: %d/%d, Eta: %d/%d, New generated yield: %g +- %g (%g %%)\n", iPt, iNBinsPtV0Tmp, iEta, hisPtEtaMeasured->GetNbinsY(), dGenNew, dGenNewErr, 100 * dGenNewErr / (dGenNew + dEpsilon));
      hisPtEtaGenScaled->SetBinContent(iPt, iEta, dGenNew);
      hisPtEtaGenScaled->SetBinError(iPt, iEta, dGenNewErr);
    }
    // Compare associated spectra
  }
  TH1D* hisPtRecScaled = (TH1D*)hisPtEtaRecScaled->ProjectionX(Form("%s_px", hisPtEtaRecScaled->GetName()), 0, -1, "e"); // sum up along eta axis
  TH1D* hisPtGenScaled = (TH1D*)hisPtEtaGenScaled->ProjectionX(Form("%s_px", hisPtEtaGenScaled->GetName()), 0, -1, "e"); // sum up along eta axis
  TH1D* hisPtEffScaled = (TH1D*)DivideHistograms1D(hisPtRecScaled, hisPtGenScaled); // calculate eff as function of pT
//          dirOutSpectra->cd();
//          hisPtEtaRecScaled->Write(Form("fh2PtEtaInJetsRecScale%s_C%d-J%d",V0Name[V0Type].Data(),iCent,iJet));
//          hisPtEtaGenScaled->Write(Form("fh2PtEtaInJetsGenScale%s_C%d-J%d",V0Name[V0Type].Data(),iCent,iJet));
  printf("GetScaledEfficiency: End\n");
  return hisPtEffScaled;
}

void DrawResultsLoopSys(TString sNameFileIn, Int_t V0Type, TString sNameFileEff, TString sNameFileFD, TString sFlag, Int_t iMode)
{
  TString sCutName = "";
  Double_t dCutVal = 0;

  sFlag = "Default";
  DrawResults(sNameFileIn, V0Type, sNameFileEff, sNameFileFD, sFlag, iMode);

  for(Int_t iCut = 0; iCut < iNArrayVarSys; iCut++)
//  for(Int_t iCut = 0; iCut < 2; iCut++)
  {
    sCutName = sArrayVarSysName[iCut];

    printf("Variations for %s\n", sCutName.Data());
    for(Int_t iVar = 0; iVar < iArrayNVarSys[iCut]; iVar++)
//    for(Int_t iVar = 0; iVar < 1; iVar++)
    {
      dCutVal = dArrayVarSys[V0Type][iCut][iVar];
      printf("Variation for %s: %g\n", sCutName.Data(), dCutVal);
      sFlag = Form("%s%d", sCutName.Data(), iVar + 1);

      DrawResults(sNameFileIn, V0Type, sNameFileEff, sNameFileFD, sFlag, iMode);
    }
  }
}
