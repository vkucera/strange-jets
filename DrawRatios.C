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

#include "AuxROOTFunctions.h"
#include "AuxFunctions.h"

TString sEnergy = "#sqrt{#it{s}_{NN}} = 2.76 TeV";
TString sSystem = "Pb#minusPb";
Float_t nCollCentBins[6] = {1502.7, 923.26, 438.8, 128.2, 26.82, 450.33}; // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies#Tables_with_centrality_bins_for
Float_t rcpErrorCent[4] = {0.1455895292, 0.1416009451, 0.1357282295, 0.1350063077}; // syst. error of ncoll ratios
Float_t maxYRho[3] = {200., 100., 40.};
Float_t etaNormalisation = 1.8;

// ================================
Float_t fRadiusJet = 0.2; // R
Float_t fDistanceV0Jet = 0.2; // D
// ================================

TString V0Name[3] = {"K0s", "Lambda", "ALambda"};
TString V0Symbol[3] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}"};
TString V0LabelM[3] = {"#it{m}_{inv} K_{S}^{0} (GeV/#it{c}^{2})", "#it{m}_{inv} #Lambda (GeV/#it{c}^{2})", "#it{m}_{inv} #bar{#Lambda} (GeV/#it{c}^{2})"};
Float_t fMassV0[3] = {0.497614, 1.115680, 1.115680};

// Global ranges
// centrality bins
Int_t iCentMin = 0;
//Int_t iCentMax = iNCentBins-1;
Int_t iCentMax = 0;
// jet pT bins
Int_t iJetMin = 1;
//Int_t iJetMax = iNBinsPtJet-1;
Int_t iJetMax = 2;

Bool_t bInclusive = 1;
Bool_t bInJets = 0;
Bool_t bInBulk = 0;
Bool_t bCorrel = 0;
Bool_t bMC = 0;
Bool_t bALambda = 1;
//{
Bool_t bCombineLaL = 1; // combine Lambda and anti-Lambda
//}

TString sLabelRatio = "";

void DrawRatios(TString sFlag = "", Int_t iMode = 0)
{
  gStyle->SetCanvasColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetOptStat(0);

  // modes: 0 uncorrected, 1 simulated, 2 corrected
  switch (iMode)
  {
    case 0:
      bMC = 0;
      break;
    case 1:
      bMC = 1;
      break;
    case 2:
      bMC = 0;
      break;
    default:
      printf("Error: Wrong mode.\n");
      return;
      break;
  }

  if(bMC)
  {
    bInclusive = 0;
    bInJets = 0;
    bInBulk = 0;
  }

  if (bALambda && bCombineLaL)
    sLabelRatio = "(#Lambda + #bar{#Lambda})/2K_{S}^{0}";
  else
    sLabelRatio = "#Lambda/K_{S}^{0}";

  TString sLabelCollisionText = Form("#splitline{%s}{%s}", sSystem.Data(), sEnergy.Data());
  TLatex* labelCollision = new TLatex();
  labelCollision->SetTextFont(42);
  labelCollision->SetTextSize(0.04);
  labelCollision->SetTextAlign(23);
  TLatex* labelPt;
  TLegend* legend;

  // canvas size
  Int_t iCanHeight = 600;
  Int_t iCanWidth = 600;
  TString kImageSuf = "png";
  //TString kImageSuf = "pdf";

  TString sPath = "";
//  sPath = (sFlag.Length() ? Form("%s/",sFlag.Data()) : "");
  TString sNameFileInK = "OutputK0s.root";
  TString sNameFileInL = "OutputLambda.root";
  TString sNameFileInAL = "OutputALambda.root";
  TString sNameDirPt = "Spectra";
  TString sNameFileOut = "OutputRatio.root";

  // input
  TString sNameHisInclusive = "fh1PtInclusive%s_C%d";
  TString sNameHisInJets = "fh1PtInJets%s_C%d-J%d";
  TString sNameHisInBulk = "fh1PtInBulk%s_C%d";
  TString sNameHisCorrel = "hisPtCorrelYield-C%d-J%d";
//  TString sNameHisGenInclusive = "fh1PtGenInclusive%s_C%d";
  TString sNameHisGenInclusive = "fh1GenPtIncl%s_%d";
  // output
  TString sNameHisRatioInclusive = "fh1PtRatioInclusive_C%d";
  TString sNameHisRatioInJets = "fh1PtRatioInJets_C%d-J%d";
  TString sNameHisRatioInBulk = "fh1PtRatioInBulk_C%d";
  TString sNameHisRatioCorrel = "fh1PtRatioCorrel_C%d_J%d";

  printf("Loading file %s ", sNameFileInK.Data());
  TFile* fileK = 0;
  fileK = new TFile(Form("%s%s",sPath.Data(), sNameFileInK.Data()), "OPEN");
  if(fileK->IsZombie())
  {
    printf("failed\nDrawRatios: Error: Cannot load file %s\n", sNameFileInK.Data());
    return;
  }
  printf("OK\n");

  printf("Loading directory %s ", sNameDirPt.Data());
  TDirectoryFile* dirK = (TDirectoryFile*)fileK->Get(sNameDirPt.Data());
  if(!dirK)
  {
    printf("failed\nDrawRatios: Error: no dir\n");
    return;
  }
  printf("OK\n");

  printf("Loading file %s ", sNameFileInL.Data());
  TFile* fileL = 0;
  fileL = new TFile(Form("%s%s",sPath.Data(), sNameFileInL.Data()), "OPEN");
  if(fileL->IsZombie())
  {
    printf("failed\nDrawRatios: Error: Cannot load file %s\n", sNameFileInL.Data());
    return;
  }
  printf("OK\n");

  printf("Loading directory %s ", sNameDirPt.Data());
  TDirectoryFile* dirL = (TDirectoryFile*)fileL->Get(sNameDirPt.Data());
  if(!dirL)
  {
    printf("failed\nDrawRatios: Error: no dir\n");
    return;
  }
  printf("OK\n");

  TFile* fileAL = 0;
  if(bALambda)
  {
    printf("Loading file %s ", sNameFileInAL.Data());
    fileAL = new TFile(Form("%s%s",sPath.Data(), sNameFileInAL.Data()), "OPEN");
    if(fileAL->IsZombie())
    {
      printf("failed\nDrawRatios: Error: Cannot load file %s\n", sNameFileInAL.Data());
      return;
    }
    printf("OK\n");
  }

  TDirectoryFile* dirAL;
  if(bALambda)
  {
    printf("Loading directory %s ", sNameDirPt.Data());
    dirAL = (TDirectoryFile*)fileAL->Get(sNameDirPt.Data());
    if(!dirAL)
    {
      printf("failed\nDrawRatios: Error: no dir\n");
      return;
    }
    printf("OK\n");
  }

  TFile* fileOut = new TFile(Form("%s%s",sPath.Data(), sNameFileOut.Data()), "RECREATE");

  Float_t fPtInJetsXMin = 0; // minimum x value
  Float_t fPtInJetsXMax = 12; // maximum x value

  gSystem->cd(sPath.Data());

  TCanvas* canRatioLambdaMC = new TCanvas("RatioLambdaMC", "", iCanWidth, iCanHeight);
  TMultiGraph* mgrRatioLambdaMC = new TMultiGraph();
  if(bMC && bALambda)
  {
    for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
    {
      TH1D* hisPtGenInclL = (TH1D*)dirL->Get(Form(sNameHisGenInclusive.Data(), V0Name[1].Data(), iCent));
      if(!hisPtGenInclL)
      {
        printf("DrawRatios: Error: Cannot load histogram %s\n", Form(sNameHisGenInclusive.Data(), V0Name[1].Data(), iCent));
        return;
      }
      TH1D* hisPtGenInclAL = (TH1D*)dirAL->Get(Form(sNameHisGenInclusive.Data(), V0Name[2].Data(), iCent));
      if(!hisPtGenInclAL)
      {
        printf("DrawRatios: Error: Cannot load histogram %s\n", Form(sNameHisGenInclusive.Data(), V0Name[2].Data(), iCent));
        return;
      }
      TH1D* hisRatioPtGenInclLambda = DivideHistograms1D(hisPtGenInclAL, hisPtGenInclL);
//          TH1D* hisRatioPtGenInclLambda = DivideHistograms1D(hisPtGenInclL,hisPtGenInclAL);
      TGraphErrors* grRatioPtGenInclLambda = MakeGraphErrors(hisRatioPtGenInclLambda, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersFull[iCent]);
      mgrRatioLambdaMC->Add(grRatioPtGenInclLambda);
    }
    canRatioLambdaMC->cd();
    canRatioLambdaMC->SetLeftMargin(0.15);
    mgrRatioLambdaMC->SetTitle("#bar{#Lambda}/#Lambda, generated;#it{p}_{T}^{h} (GeV/#it{c});ratio");
    mgrRatioLambdaMC->SetMinimum(0.);
    mgrRatioLambdaMC->SetMaximum(1.5);
    mgrRatioLambdaMC->Draw("AP0");
    mgrRatioLambdaMC->GetYaxis()->SetTitleOffset(2);
    mgrRatioLambdaMC->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
    legend = canRatioLambdaMC->BuildLegend(0.7, 0.6, 0.85, 0.85);
    SetLegend(legend);
//      labelPt=labelCollision->DrawLatex(2.5,1.9,sLabelCollisionText.Data());
    labelPt = labelCollision->DrawLatex(2, 1.3, sLabelCollisionText.Data());
    labelPt->SetTextSize(0.04);
    labelPt->SetTextAlign(23);
    canRatioLambdaMC->SaveAs(Form("canRatioLambdaMC.%s", kImageSuf.Data()));
    delete canRatioLambdaMC;
    delete mgrRatioLambdaMC;
  }

  TCanvas* canRatioCompare[iNCentBins];
  TMultiGraph* mgrRatioCompare[iNCentBins];
  TLegend* legRatioCompare[iNCentBins];
  for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
  {
    canRatioCompare[iCent] = new TCanvas(Form("canRatioCompare%d", iCent), "", iCanWidth, iCanHeight);
    mgrRatioCompare[iCent] = new TMultiGraph();
    legRatioCompare[iCent] = new TLegend(0.5, 0.65, 0.85, 0.87);
  }

  if(bInclusive)
  {
    TCanvas* canRatioInclusive = new TCanvas("RatioInclusive", "", iCanWidth, iCanHeight);
    TMultiGraph* mgrRatioInclusive = new TMultiGraph();
    TCanvas* canRatioLambdaInclusive = new TCanvas("RatioLambdaInclusive", "", iCanWidth, iCanHeight);
    TMultiGraph* mgrRatioLambdaInclusive = new TMultiGraph();
    for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
    {
      TH1D* hisPtInclusiveK = (TH1D*)dirK->Get(Form(sNameHisInclusive.Data(), V0Name[0].Data(), iCent));
      if(!hisPtInclusiveK)
      {
        printf("DrawRatios: Error: Cannot load histogram %s\n", Form(sNameHisInclusive.Data(), V0Name[0].Data(), iCent));
        return;
      }
      TH1D* hisPtInclusiveL = (TH1D*)dirL->Get(Form(sNameHisInclusive.Data(), V0Name[1].Data(), iCent));
      if(!hisPtInclusiveL)
      {
        printf("DrawRatios: Error: Cannot load histogram %s\n", Form(sNameHisInclusive.Data(), V0Name[1].Data(), iCent));
        return;
      }
      TH1D* hisPtInclusiveAL = 0;
      if(bALambda)
      {
        hisPtInclusiveAL = (TH1D*)dirAL->Get(Form(sNameHisInclusive.Data(), V0Name[2].Data(), iCent));
        if(!hisPtInclusiveAL)
        {
          printf("DrawRatios: Error: Cannot load histogram %s\n", Form(sNameHisInclusive.Data(), V0Name[2].Data(), iCent));
          return;
        }
        TH1D* hisRatioLambdaInclusive = DivideHistograms1D(hisPtInclusiveAL, hisPtInclusiveL);
        TGraphErrors* grRatioLambdaInclusive = MakeGraphErrors(hisRatioLambdaInclusive, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersFull[iCent]);
        mgrRatioLambdaInclusive->Add(grRatioLambdaInclusive);
      }
      // combine Lambda + anti-Lambda
      if(bALambda && bCombineLaL)
      {
        hisPtInclusiveL->Add(hisPtInclusiveAL);
        hisPtInclusiveL = DivideHistogram(hisPtInclusiveL,2,0);
      }
      // calculate inclusive ratio
      TH1D* hisRatioInclusive = DivideHistograms1D(hisPtInclusiveL, hisPtInclusiveK);
      fileOut->cd();
      hisRatioInclusive->Write(Form(sNameHisRatioInclusive.Data(), iCent));
      TGraphErrors* grRatioInclusive = MakeGraphErrors(hisRatioInclusive, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersFull[iCent]);
      mgrRatioInclusive->Add(grRatioInclusive);
      TGraphErrors* grRatioInclusiveC = MakeGraphErrors(hisRatioInclusive, "inclusive", iMyColors[0], iMyMarkersFull[iCent]);
      mgrRatioCompare[iCent]->Add(grRatioInclusiveC);
      legRatioCompare[iCent]->AddEntry(grRatioInclusiveC);
    }

    canRatioInclusive->cd();
    canRatioInclusive->SetLeftMargin(0.15);
    mgrRatioInclusive->SetTitle(Form("Baryon-meson ratio, inclusive;#it{p}_{T}^{h} (GeV/#it{c});%s",sLabelRatio.Data()));
    mgrRatioInclusive->SetMinimum(0);
    mgrRatioInclusive->SetMaximum(2);
    mgrRatioInclusive->Draw("AP0");
    mgrRatioInclusive->GetYaxis()->SetTitleOffset(2);
    mgrRatioInclusive->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
    legend = canRatioInclusive->BuildLegend(0.7, 0.6, 0.85, 0.85);
    SetLegend(legend);
//      labelPt=labelCollision->DrawLatex(2.5,1.9,sLabelCollisionText.Data());
    labelPt = labelCollision->DrawLatex(8, 1.1, sLabelCollisionText.Data());
    labelPt->SetTextSize(0.04);
    labelPt->SetTextAlign(23);
    canRatioInclusive->SaveAs(Form("canRatioInclusive.%s", kImageSuf.Data()));
    delete canRatioInclusive;
    delete mgrRatioInclusive;

    if(bALambda)
    {
      canRatioLambdaInclusive->cd();
      canRatioLambdaInclusive->SetLeftMargin(0.15);
      mgrRatioLambdaInclusive->SetTitle("#bar{#Lambda}/#Lambda, inclusive;#it{p}_{T}^{h} (GeV/#it{c});ratio");
      mgrRatioLambdaInclusive->SetMinimum(0.5);
      mgrRatioLambdaInclusive->SetMaximum(1.5);
      mgrRatioLambdaInclusive->Draw("AP0");
      mgrRatioLambdaInclusive->GetYaxis()->SetTitleOffset(2);
      mgrRatioLambdaInclusive->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
      legend = canRatioLambdaInclusive->BuildLegend(0.7, 0.6, 0.85, 0.85);
      SetLegend(legend);
//      labelPt=labelCollision->DrawLatex(2.5,1.9,sLabelCollisionText.Data());
      labelPt = labelCollision->DrawLatex(2, 1.3, sLabelCollisionText.Data());
      labelPt->SetTextSize(0.04);
      labelPt->SetTextAlign(23);
      canRatioLambdaInclusive->SaveAs(Form("canRatioLambdaInclusive.%s", kImageSuf.Data()));
      delete canRatioLambdaInclusive;
      delete mgrRatioLambdaInclusive;
    }
  }
  if(bInJets)
  {
    TCanvas* canRatioInJets[iNBinsPtJet];
    TMultiGraph* mgrRatioInJets[iNBinsPtJet];
    TCanvas* canRatioLambdaInJets[iNBinsPtJet];
    TMultiGraph* mgrRatioLambdaInJets[iNBinsPtJet];
    for(Int_t iJet = iJetMin; iJet <= iJetMax; iJet++)
    {
      canRatioInJets[iJet] = new TCanvas(Form("RatioInJets%d", iJet), "", iCanWidth, iCanHeight);
      mgrRatioInJets[iJet] = new TMultiGraph();
      canRatioLambdaInJets[iJet] = new TCanvas(Form("RatioLambdaInJets%d", iJet), "", iCanWidth, iCanHeight);
      mgrRatioLambdaInJets[iJet] = new TMultiGraph();

      for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
      {
        TH1D* hisPtKInJets = (TH1D*)dirK->Get(Form(sNameHisInJets.Data(), V0Name[0].Data(), iCent, iJet));
        if(!hisPtKInJets)
        {
          printf("DrawRatios: Error: Cannot load histogram %s\n", Form(sNameHisInJets.Data(), V0Name[0].Data(), iCent, iJet));
          return;
        }
        TH1D* hisPtLInJets = (TH1D*)dirL->Get(Form(sNameHisInJets.Data(), V0Name[1].Data(), iCent, iJet));
        if(!hisPtLInJets)
        {
          printf("DrawRatios: Error: Cannot load histogram %s\n", Form(sNameHisInJets.Data(), V0Name[1].Data(), iCent, iJet));
          return;
        }
        TH1D* hisPtALInJets;
        if(bALambda)
        {
          hisPtALInJets = (TH1D*)dirAL->Get(Form(sNameHisInJets.Data(), V0Name[2].Data(), iCent, iJet));
          if(!hisPtALInJets)
          {
            printf("DrawRatios: Error: Cannot load histogram %s\n", Form(sNameHisInJets.Data(), V0Name[2].Data(), iCent, iJet));
            return;
          }
        }
        TH1* hisRatioLambdaInJets;
        if(bALambda)
          hisRatioLambdaInJets = DivideHistograms1D(hisPtALInJets, hisPtLInJets);
        // combine Lambda + anti-Lambda
        if(bALambda && bCombineLaL)
        {
          hisPtLInJets->Add(hisPtALInJets);
          hisPtLInJets = DivideHistogram(hisPtLInJets,2,0);
        }
        TH1* hisRatioInJets = DivideHistograms1D(hisPtLInJets, hisPtKInJets);
        fileOut->cd();
        hisRatioInJets->Write(Form(sNameHisRatioInJets.Data(), iCent, iJet));
        TGraphErrors* grRatioInJets = MakeGraphErrors(hisRatioInJets, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersFull[iCent]);
        mgrRatioInJets[iJet]->Add(grRatioInJets);
        if(bALambda)
        {
          TGraphErrors* grRatioLambdaInJets = MakeGraphErrors(hisRatioLambdaInJets, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersFull[iCent]);
          mgrRatioLambdaInJets[iJet]->Add(grRatioLambdaInJets);
        }
        TGraphErrors* grRatioInJetsC = MakeGraphErrors(hisRatioInJets, Form("jets: #it{p}_{T}^{jet,ch} > %g GeV/#it{c}", dBinsPtJet[iJet]), iMyColors[iJet + 1], iMyMarkersFull[iCent]);
        TGraphErrors* grRatioInJetsCSys = MakeGraphErrors(hisRatioInJets, Form("jets: #it{p}_{T}^{jet,ch} > %g GeV/#it{c}", dBinsPtJet[iJet]), iMyColors[iJet + 1], iMyMarkersFull[iCent]);
        mgrRatioCompare[iCent]->Add(grRatioInJetsC);
        legRatioCompare[iCent]->AddEntry(grRatioInJetsC);
        grRatioInJetsCSys->SetFillStyle(0);
        grRatioInJetsCSys->SetLineColor(iMyColors[iJet + 1]);
//              mgrRatioCompare[iCent]->Add(grRatioInJetsCSys,"2");
      }
      canRatioInJets[iJet]->cd();
      canRatioInJets[iJet]->SetLeftMargin(0.15);
//          mgrRatioInJets[iJet]->SetTitle(Form("Baryon/meson ratio, in jets %1.f-%1.f GeV/#it{c};#it{p}_{T}^{h} (GeV/#it{c});#Lambda/K_{S}^{0}",dBinsPtJet[iJet],dBinsPtJet[iJet+1]));
      mgrRatioInJets[iJet]->SetTitle(Form("Baryon-meson ratio, in jets > %1.f GeV/#it{c};#it{p}_{T}^{h} (GeV/#it{c});%s", dBinsPtJet[iJet], sLabelRatio.Data()));
      mgrRatioInJets[iJet]->SetMinimum(0);
      mgrRatioInJets[iJet]->SetMaximum(2);
      mgrRatioInJets[iJet]->Draw("AP0");
      mgrRatioInJets[iJet]->GetYaxis()->SetTitleOffset(2);
      mgrRatioInJets[iJet]->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
      legend = canRatioInJets[iJet]->BuildLegend(0.7, 0.6, 0.85, 0.85);
      SetLegend(legend);
//          labelPt=labelCollision->DrawLatex(2.5,1.9,sLabelCollisionText.Data());
      labelPt = labelCollision->DrawLatex(8, 1.1, sLabelCollisionText.Data());
      labelPt->SetTextSize(0.04);
      labelPt->SetTextAlign(23);
      canRatioInJets[iJet]->SaveAs(Form("canRatioInJets_J%d.%s", iJet, kImageSuf.Data()));
      delete canRatioInJets[iJet];
      delete mgrRatioInJets[iJet];

      if(bALambda)
      {
        canRatioLambdaInJets[iJet]->cd();
        canRatioLambdaInJets[iJet]->SetLeftMargin(0.15);
//          mgrRatioLambdaInJets[iJet]->SetTitle(Form("#Lambda/#bar{#Lambda}, in jets %1.f-%1.f GeV/#it{c};#it{p}_{T}^{h} (GeV/#it{c});ratio",dBinsPtJet[iJet],dBinsPtJet[iJet+1]));
        mgrRatioLambdaInJets[iJet]->SetTitle(Form("#Lambda/#bar{#Lambda}, in jets > %1.f GeV/#it{c};#it{p}_{T}^{h} (GeV/#it{c});ratio", dBinsPtJet[iJet]));
        mgrRatioLambdaInJets[iJet]->SetMinimum(0.5);
        mgrRatioLambdaInJets[iJet]->SetMaximum(1.5);
        mgrRatioLambdaInJets[iJet]->Draw("AP0");
        mgrRatioLambdaInJets[iJet]->GetYaxis()->SetTitleOffset(2);
        mgrRatioLambdaInJets[iJet]->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
        legend = canRatioLambdaInJets[iJet]->BuildLegend(0.7, 0.6, 0.85, 0.85);
        SetLegend(legend);
//          labelPt=labelCollision->DrawLatex(2.5,1.9,sLabelCollisionText.Data());
        labelPt = labelCollision->DrawLatex(2, 1.3, sLabelCollisionText.Data());
        labelPt->SetTextSize(0.04);
        labelPt->SetTextAlign(23);
        canRatioLambdaInJets[iJet]->SaveAs(Form("canRatioLambdaInJets_J%d.%s", iJet, kImageSuf.Data()));
        delete canRatioLambdaInJets[iJet];
        delete mgrRatioLambdaInJets[iJet];
      }
    }
  }
  if(bInBulk)
  {
    TCanvas* canRatioInBulk = new TCanvas("RatioInBulk", "", iCanWidth, iCanHeight);
    TMultiGraph* mgrRatioInBulk = new TMultiGraph();
    TCanvas* canRatioLambdaInBulk = new TCanvas("RatioLambdaInBulk", "", iCanWidth, iCanHeight);
    TMultiGraph* mgrRatioLambdaInBulk = new TMultiGraph();
    for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
    {
      TH1D* hisPtKInBulk = (TH1D*)dirK->Get(Form(sNameHisInBulk.Data(), V0Name[0].Data(), iCent));
      if(!hisPtKInBulk)
      {
        printf("DrawRatios: Error: Cannot load histogram %s\n", Form(sNameHisInBulk.Data(), V0Name[0].Data(), iCent));
        return;
      }
      TH1D* hisPtLInBulk = (TH1D*)dirL->Get(Form(sNameHisInBulk.Data(), V0Name[1].Data(), iCent));
      if(!hisPtLInBulk)
      {
        printf("DrawRatios: Error: Cannot load histogram %s\n", Form(sNameHisInBulk.Data(), V0Name[1].Data(), iCent));
        return;
      }
      TH1D* hisPtALInBulk;
      if(bALambda)
      {
        hisPtALInBulk = (TH1D*)dirAL->Get(Form(sNameHisInBulk.Data(), V0Name[2].Data(), iCent));
        if(!hisPtALInBulk)
        {
          printf("DrawRatios: Error: Cannot load histogram %s\n", Form(sNameHisInBulk.Data(), V0Name[2].Data(), iCent));
          return;
        }
      }
      TH1* hisRatioLambdaInBulk;
      if(bALambda)
        hisRatioLambdaInBulk = DivideHistograms1D(hisPtALInBulk, hisPtLInBulk);
      // combine Lambda + anti-Lambda
      if(bALambda && bCombineLaL)
      {
        hisPtLInBulk->Add(hisPtALInBulk);
        hisPtLInBulk = DivideHistogram(hisPtLInBulk,2,0);
      }
      TH1* hisRatioInBulk = DivideHistograms1D(hisPtLInBulk, hisPtKInBulk);
      fileOut->cd();
      hisRatioInBulk->Write(Form(sNameHisRatioInBulk.Data(), iCent));
      TGraphErrors* grRatioInBulk = MakeGraphErrors(hisRatioInBulk, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersEmpty[iCent]);
      mgrRatioInBulk->Add(grRatioInBulk);
      if(bALambda)
      {
        TGraphErrors* grRatioLambdaInBulk = MakeGraphErrors(hisRatioLambdaInBulk, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersEmpty[iCent]);
        mgrRatioLambdaInBulk->Add(grRatioLambdaInBulk);
      }
      TGraphErrors* grRatioInBulkC = MakeGraphErrors(hisRatioInBulk, Form("underlying event"), iMyColors[1], iMyMarkersEmpty[iCent]);
      mgrRatioCompare[iCent]->Add(grRatioInBulkC);
      legRatioCompare[iCent]->AddEntry(grRatioInBulkC);
    }
    canRatioInBulk->cd();
    canRatioInBulk->SetLeftMargin(0.15);
    mgrRatioInBulk->SetTitle(Form("Baryon-meson ratio, in underlying event;#it{p}_{T}^{h} (GeV/#it{c});%s", sLabelRatio.Data()));
    mgrRatioInBulk->SetMinimum(0);
    mgrRatioInBulk->SetMaximum(2);
    mgrRatioInBulk->Draw("AP0");
    mgrRatioInBulk->GetYaxis()->SetTitleOffset(2);
    mgrRatioInBulk->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
    legend = canRatioInBulk->BuildLegend(0.7, 0.6, 0.85, 0.85);
    SetLegend(legend);
//          labelPt=labelCollision->DrawLatex(2.5,1.9,sLabelCollisionText.Data());
    labelPt = labelCollision->DrawLatex(8, 1.1, sLabelCollisionText.Data());
    labelPt->SetTextSize(0.04);
    labelPt->SetTextAlign(23);
    canRatioInBulk->SaveAs(Form("canRatioInBulk.%s", kImageSuf.Data()));
    delete canRatioInBulk;
    delete mgrRatioInBulk;

    if(bALambda)
    {
      canRatioLambdaInBulk->cd();
      canRatioLambdaInBulk->SetLeftMargin(0.15);
      mgrRatioLambdaInBulk->SetTitle(Form("#bar{#Lambda}/#Lambda, in underlying event;#it{p}_{T}^{h} (GeV/#it{c});ratio"));
      mgrRatioLambdaInBulk->SetMinimum(0.5);
      mgrRatioLambdaInBulk->SetMaximum(1.5);
      mgrRatioLambdaInBulk->Draw("AP0");
      mgrRatioLambdaInBulk->GetYaxis()->SetTitleOffset(2);
      mgrRatioLambdaInBulk->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
      legend = canRatioLambdaInBulk->BuildLegend(0.7, 0.6, 0.85, 0.85);
      SetLegend(legend);
//          labelPt=labelCollision->DrawLatex(2.5,1.9,sLabelCollisionText.Data());
      labelPt = labelCollision->DrawLatex(2, 1.3, sLabelCollisionText.Data());
      labelPt->SetTextSize(0.04);
      labelPt->SetTextAlign(23);
      canRatioLambdaInBulk->SaveAs(Form("canRatioLambdaInBulk.%s", kImageSuf.Data()));
      delete canRatioLambdaInBulk;
      delete mgrRatioLambdaInBulk;
    }
  }
  if(bCorrel)
  {
    TCanvas* canRatioCorrel[iNBinsPtJet];
    TMultiGraph* mgrRatioCorrel[iNBinsPtJet];
    for(Int_t iJet = iJetMin; iJet <= iJetMax; iJet++)
    {
      canRatioCorrel[iJet] = new TCanvas(Form("RatioCorrel%d", iJet), "", iCanWidth, iCanHeight);
      mgrRatioCorrel[iJet] = new TMultiGraph();

      for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
      {
        TH1D* hisPtKCorrel = (TH1D*)dirK->Get(Form(sNameHisCorrel.Data(), iCent, iJet));
        if(!hisPtKCorrel)
        {
          printf("DrawRatios: Error: Cannot load histogram %s\n", Form(sNameHisCorrel.Data(), iCent, iJet));
          return;
        }
        TH1D* hisPtLCorrel = (TH1D*)dirL->Get(Form(sNameHisCorrel.Data(), iCent, iJet));
        if(!hisPtLCorrel)
        {
          printf("DrawRatios: Error: Cannot load histogram %s\n", Form(sNameHisCorrel.Data(), iCent, iJet));
          return;
        }
        TH1* hisRatioCorrel = DivideHistograms1D(hisPtLCorrel, hisPtKCorrel);
        fileOut->cd();
        hisRatioCorrel->Write(Form(sNameHisRatioCorrel.Data(), iCent, iJet));
        TGraphErrors* grRatioCorrel = MakeGraphErrors(hisRatioCorrel, GetCentBinLabel(iCent).Data(), iMyColors[iCent], iMyMarkersFull[iCent]);
        mgrRatioCorrel[iJet]->Add(grRatioCorrel);
        TGraphErrors* grRatioCorrelC = MakeGraphErrors(hisRatioCorrel, Form("Cjets: #it{p}_{T}^{jet,ch} > %g GeV/#it{c}", dBinsPtJet[iJet]), iMyColors[iJet + 1], iMyMarkersEmpty[iCent]);
        TGraphErrors* grRatioCorrelCSys = MakeGraphErrors(hisRatioCorrel, Form("Cjets: #it{p}_{T}^{jet,ch} > %g GeV/#it{c}", dBinsPtJet[iJet]), iMyColors[iJet + 1], iMyMarkersEmpty[iCent]);
        mgrRatioCompare[iCent]->Add(grRatioCorrelC);
        legRatioCompare[iCent]->AddEntry(grRatioCorrelC);
        grRatioCorrelCSys->SetFillStyle(0);
        grRatioCorrelCSys->SetLineColor(iMyColors[iJet + 1]);
      }
      canRatioCorrel[iJet]->cd();
      canRatioCorrel[iJet]->SetLeftMargin(0.15);
      mgrRatioCorrel[iJet]->SetTitle(Form("Baryon-meson ratio, in jets > %g GeV/#it{c};#it{p}_{T}^{h} (GeV/#it{c});%s", dBinsPtJet[iJet],sLabelRatio.Data()));
      mgrRatioCorrel[iJet]->SetMinimum(0);
      mgrRatioCorrel[iJet]->SetMaximum(2);
      mgrRatioCorrel[iJet]->Draw("AP0");
      mgrRatioCorrel[iJet]->GetYaxis()->SetTitleOffset(2);
      mgrRatioCorrel[iJet]->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
      legend = canRatioCorrel[iJet]->BuildLegend(0.7, 0.6, 0.85, 0.85);
      SetLegend(legend);
      labelPt = labelCollision->DrawLatex(8, 1.1, sLabelCollisionText.Data());
      labelPt->SetTextSize(0.04);
      labelPt->SetTextAlign(23);
      canRatioCorrel[iJet]->SaveAs(Form("canRatioCorrel_J%d.%s", iJet, kImageSuf.Data()));
      delete canRatioCorrel[iJet];
      delete mgrRatioCorrel[iJet];

    }
  }

  if(bInclusive || bInJets || bInBulk || bCorrel)
  {
    for(Int_t iCent = iCentMin; iCent <= iCentMax; iCent++)
    {
      canRatioCompare[iCent]->cd();
      canRatioCompare[iCent]->SetLeftMargin(0.15);
      canRatioCompare[iCent]->SetLeftMargin(0.15);

      mgrRatioCompare[iCent]->SetTitle(Form("Baryon-meson ratios, c. %s;#it{p}_{T}^{hadron} (GeV/#it{c});%s", GetCentBinLabel(iCent).Data(),sLabelRatio.Data()));
      mgrRatioCompare[iCent]->SetMinimum(0);
      mgrRatioCompare[iCent]->SetMaximum(2);
      mgrRatioCompare[iCent]->Draw("AP0");
      mgrRatioCompare[iCent]->GetYaxis()->SetTitleOffset(2);
      mgrRatioCompare[iCent]->GetXaxis()->SetLimits(fPtInJetsXMin, fPtInJetsXMax);
//          legend = canRatioCompare[iCent]->BuildLegend(0.5,0.65,0.85,0.87);
//          SetLegend(legend);
      legRatioCompare[iCent]->Draw();
      SetLegend(legRatioCompare[iCent]);

//      labelPt=labelCollision->DrawLatex(2.5,1.9,sLabelCollisionText.Data());
      labelPt = labelCollision->DrawLatex(9, 1.3, sLabelCollisionText.Data());
      labelPt->SetTextSize(0.04);
      labelPt->SetTextAlign(23);
      labelPt = labelCollision->DrawLatex(9, 1, Form("#it{R} = %g, #it{D} = %g", fRadiusJet, fDistanceV0Jet));
      labelPt->SetTextSize(0.04);
      labelPt->SetTextAlign(23);
      canRatioCompare[iCent]->SaveAs(Form("canRatioCompare_C%d.%s", iCent, kImageSuf.Data()));
      delete canRatioCompare[iCent];
      delete mgrRatioCompare[iCent];
    }
  }
  fileOut->Close();
  fileK->Close();
  fileL->Close();
}
