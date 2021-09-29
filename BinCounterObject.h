#ifndef BINCOUNTEROBJECT_H
#define BINCOUNTEROBJECT_H

#include "TMath.h"

#include "TH1D.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"

#include "TStyle.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"

class TFitResultPtr;

class BinCounterObject
{
public:
  BinCounterObject(TH1* hisIn, Double_t dMeanIn, Double_t dSigmaIn); // Default constructor
  virtual ~BinCounterObject(); // Default destructor
//    BinCounterObject(const BinCounterObject& other); // Copy constructor
//    BinCounterObject& operator=(const BinCounterObject& other); // Assignment operator

  Double_t GetMean() { return dMean; }
  Double_t GetMeanError() { return dMeanErr; }
  void SetMean(Double_t val) { dMean = val; }
  Double_t GetSigma() { return dSigma; }
  Double_t GetSigmaError() { return dSigmaErr; }
  void SetSigma(Double_t val) { dSigma = val; }
  Bool_t SetNSigmas(Double_t nSig, Double_t nBgIn, Double_t nBgOut);
  Bool_t FixRegions(Double_t dSigMin, Double_t dSigMax, Double_t dBgInMin, Double_t dBgInMax, Double_t dBgOutMin, Double_t dBgOutMax);
  void SetConsiderBackground(Bool_t val = kTRUE) {bConsiderBg = val; bSubtract = val;}
  void SetBackgroundSubtraction(Bool_t val = kTRUE) {bSubtract = val;}
  void SetFitOptionGlobal(TString val) {sOptionFitGlob = val;}
  void SetFitOptionSideBands(TString val) {sOptionFitSB = val;}
  Bool_t CheckHistogram();
  Bool_t SetRegions(Bool_t bSideBands = kTRUE);
  Bool_t CheckRegions(Bool_t bSideBands = kTRUE);
  Bool_t SetXLimits(Double_t xmin, Double_t xmax);
  Bool_t EstimateParameters(Int_t iDegPol = 0);
  Bool_t Run();
  void RunFitGlobal(Int_t iDegreePol, TString sOption);
  void RunFitSideBands(Int_t iDegreePol, TString sOption);
  void PrintFitResults(TFitResultPtr result, TF1* funFit = NULL);
  static Bool_t IsFitOK(TFitResultPtr result);
  Double_t GetSignalAndError(Double_t* dSigErr = NULL) {*dSigErr = dIntegralSignalErr; return dIntegralSignal;}
  Double_t GetBgAndError(Double_t* dBgErr = NULL) {*dBgErr = dIntegralBgErr; return dIntegralBg;}
  Double_t GetPurityAndError(Double_t* dPurErr = NULL) {*dPurErr = dPurityErr; return dPurity;}
  void Plot(TString sNameFig, TString sTitle = "", Bool_t bLogY = kTRUE, TString sSuffix = "png");
  void SetVerbose(Bool_t val = kTRUE) {bVerbose = val;}
  TF1* GetFunctionGlobalFit() {return funFitGlob;}
  TF1* GetFunctionSideBands() {return funFitBg;}
  TF1* GetFunctionBackground() {return funFitBgInteg;}
  void SetDegreePolMax(Int_t val) {iDegreePolMax = val;}
protected:
private:
  Bool_t bVerbose; // print all messages?
  // Input
  TH1D* hisInput; // input histogram
  TString sNameHis; // name of histogram
  Double_t dWidthBin; // width of one histogram bin (assume uniform binning)
  Double_t dXMin; // min of histogram axis
  Double_t dXMax; // max of histogram axis
  Double_t dXLimitMin; // min of global range
  Double_t dXLimitMax; // max of global range
  Double_t dNEntriesSigMin; // minimum average number of entries in signal region for fitting
  Double_t dNEntriesBg1Min; // minimum average number of entries in background region for fitting with linear function
  Double_t dNEntriesBg2Min; // minimum average number of entries in background region for fitting with quadratic function
  // Formulas
  TString sFormulaGlob; // pol3+Gauss, formula for the combined fit
  TString sFormulaBgInteg; // pol3, formula for bg function to integrate bg under peak
  // Functions
  static const Int_t iNParGlob = 7; // number of parameters of the global fit
  Int_t iDegreePolMax; // default degree of polynomial to fit background
  TF1* funFitGlob; // function for the combined fit
  TF1* funFitBg; // function for fitting side bands
  TF1* funFitBgInteg; // function for integrating bg under peak
  // Fit options
  TString sOptionFitGlob;
  TString sOptionFitSB;
  // Gauss parameters
  Double_t dMean; // Member variable "dMean"
  Double_t dMeanErr; // Member variable "dMean"
  Double_t dSigma; // Member variable "dSigma"
  Double_t dSigmaErr; // Member variable "dSigma"
  Double_t dArea; // Member variable "dArea"
  // Bg parameters
  Bool_t bConsiderBg; // should be background subtracted considered?
  Bool_t bSubtract; // should be background subtracted from integral in signal region?
  Int_t iDegreePolInit; // initial degree of polynomial for fitting bg (3 by default)
  Double_t dPar0; // constant
  Double_t dPar1; // linear
  Double_t dPar2; // quadratic
  Double_t dPar3; // cubic
  // Edges of regions expressed as multiples of sigma
  Double_t dNSigmaSig; // signal
  Double_t dNSigmaBgIn; // inner edge of side band
  Double_t dNSigmaBgOut; // outer edge of side band
  // Edges of regions
  Bool_t bFixedRegions; // indicator of manually set regions
  Double_t dXSigMin; // signal min
  Double_t dXSigMax; // signal max
  Double_t dXBgInMin; // left side band max
  Double_t dXBgInMax; // right side band min
  Double_t dXBgOutMin; // left side band min
  Double_t dXBgOutMax; // right side band max
  Int_t iBinSigMin; // bin corresponding to signal min
  Int_t iBinSigMax; // bin corresponding to signal max
  Int_t iBinBgInMin;
  Int_t iBinBgInMax;
  Int_t iBinBgOutMin;
  Int_t iBinBgOutMax;
  // Fit results
  TFitResultPtr resultFitGlobal; // global fit
  TFitResultPtr resultFitSideBands; // side band fit
  // Results
  Double_t dIntegralSignal;
  Double_t dIntegralSignalErr;
  Double_t dIntegralBg;
  Double_t dIntegralBgErr;
  Double_t dPurity;
  Double_t dPurityErr;
};

#endif // BINCOUNTEROBJECT_H
