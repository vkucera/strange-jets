#include "BinCounterObject.h"

Double_t BgFit(Double_t* x, Double_t* par);
void PrintMatrix(TMatrixDSym matrix);

BinCounterObject::BinCounterObject(TH1* hisIn, Double_t dMeanIn, Double_t dSigmaIn)
{
  //ctor
  bVerbose = kFALSE;
  hisInput = 0; // input histogram
  sNameHis = ""; // name of histogram
  dWidthBin = 0; // width of one histogram bin (assume uniform binning)
  dXMin = 0; // min of histogram axis
  dXMax = 0; // max of histogram axis
  dXLimitMin = 0; // min of range
  dXLimitMax = 0; // max of range
  dNEntriesSigMin = 1; // minimum average number of entries in signal region for fitting
  dNEntriesBg1Min = 1; // minimum average number of entries in background region for fitting with linear function
  dNEntriesBg2Min = 5; // minimum average number of entries in background region for fitting with quadratic function
  // Formulas
  sFormulaGlob = "(([4]/([6]*pow(6.2831853,0.5)))*exp(-0.5*pow((x-[5])/[6],2))) + [0] + [1]*x + [2]*x*x + [3]*x*x*x"; // pol3+Gauss, formula for the combined fit
  sFormulaBgInteg = "[0] + [1]*x + [2]*x*x + [3]*x*x*x + 0*[4] + 0*[5]"; // pol3, formula for bg function to integrate bg under peak
  // Functions
  funFitGlob = 0; // function for the combined fit
  funFitBg = 0; // function for fitting side bands
  funFitBgInteg = 0; // function for integrating bg under peak
  // Fit options
  sOptionFitGlob = "SLRI";
  sOptionFitSB = "SLRINE+"; // "SLRNFE+"; // crash
  // Gauss parameters
  dMean = 0;
  dMeanErr = 0;
  dSigma = 0;
  dSigmaErr = 0;
  dArea = 1;
  // Bg parameters
  bConsiderBg = kTRUE;
  bSubtract = kTRUE;
  iDegreePolMax = 3;
  iDegreePolInit = iDegreePolMax; // initial degree of polynomial for fitting bg (3 by default)
  dPar0 = 0; // constant
  dPar1 = 0; // linear
  dPar2 = 0; // quadratic
  dPar3 = 0; // cubic
  // Edges of regions expressed as multiples of sigma
  dNSigmaSig = 5;
  dNSigmaBgIn = 6;
  dNSigmaBgOut = 8;
  // Edges of regions
  bFixedRegions = kFALSE; // indicator of manually set regions
  dXSigMin = 0; // signal min
  dXSigMax = 0; // signal max
  dXBgInMin = 0; // left side band max
  dXBgInMax = 0; // right side band min
  dXBgOutMin = 0; // left side band min
  dXBgOutMax = 0; // right side band max
  // Fit results
  resultFitGlobal = 0; // global fit
  resultFitSideBands = 0; // side band fit
  // Results
  dIntegralSignal = 0;
  dIntegralSignalErr = 0;
  dIntegralBg = 0;
  dIntegralBgErr = 0;
  dPurity = 0;
  dPurityErr = 0;

  if(!hisIn)
  {
    CheckHistogram();
    return;
  }
  sNameHis = hisIn->GetName();
  printf("BinCounterObject::BinCounterObject: %s\n", sNameHis.Data());
  hisInput = (TH1D*)hisIn->Clone(Form("hisBC-%s", sNameHis.Data())); // copy of input histogram
  Int_t iNBinsX = hisInput->GetNbinsX();
  dXMin = hisInput->GetXaxis()->GetBinLowEdge(1);
  dXMax = hisInput->GetXaxis()->GetBinUpEdge(iNBinsX);
  dWidthBin = hisInput->GetBinWidth(1);
  dXLimitMin = dXMin;
  dXLimitMax = dXMax;
  funFitGlob = new TF1("funFitGlob", sFormulaGlob.Data(), dXLimitMin, dXLimitMax);
  funFitBg = new TF1("funFitBg", BgFit, dXLimitMin, dXLimitMax, 6);
  funFitBgInteg = new TF1("funFitBgInteg", sFormulaBgInteg.Data(), dXLimitMin, dXLimitMax);
  funFitGlob->SetParNames("Constant", "Linear", "Quadratic", "Cubic", "Area", "Mean", "Sigma");
  funFitBg->SetParNames("Constant", "Linear", "Quadratic", "Cubic", "ExcludeMin", "ExcludeMax");

  dMean = dMeanIn;
  if(dMeanIn < dXMin || dMeanIn > dXMax)
  {
    printf("BinCounterObject::BinCounterObject: Error: Mean outside range of histogram, setting middle of axis\n");
    dMean = (dXMax + dXMin) / 2.;
  }
  dSigma = dSigmaIn;
  Double_t dSigmaMax = (dXMax - dXMin) / 2.;
  if(dSigma > dSigmaMax)
  {
    printf("BinCounterObject::BinCounterObject: Error: Sigma too large, setting half of the axis length\n");
    dSigma = dSigmaMax;
  }
  dPar0 = hisInput->GetBinContent(iNBinsX / 2) / 2.;
  funFitGlob->SetParameters(dPar0, dPar1, dPar2, dPar3, dArea, dMean, dSigma);
  funFitBg->SetParameters(dPar0, dPar1, dPar2, dPar3, dXBgInMin, dXBgInMax);
  funFitBgInteg->SetParameters(dPar0, dPar1, dPar2, dPar3);
  SetRegions();
}

BinCounterObject::~BinCounterObject()
{
  //dtor
  if(hisInput)
    delete hisInput;
  if(funFitGlob)
    delete funFitGlob;
  if(funFitBg)
    delete funFitBg;
  if(funFitBgInteg)
    delete funFitBgInteg;
}
/*
BinCounterObject::BinCounterObject(const BinCounterObject& other)
{
  //copy ctor
}

BinCounterObject& BinCounterObject::operator=(const BinCounterObject& rhs)
{
  if (this == &rhs) return *this; // handle self assignment
  //assignment operator
  return *this;
}
*/
Bool_t BinCounterObject::CheckHistogram()
{
  if(!hisInput)
  {
    printf("BinCounterObject::CheckHistogram: Error: Invalid histogram!\n");
    hisInput = 0;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t BinCounterObject::SetXLimits(Double_t xmin, Double_t xmax)
{
  // check values
  if(bVerbose)
    printf("BinCounterObject::SetXLimits: Start\n");
  if(xmin >= xmax)
  {
    printf("BinCounterObject::SetXLimits: Error: Wrong range: %g - %g\n", xmin, xmax);
    return kFALSE;
  }
  if((xmin < dXMin) || (xmax > dXMax))
    printf("BinCounterObject::SetXLimits: Warning: range exceeds range of histogram: %g - %g, %g - %g, trimming\n", dXMin, xmin, xmax, dXMax);
  dXLimitMin = TMath::Max(xmin, dXMin);
  dXLimitMax = TMath::Min(xmax, dXMax);
  if(bVerbose)
    printf("BinCounterObject::SetXLimits: End\n");
  return kTRUE;
}

Bool_t BinCounterObject::SetNSigmas(Double_t nSig, Double_t nBgIn, Double_t nBgOut)
{
  if(bVerbose)
    printf("BinCounterObject::SetNSigmas: Start\n");
  Bool_t bStatus = kTRUE;
  if(nSig <= 0)
  {
    printf("BinCounterObject::SetNSigmas: Error: Wrong value of Sigma Signal: %g\n", nSig);
    bStatus = kFALSE;
  }
  if(nBgIn <= 0)
  {
    printf("BinCounterObject::SetNSigmas: Error: Wrong value of Sigma BgIn: %g\n", nBgIn);
    bStatus = kFALSE;
  }
  if(nBgOut <= 0)
  {
    printf("BinCounterObject::SetNSigmas: Error: Wrong value of Sigma BgOut: %g\n", nBgOut);
    bStatus = kFALSE;
  }
  if(!bStatus)
    return kFALSE;
  dNSigmaSig = nSig;
  dNSigmaBgIn = nBgIn;
  dNSigmaBgOut = nBgOut;
  bFixedRegions = kFALSE;
  if(!SetRegions())
    return kFALSE;
  if(bVerbose)
    printf("BinCounterObject::SetNSigmas: End\n");
  return kTRUE;
}

Bool_t BinCounterObject::FixRegions(Double_t dSigMin, Double_t dSigMax, Double_t dBgInMin, Double_t dBgInMax, Double_t dBgOutMin, Double_t dBgOutMax)
{
  if(bVerbose)
    printf("BinCounterObject::FixRegions: Start\n");
  dXSigMin = dSigMin;
  dXSigMax = dSigMax;
  dXBgInMin = dBgInMin;
  dXBgInMax = dBgInMax;
  dXBgOutMin = dBgOutMin;
  dXBgOutMax = dBgOutMax;
  bFixedRegions = kTRUE;
  dNSigmaSig = -1;
  dNSigmaBgIn = -1;
  dNSigmaBgOut = -1;
  if(!SetRegions())
    return kFALSE;
  if(bVerbose)
    printf("BinCounterObject::FixRegions: End\n");
  return kTRUE;
}

Bool_t BinCounterObject::CheckRegions(Bool_t bSideBands)
{
  if(bVerbose)
    printf("BinCounterObject::CheckRegions: Start\n");
  if(!CheckHistogram())
    return kFALSE;
  Bool_t bStatus = kTRUE;
  // trimm if outside of histogram
  if((dXSigMin < dXLimitMin) || (dXSigMax > dXLimitMax))
    printf("BinCounterObject::CheckRegions: Warning: Range Sig exceeds allowed range: %g - %g, %g - %g, trimming\n", dXLimitMin, dXSigMin, dXSigMax, dXLimitMax);
  dXSigMin = TMath::Max(dXSigMin, dXLimitMin);
  dXSigMax = TMath::Min(dXSigMax, dXLimitMax);
  if((dXBgInMin < dXLimitMin) || (dXBgInMax > dXLimitMax))
    printf("BinCounterObject::CheckRegions: Warning: Range BgIn exceeds allowed range: %g - %g, %g - %g, trimming\n", dXLimitMin, dXBgInMin, dXBgInMax, dXLimitMax);
  dXBgInMin = TMath::Max(dXBgInMin, dXLimitMin);
  dXBgInMax = TMath::Min(dXBgInMax, dXLimitMax);
  if((dXBgOutMin < dXLimitMin) || (dXBgOutMax > dXLimitMax))
    printf("BinCounterObject::CheckRegions: Warning: Range BgOut exceeds allowed range: %g - %g, %g - %g, trimming\n", dXLimitMin, dXBgOutMin, dXBgOutMax, dXLimitMax);
  dXBgOutMin = TMath::Max(dXBgOutMin, dXLimitMin);
  dXBgOutMax = TMath::Min(dXBgOutMax, dXLimitMax);
  if(bSideBands)
  {
    if(dXBgOutMin >= dXBgInMin)
    {
      printf("BinCounterObject::CheckRegions: Error: Wrong range of left side band: %g - %g\n", dXBgOutMin, dXBgInMin);
      bStatus = kFALSE;
    }
    if(dXBgInMax >= dXBgOutMax)
    {
      printf("BinCounterObject::CheckRegions: Error: Wrong range of right side band: %g - %g\n", dXBgInMax, dXBgOutMax);
      bStatus = kFALSE;
    }
  }
  if(dXSigMin >= dXSigMax)
  {
    printf("BinCounterObject::CheckRegions: Error: Wrong range of signal: %g - %g\n", dXSigMin, dXSigMax);
    bStatus = kFALSE;
  }
  if((dMean <= dXSigMin) || (dMean >= dXSigMax))
  {
    printf("BinCounterObject::CheckRegions: Error: Mean outside the signal range: %g - %g - %g\n", dXSigMin, dMean, dXSigMax);
    bStatus = kFALSE;
  }
  if(dSigma <= 0)
  {
    printf("BinCounterObject::CheckRegions: Error: Invalid value of sigma: %g\n", dSigma);
    bStatus = kFALSE;
  }
  if(dArea < 0)
  {
    printf("BinCounterObject::CheckRegions: Error: Invalid value of area: %g\n", dArea);
    bStatus = kFALSE;
  }
  if(bVerbose)
    printf("BinCounterObject::CheckRegions: End\n");
  return bStatus;
}

Bool_t BinCounterObject::SetRegions(Bool_t bSideBands)
{
  Double_t dEpsilon = dWidthBin * 1e-2;
  if(bVerbose)
    printf("BinCounterObject::SetRegions: Start\n");
  if(!bFixedRegions)
  {
    // set region edges
    dXSigMin = dMean - dNSigmaSig * dSigma; // signal min
    dXSigMax = dMean + dNSigmaSig * dSigma; // signal max
    dXBgInMin = dMean - dNSigmaBgIn * dSigma; // left side band max
    dXBgInMax = dMean + dNSigmaBgIn * dSigma; // right side band min
    dXBgOutMin = dMean - dNSigmaBgOut * dSigma; // left side band min
    dXBgOutMax = dMean + dNSigmaBgOut * dSigma; // right side band max
  }
  // check if new regions are OK
  if(!CheckRegions(bSideBands))
    return kFALSE;
  // correct values to match bin edges
  iBinSigMin = hisInput->FindBin(dXSigMin + dEpsilon);
  dXSigMin = hisInput->GetBinLowEdge(iBinSigMin);
  iBinSigMax = hisInput->FindBin(dXSigMax - dEpsilon);
  dXSigMax = hisInput->GetBinLowEdge(iBinSigMax + 1);
  iBinBgInMin = hisInput->FindBin(dXBgInMin - dEpsilon);
//  printf("Find: %g -> Bin: %d -> Low Edge: %g\n",dXBgInMin-dEpsilon,iBinBgInMin,hisInput->GetBinLowEdge(iBinBgInMin+1));
  dXBgInMin = hisInput->GetBinLowEdge(iBinBgInMin + 1);
  iBinBgInMax = hisInput->FindBin(dXBgInMax + dEpsilon);
//  printf("Find: %g -> Bin: %d -> Low Edge: %g\n",dXBgInMax+dEpsilon,iBinBgInMax,hisInput->GetBinLowEdge(iBinBgInMax));
  dXBgInMax = hisInput->GetBinLowEdge(iBinBgInMax);
  iBinBgOutMin = hisInput->FindBin(dXBgOutMin + dEpsilon);
  dXBgOutMin = hisInput->GetBinLowEdge(iBinBgOutMin);
  iBinBgOutMax = hisInput->FindBin(dXBgOutMax - dEpsilon);
  dXBgOutMax = hisInput->GetBinLowEdge(iBinBgOutMax + 1);
  printf("BinCounterObject::SetRegions: Side bands: %g (bin %d) - %g (bin %d), %g (bin %d) - %g (bin %d) (%g - %g sigma)\n", dXBgOutMin, iBinBgOutMin, dXBgInMin, iBinBgInMin, dXBgInMax, iBinBgInMax, dXBgOutMax, iBinBgOutMax, dNSigmaBgIn, dNSigmaBgOut);
  printf("BinCounterObject::SetRegions: Signal: %g (bin %d) - %g (bin %d) (%g sigma)\n", dXSigMin, iBinSigMin, dXSigMax, iBinSigMax, dNSigmaSig);
  // check if new regions are OK
  if(!CheckRegions(bSideBands))
    return kFALSE;
  // set regions in functions
  funFitGlob->SetRange(dXLimitMin, dXLimitMax);
  funFitBg->SetRange(dXBgOutMin, dXBgOutMax);
  funFitBgInteg->SetRange(dXSigMin, dXSigMax);
  funFitBg->FixParameter(4, dXBgInMin);
  funFitBg->FixParameter(5, dXBgInMax);
  Bool_t bStatus = kTRUE;
  if(bConsiderBg && bSideBands)
  {
    if(dNEntriesBg2Min < dNEntriesBg1Min)
      printf("BinCounterObject::SetRegions: Warning: Minimum for fitting pol2 is lower than minimum for fitting pol1\n");
//      iDegreePolInit = TMath::Min(iDegreePolInit,iDegreePolMax);
    iDegreePolInit = iDegreePolMax;
    Double_t dIntegralTmp = hisInput->Integral(iBinBgOutMin, iBinBgInMin) / (iBinBgInMin - iBinBgOutMin + 1); // average number of entries in the left side band
    if(iDegreePolInit > 0 && dIntegralTmp < dNEntriesBg1Min)
    {
      printf("BinCounterObject::SetRegions: Warning: Too few entries in the left side band for fitting with pol1,2,3: %g < %g, setting to pol0\n", dIntegralTmp, dNEntriesBg1Min);
      iDegreePolInit = 0;
    }
    if(iDegreePolInit > 1 && dIntegralTmp < dNEntriesBg2Min)
    {
      printf("BinCounterObject::SetRegions: Warning: Too few entries in the left side band for fitting with pol2,3: %g < %g, setting to pol1\n", dIntegralTmp, dNEntriesBg2Min);
      iDegreePolInit = 1;
    }
    dIntegralTmp = hisInput->Integral(iBinBgInMax, iBinBgOutMax) / (iBinBgOutMax - iBinBgInMax + 1); // average number of entries in the right side band
    if(iDegreePolInit > 0 && dIntegralTmp < dNEntriesBg1Min)
    {
      printf("BinCounterObject::SetRegions: Warning: Too few entries in the right side band for fitting with pol1,2,3: %g < %g, setting to pol0\n", dIntegralTmp, dNEntriesBg1Min);
      iDegreePolInit = 0;
    }
    if(iDegreePolInit > 1 && dIntegralTmp < dNEntriesBg2Min)
    {
      printf("BinCounterObject::SetRegions: Warning: Too few entries in the right side band for fitting with pol2,3: %g < %g, setting to pol1\n", dIntegralTmp, dNEntriesBg2Min);
      iDegreePolInit = 1;
    }
    /*
      dIntegralTmp = hisInput->Integral(iBinSigMin,iBinSigMax)/(iBinSigMax-iBinSigMin+1);
      if (dIntegralTmp<dNEntriesSigMin)
        {
          printf("BinCounterObject::SetRegions: Warning: Too few entries in the signal region: %g < %g\n",dIntegralTmp,dNEntriesSigMin);
          iDegreePolInit = 0;
        }
    */
  }
  if(!bStatus)
    return kFALSE;
  if(bVerbose)
    printf("BinCounterObject::SetRegions: End\n");
  return kTRUE;
}

Bool_t BinCounterObject::EstimateParameters(Int_t iDegPol)
{
  if(bVerbose)
    printf("BinCounterObject::EstimateParameters: Start\n");
  if(!SetRegions())
    return kFALSE;
  if(iDegPol == 0) // estimate constant
  {
    Double_t dIntegralSideBands = hisInput->Integral(iBinBgOutMin, iBinBgInMin) + hisInput->Integral(iBinBgInMax, iBinBgOutMax);
    Double_t dLengthSideBands = (dXBgInMin - dXBgOutMin) + (dXBgOutMax - dXBgInMax);
//  printf("Division 3 by %g\n",dLengthSideBands);
    Double_t dBgAvgSideBands = (bConsiderBg ? dIntegralSideBands / dLengthSideBands* dWidthBin : 0.);
//  printf("%d-%d + %d-%d: %g %g %g\n",iBinBgOutMin,iBinBgInMin,iBinBgInMax,iBinBgOutMax,dIntegralSideBands,dLengthSideBands,dBgAvgSideBands);
    dPar0 = dBgAvgSideBands;
    dPar1 = 0;
    dPar2 = 0;
    dPar3 = 0;
    Double_t dIntegralTotPeak = hisInput->Integral(iBinSigMin, iBinSigMax) * dWidthBin;
    Double_t dIntegralBgPeak = dBgAvgSideBands * (dXSigMax - dXSigMin);
    Double_t dIntegralSigPeak = dIntegralTotPeak - dIntegralBgPeak;
    dArea = TMath::Max(0., dIntegralSigPeak);
  }
  /*
    if (iDegPol==1) // estimate linear function
    {
      Double_t dIntegralLeft = hisInput->Integral(iBinBgOutMin,iBinBgInMin);
      Double_t dIntegralRight = hisInput->Integral(iBinBgInMax,iBinBgOutMax);
      Double_t dMeanLeft = 0;
      Double_t dMeanRight = 0;
      if (dIntegralLeft>0)
      {
      for (Int_t i=iBinBgOutMin; i<=iBinBgInMin; i++)
        dMeanLeft += hisInput->GetBinContent(i)*hisInput->GetXaxis()->GetBinCenter(i);
      dMeanLeft /= dIntegralLeft;
      }
      else
        dMeanLeft = (dXBgInMin-dXBgOutMin)/2.;
      if (dIntegralRight>0)
      {
      for (Int_t i=iBinBgInMax; i<=iBinBgOutMax; i++)
        dMeanRight += hisInput->GetBinContent(i)*hisInput->GetXaxis()->GetBinCenter(i);
      dMeanRight /= dIntegralRight;
      }
      else
        dMeanRight = (dXBgOutMax-dXBgInMax)/2.;
      dPar2 = 0;
      dPar3 = 0;
    }
  */
  funFitGlob->SetParameters(dPar0, dPar1, dPar2, dPar3, dArea, dMean, dSigma);
  funFitBg->SetParameters(dPar0, dPar1, dPar2, dPar3, dXBgInMin, dXBgInMax);
  funFitBgInteg->SetParameters(dPar0, dPar1, dPar2, dPar3);
  if(bVerbose)
    printf("BinCounterObject::EstimateParameters: End\n");
  return kTRUE;
}

Bool_t BinCounterObject::Run()
{
  // To do before calling this function:
  // Constructor
  // (SetXLimits)
  // (SetNSigmas)
  // (Plot) (initial conditions)
  // EstimateParameters
  // (Plot) (check the estimated parameters)

  if(bVerbose)printf("BinCounterObject::Run: Start\n");

  if(hisInput->Integral(hisInput->FindBin(dXLimitMin), hisInput->FindBin(dXLimitMax)) == 0.)
  {
    printf("BinCounterObject::Run: Empty histogram\n");
    dIntegralSignal = 0;
    dIntegralSignalErr = 0;
    dIntegralBg = 0;
    dIntegralBgErr = 0;
    return kTRUE;
  }

  Double_t dIntegralTmp = hisInput->Integral(iBinSigMin, iBinSigMax) / (iBinSigMax - iBinSigMin + 1);
  if(dIntegralTmp < dNEntriesSigMin)
  {
    printf("BinCounterObject::Run: Warning: Too few entries in the signal region for fitting: %g < %g, using fixed regions\n", dIntegralTmp, dNEntriesSigMin);
    bFixedRegions = kTRUE;
  }

  Int_t iDGlob = iDegreePolInit;
  if(!bConsiderBg)
    iDegreePolInit = -1;
//  printf("Consider bg %d\n",bConsiderBg);
//  printf("We start with pol %d\n",iDegreePolInit);
  if(!bFixedRegions)
  {
    printf("BinCounterObject::Run: Global fit: Start\n");
    if(!SetRegions())
      return kFALSE;
    // Fit bg+sig to estimate sigma
    printf("BinCounterObject::Run: Fitting background+signal in %g - %g\n", dXLimitMin, dXLimitMax);
    printf("BinCounterObject::Run: Consider bg %d\n", bConsiderBg);
    printf("BinCounterObject::Run: We start with pol %d\n", iDegreePolInit);
    for(iDGlob = iDegreePolInit; iDGlob >= -1; iDGlob--) // decrease degree of polynomial in case of failure and try again
    {
      RunFitGlobal(iDGlob, sOptionFitGlob.Data()); // pol-i
      if(IsFitOK(resultFitGlobal))
        break;
      if(bConsiderBg && iDGlob == 0)
        break;
    }
    if(!IsFitOK(resultFitGlobal))
    {
      printf("BinCounterObject::Run: Error: Failed to fit background+signal in histogram %s\n", sNameHis.Data());
      return kFALSE;
    }
    printf("BinCounterObject::Run: Info: Success to fit background+signal in histogram %s with Gauss+pol%d\n", sNameHis.Data(), iDGlob);

    dPar0 = funFitGlob->GetParameter(0);
    dPar1 = funFitGlob->GetParameter(1);
    dPar2 = funFitGlob->GetParameter(2);
    dPar3 = funFitGlob->GetParameter(3);
    dArea = funFitGlob->GetParameter(4);
    dMean = funFitGlob->GetParameter(5);
    dSigma = funFitGlob->GetParameter(6);
    dMeanErr = funFitGlob->GetParError(5);
    dSigmaErr = funFitGlob->GetParError(6);
    printf("BinCounterObject::Run: Global fit result: mean = %g +- %g, sigma = %g +- %g\n", dMean, dMeanErr, dSigma, dSigmaErr);
//  if (!CheckRegions())
//    return kFALSE;
  }
  else
  {
    printf("BinCounterObject::Run: No global fit\n");
  }

  // Set parameters for fitting side bands
  funFitBg->SetParameters(dPar0, dPar1, dPar2, dPar3, dXBgInMin, dXBgInMax);
  if(!SetRegions(bSubtract))
    return kFALSE;
//  Plot("canTmp");

  if(bSubtract)
  {
    dIntegralTmp = hisInput->Integral(iBinBgOutMin, iBinBgInMin);
    if(dIntegralTmp == 0)
    {
      printf("BinCounterObject::Run: Warning: No entries in the left side band: Turning off bg subtraction for histogram %s\n", sNameHis.Data());
      bSubtract = kFALSE;
    }
    dIntegralTmp = hisInput->Integral(iBinBgInMax, iBinBgOutMax);
    if(dIntegralTmp == 0)
    {
      printf("BinCounterObject::Run: Warning: No entries in the right side band: Turning off bg subtraction for histogram %s\n", sNameHis.Data());
      bSubtract = kFALSE;
    }
  }

  if(bSubtract)
  {
    // Fit bg in side bands and get the result info
    printf("BinCounterObject::Run: Fitting side bands in %g - %g, %g - %g\n", dXBgOutMin, dXBgInMin, dXBgInMax, dXBgOutMax);
    Int_t iDSB;
    for(iDSB = TMath::Min(iDGlob, iDegreePolInit); iDSB >= 0; iDSB--) // decrease degree of polynomial in case of failure and try again
    {
      RunFitSideBands(iDSB, sOptionFitSB); // pol-i
      if(!bFixedRegions)  // don't allow to decrease degree if the global fit succeeded with the current one
        break;
      if(IsFitOK(resultFitSideBands))  // uncomment if you want to allow decrease of degree when fixed regions
        break;
    }
    if(!IsFitOK(resultFitSideBands))
    {
      printf("BinCounterObject::Run: Error: Failed to fit background in side bands in histogram %s\n", sNameHis.Data());
      return kFALSE;
    }
    printf("BinCounterObject::Run: Info: Success to fit background in side bands in histogram %s with pol%d\n", sNameHis.Data(), iDSB);
  }
  else
  {
    printf("BinCounterObject::Run: No background subtraction\n");
  }

  Double_t dEpsilon = 1e-3;
  // get bg integral and error
  printf("BinCounterObject::Run: Extracting signal in %g - %g\n", dXSigMin, dXSigMax);
  dIntegralBg = 0;
  dIntegralBgErr = 0;
  funFitBgInteg->SetParameters(funFitBg->GetParameter(0), funFitBg->GetParameter(1), funFitBg->GetParameter(2), funFitBg->GetParameter(3), 0, 0);
  if(bSubtract)
  {
    TMatrixDSym cov = resultFitSideBands->GetCovarianceMatrix();
    PrintMatrix(cov);
//      printf("Division 4 by %g\n",dWidthBin);
    //dIntegralBg = funFitBgInteg->Integral(dXSigMin, dXSigMax, resultFitSideBands->GetParams()) / dWidthBin;
    dIntegralBg = funFitBgInteg->Integral(dXSigMin, dXSigMax) / dWidthBin;
    dIntegralBgErr = funFitBgInteg->IntegralError(dXSigMin, dXSigMax, resultFitSideBands->GetParams(), cov.GetMatrixArray()) / dWidthBin;
  }
  printf("BinCounterObject::Run: Bg: Integral and Err: %g +- %g (%g %%)\n", dIntegralBg, dIntegralBgErr, 100 * dIntegralBgErr / TMath::Abs(dIntegralBg + dEpsilon));
  printf("%s Bg: Integral and Err: %g +- %g (%g %%)\n", sNameHis.Data(), dIntegralBg, dIntegralBgErr, 100 * dIntegralBgErr / TMath::Abs(dIntegralBg + dEpsilon));
  if(dIntegralBg < 0 && bSubtract)
  {
    printf("BinCounterObject::Run: Error: Negative Bg integral: %g\n", dIntegralBg);
    return kFALSE;
  }
  // get bg+sig integral and error
  Double_t dIntegralTotErr;
  Double_t dIntegralTot = hisInput->IntegralAndError(iBinSigMin, iBinSigMax, dIntegralTotErr);
  printf("BinCounterObject::Run: All: Integral and Err: %g +- %g (%g %%)\n", dIntegralTot, dIntegralTotErr, 100 * dIntegralTotErr / TMath::Abs(dIntegralTot + dEpsilon));
  printf("%s All: Integral and Err: %g +- %g (%g %%)\n", sNameHis.Data(), dIntegralTot, dIntegralTotErr, 100 * dIntegralTotErr / TMath::Abs(dIntegralTot + dEpsilon));
  // get signal integral and error
  dIntegralSignal = dIntegralTot - dIntegralBg;
//  dIntegralBgErr = 0;
  dIntegralSignalErr = TMath::Sqrt(dIntegralTotErr * dIntegralTotErr + dIntegralBgErr * dIntegralBgErr);
  printf("BinCounterObject::Run: Sig: Integral and Err: %g +- %g (%g %%)\n", dIntegralSignal, dIntegralSignalErr, 100 * dIntegralSignalErr / TMath::Abs(dIntegralSignal + dEpsilon));
  printf("%s Sig: Integral and Err: %g +- %g (%g %%)\n", sNameHis.Data(), dIntegralSignal, dIntegralSignalErr, 100 * dIntegralSignalErr / TMath::Abs(dIntegralSignal + dEpsilon));
  // calculate sample purity
  if(dIntegralTot)
  {
    dPurity = dIntegralSignal / dIntegralTot;
    Double_t dErrTot = dIntegralTotErr * dIntegralBg / (dIntegralTot * dIntegralTot);
    Double_t dErrBg = dIntegralBgErr / dIntegralTot;
    dPurityErr = TMath::Sqrt(dErrTot * dErrTot + dErrBg * dErrBg);
  }
  if(bVerbose)
    printf("BinCounterObject::Run: End\n");
  return kTRUE;

  // To do after calling this function:
  // (Plot) (plot result of fitting)
  // GetSignalAndError
  // (GetMean, GetSigma, GetMeanError, GetSigmaError, GetBgAndError)
}

void BinCounterObject::RunFitGlobal(Int_t iDegreePol, TString sOption)
{
  if(bVerbose)
  {
    printf("BinCounterObject::RunFitGlobal: Start\n");
    printf("BinCounterObject::RunFitGlobal: Gauss+pol%d\n", iDegreePol);
  }
  funFitGlob->SetParameters(dPar0, dPar1, dPar2, dPar3, dArea, dMean, dSigma);
  for(Int_t i = 3; i > iDegreePol; i--)
  {
    printf("BinCounterObject::RunFitGlobal: Fixing parameter %d\n", i);
    funFitGlob->FixParameter(i, 0);
  }
  resultFitGlobal = hisInput->Fit(funFitGlob, sOption.Data());
  if(bVerbose)
  {
    PrintFitResults(resultFitGlobal, funFitGlob);
    printf("BinCounterObject::RunFitGlobal: End\n");
  }
}

void BinCounterObject::RunFitSideBands(Int_t iDegreePol, TString sOption)
{
  if(bVerbose)
  {
    printf("BinCounterObject::RunFitSideBands: Start\n");
    printf("BinCounterObject::RunFitSideBands: pol%d\n", iDegreePol);
  }
  funFitBg->SetParameters(dPar0, dPar1, dPar2, dPar3, dXBgInMin, dXBgInMax);
  for(Int_t i = 3; i > iDegreePol; i--)
  {
    printf("BinCounterObject::RunFitSideBands: Fixing parameter %d\n", i);
    funFitBg->FixParameter(i, 0);
  }
//  printf("BinCounterObject::RunFitSideBands: We go for fitting now\n");
  resultFitSideBands = hisInput->Fit(funFitBg, sOption.Data());
  if(bVerbose)
  {
    PrintFitResults(resultFitSideBands, funFitBg);
    printf("BinCounterObject::RunFitSideBands: End\n");
  }
}

void BinCounterObject::Plot(TString sNameFig, TString sTitle, Bool_t bLogY, TString sSuffix)
{
  if(bVerbose)
    printf("BinCounterObject::Plot: Start\n");
  if(!CheckHistogram())
    return;
  TCanvas* can = new TCanvas(Form("canTmp-%s", sNameHis.Data()), "", 800, 600);
  can->cd();
  if(sTitle.Length())
    hisInput->SetTitle(sTitle.Data());
  hisInput->Draw();
  if(!bFixedRegions)
  {
    funFitGlob->SetLineColor(kRed);
    funFitGlob->SetLineStyle(3);
    funFitGlob->DrawCopy("SAME");
  }
  if(bSubtract)
  {
    funFitBg->SetLineColor(kBlack);
    funFitBg->SetFillColor(kBlack);
    funFitBg->SetFillStyle(3013);
    funFitBg->DrawCopy("SAME");
  }
  funFitBgInteg->SetLineColor(kBlue);
  funFitBgInteg->SetLineStyle(2);
  funFitBgInteg->SetFillColor(kBlue);
  funFitBgInteg->SetFillStyle(3013);
  funFitBgInteg->DrawCopy("SAME");
  if(bLogY)
    can->SetLogy();
  can->SaveAs(Form("%s-%s.%s", sNameFig.Data(), sNameHis.Data(), sSuffix.Data()));
  delete can;
  if(bVerbose)
    printf("BinCounterObject::Plot: End\n");
  return;
}

void BinCounterObject::PrintFitResults(TFitResultPtr result, TF1* funFit)
{
//  if (!result)
//    {
//      printf("BinCounterObject::PrintFitResults: Error: No object TFitResultPtr given\n");
//      return;
//    }
  Double_t dEpsilon = 1e-4;
  Bool_t bIsValid = result->IsValid();
  Int_t iStatusFit = result->Status();
  Int_t iStatusCov = result->CovMatrixStatus();
  Double_t dChi2Result = result->Chi2();
  Int_t dNdfResult = result->Ndf();
  printf("BinCounterObject::PrintFitResults: IsValid = %d, Status = %d, CovMatrixStatus = %d, Chi2/NDF = %g/%d = %g\n", int(bIsValid), iStatusFit, iStatusCov, dChi2Result, dNdfResult, dChi2Result / (dNdfResult + dEpsilon));
  if(funFit)
  {
    Double_t dChi2Fun = funFit->GetChisquare();
    Int_t dNdfFun = funFit->GetNDF();
    printf("BinCounterObject::PrintFitResults: Fun: Chi2/NDF = %g/%d\n", dChi2Fun, dNdfFun);
  }
  return;
}

Bool_t BinCounterObject::IsFitOK(TFitResultPtr result)
{
//  if (!result)
//    {
//      printf("BinCounterObject::IsFitOK: Error: No object TFitResultPtr given\n");
//      return kFALSE;
//    }
  Bool_t bIsOK = kTRUE;
  Bool_t bIsValid = result->IsValid();
  Int_t iStatusFit = result->Status();
  Int_t iStatusCov = result->CovMatrixStatus();
  if(!bIsValid)
    bIsOK = kFALSE;
  if(iStatusFit != 0 && iStatusFit != 4000) // 0 = converged, 4000 = unchanged
    bIsOK = kFALSE;
  if(iStatusCov != 3) // 3 = accurate
    bIsOK = kFALSE;
//  if (bVerbose)
//    printf("BinCounterObject::IsFitOK: %s\n",(bIsOK ? "Yes" : "No"));
  return bIsOK;
}

Double_t BgFit(Double_t* x, Double_t* par)
{
  if(x[0] > par[4] && x[0] < par[5])
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0]; // a+b*x+c*x^2+d*x^3
}

void PrintMatrix(TMatrixDSym matrix)
{
  Double_t* elements = matrix.GetMatrixArray();
  Int_t iNElem = matrix.GetNoElements();
  Int_t iDim = Int_t(TMath::Sqrt(iNElem));
  if(iNElem != iDim * iDim)
  {
    printf("PrintMatrix: Error: Matrix is not squared\n");
    return;
  }
  printf("Matrix %d x %d:\n", iDim, iDim);
  for(Int_t iEl = 0; iEl < iNElem; iEl++)
  {
    printf("%g", elements[iEl]);
    if((iEl + 1) % iDim == 0)
      printf("\n");
    else
      printf("\t");
  }
  return;
}
