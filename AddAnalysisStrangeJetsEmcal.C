#include <vector>
#include "TROOT.h"
#include "TSystem.h"
#include "AuxFunctions.h"

Int_t AddAnalysisStrangeJetsEmcal(
  Bool_t bIsPbPb = 1, // switch between Pb+Pb and p+p
  TString sYear = "2010", // set the run period (used on grid)
  Bool_t bMCData = 0 // 0, 1 - analysis of Monte Carlo data?
)
{
  // Event selection
  Double_t dVertexWindow = 10; // [cm] range of z of the primary vertex
  Double_t dVertexR2 = 1; // [cm^2] radius of the primary vertex
  Double_t dCentLo = 0; // minimum centrality
  Double_t dCentUp = 10; // maximum centrality
  Double_t dDeltaZMax = 0.5; // [cm] max |delta z| between nominal prim vtx and SPD vertex
  Int_t iNContribMin = 1; // minimum number of prim vtx contributors
  Bool_t bMultiplicity = 1; // use AliMultSelection to estimate centrality
  Bool_t bIonut = 1; // use Ionut's cut

  // MC selection
  TString sGenerator = ""; // "Hijing"; // MC generator name

  // Configuration of pools for event mixing
  Bool_t bCorrelations = 1; // switch for V0-jet correlations
  Int_t iSizePool = 1; // available number of events per pool
  Int_t iNJetsPerPool = 100; // required number of jets available in each pool
  Float_t fFractionMin = 1.; // minimum fraction of iNJetsPerPool at which pool is ready (default: 1.0)
  Int_t iNEventsMin = 1; // if non-zero: number of filled events after which pool is ready regardless of iNJetsPerPool (default: 0)
  Double_t dDeltaEtaMax = 0.2; // maximum delta-eta_V0-jet for angular correlations

  // Data selection
  UInt_t iPhysicsSelectionFlag = 0; // set by physics selection and passed to the task, kMB, kUserDefined etc
  TString sPathOADB = "";
  TString sPathOADB2011 = "$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/data/OADB-LHC11h.root";
  TString sPathOADBMC = "";
  TString sPathOADBMC2011 = "LHC12a17d_fix";
  const char* runPeriod = "";
  if(bIsPbPb) // Pb+Pb
  {
    if(sYear == "2010") // Pb+Pb 2010
    {
      iPhysicsSelectionFlag = AliVEvent::kMB;
      runPeriod = "LHC10h";
      dDeltaZMax = 0;
    }
    else if(sYear == "2011") // Pb+Pb 2011
    {
      iPhysicsSelectionFlag = (AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
      runPeriod = "LHC11h";
      sPathOADB = sPathOADB2011;
      if (bMCData)
        sPathOADBMC = sPathOADBMC2011;
    }
  }
  else // p+p 2010
  {
    iPhysicsSelectionFlag = AliVEvent::kMB;
    runPeriod = "LHC10d";
    dDeltaZMax = 0;
  }

  // Jet analysis & tracks
//  Double_t dRadiusJet[] = {0.2, 0.3}; // resolution parameter R of clustering jet finder
  Double_t dRadiusJet[] = {0.2, 0.3}; // resolution parameter R of clustering jet finder
  const Int_t iNRadiusJet = sizeof(dRadiusJet) / sizeof(dRadiusJet[0]);
  Double_t dDistanceV0JetMax[] = {0.2, 0.3}; // maximum distance between V0 and jet axis used for finding V0s in the jet cone
  const Int_t iNDistanceV0JetMax = sizeof(dDistanceV0JetMax) / sizeof(dDistanceV0JetMax[0]);
  Double_t dRadiusJetBg = 0.4; // resolution parameter R of clustering jet finder for background estimation
  Double_t dPtTrackJetMin = 0.15; // [GeV/c] minimum pT of charged jet constituents

  // V0 selection
  // Daughter tracks
  Bool_t bTPCRefit = 1; // TPC refit for daughter tracks
  Bool_t bRejectKinks = 1; // reject kink-like production vertices of daughter tracks
  Bool_t bFindableClusters = 0; // require positive number of findable clusters
  Double_t dCutNCrossedRowsTPCMin = 0.; // min number of crossed TPC rows
  Double_t dCutCrossedRowsOverFindMin = 0.; // min ratio crossed rows / findable clusters
  Double_t dCutCrossedRowsOverFindMax = 0.; // max ratio crossed rows / findable clusters
  Double_t dCutPtDaughterMin = 0.; // [GeV/c] min transverse momentum of daughter tracks
  Double_t dCutDCAToPrimVtxMin = 0.1; // [cm] min DCA of daughters to the prim vtx
  Double_t dCutDCADaughtersMax = 1.; // [sigma of TPC tracking] max DCA between daughters
  Double_t dCutEtaDaughterMax = 0.8; // max |eta| of daughter tracks
  Double_t dCutNSigmadEdxMax = 0.; // [sigma dE/dx] max difference between measured and expected signal of dE/dx in the TPC
  Double_t dPtProtonPIDMax = 0.; // [GeV/c] maxium pT of proton for applying PID cut
  // V0 candidate
  Bool_t bOnFly = 0; // on-the-fly (yes) or offline (no) reconstructed
  Double_t dCutCPAKMin = 0.998; // min cosine of the pointing angle, K0S
  Double_t dCutCPALMin = 0.998; // min cosine of the pointing angle, Lambda
  Double_t dCutRadiusDecayMin = 5.; // [cm] min radial distance of the decay vertex
  Double_t dCutRadiusDecayMax = 100.; // [cm] max radial distance of the decay vertex
  Double_t dCutEtaV0Max = 0.7; // max |eta| of V0
  Double_t dCutRapV0Max = 0.; // max |rapidity| of V0
  Double_t dCutNTauKMax = 5.0; // [tau] max proper lifetime in multiples of the mean lifetime, K0S
  Double_t dCutNTauLMax = 5.0; // [tau] max proper lifetime in multiples of the mean lifetime, Lambda
  Bool_t bCutArmPod = 1; // Armenteros-Podolanski for K0S
  Bool_t bCutCross = 0; // cross contamination

  if(!bIsPbPb) // pp cuts
  {
    bTPCRefit = 1;
    bRejectKinks = 1;
    bFindableClusters = 1;
    dCutNCrossedRowsTPCMin = 70.;
    dCutCrossedRowsOverFindMin = 0.8;
    dCutCrossedRowsOverFindMax = 0.;
    dCutPtDaughterMin = 0.;
    dCutDCAToPrimVtxMin = 0.06;
    dCutDCADaughtersMax = 1.;
    dCutEtaDaughterMax = 0.8;
    dCutNSigmadEdxMax = 5.;
    dPtProtonPIDMax = 0.;
    bOnFly = 0;
    dCutCPAKMin = 0.97;
    dCutCPALMin = 0.995;
    dCutRadiusDecayMin = 0.5;
    dCutRadiusDecayMax = 1000.;
    dCutEtaV0Max = 0.;
    dCutRapV0Max = 0.5;
    dCutNTauKMax = 7.450; // 20 cm
    dCutNTauLMax = 3.802; // 30 cm
    bCutArmPod = 1;
    bCutCross = 1;
  }
  // jet selection for V0 analysis
  Int_t iBgSubtraction = 1; // subtraction of rho from jet pt, 0 - no subtraction, 1 - scalar subtraction, 2 - vector subtraction
  Double_t dCutPtJetMin = 5; // [GeV/c] minimum jet pt for analysis of V0s in jets
  Double_t dCutPtTrackJetMin = 5; // [GeV/c] minimum pt of leading jet-track for analysis of V0s in jets
  Double_t dCutAreaPercJetMin = 0.6; // [pi*R^2] minimum jet area with respect to the expected value

  // Active tasks
  Bool_t bEventPlane = 0; // Event plane calculation
  Bool_t bV0Analysis = 1; // V0 analysis
  Bool_t bV0JetAnalysis = 1; // Jet reconstruction: KT, ANTIKT, analysis of V0s in jets
  Bool_t bCompareTriggers = 0; // pt correlations of jets with trigger tracks
  Bool_t bCutVariations = 0; // run task with cut variations

  // Output files
  TString kResultsFileName = "FinalOutputs.root";

  Int_t iCommonTaskDebugLevel = 1;
  // >0 - basic code structure, non-critical warnings
  // >1 - report success
  // >2 - basic details, counts of objects
  // >3 - more details
  // >4 - very detailed reports about frequent steps

  // Some pre-settings and constants
  enum AlgoType {kKT, kANTIKT};
  enum JetType  {kFULLJETS, kCHARGEDJETS, kNEUTRALJETS};

  // ################# Now: Add some basic tasks
  // =============================================================
  // Configure analysis tasks and add them to the analysis manager

  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr)
  {
    printf("AddAnalysisStrangeJetsEmcal: Error: No analysis manager\n");
    return -1;
  }

  // Check type of input and create handler for it
  AliAODInputHandler* aodH = reinterpret_cast<AliAODInputHandler*>(gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C"))));

  printf("AddAnalysisStrangeJetsEmcal: Configuring analysis tasks\n");

  //PID Response Task
  if(dCutNSigmadEdxMax > 0. && (!bIsPbPb || (bIsPbPb && dPtProtonPIDMax > 0.)))
  {
    AliAnalysisTaskPIDResponse* taskPid = reinterpret_cast<AliAnalysisTaskPIDResponse*>(gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"))));
    if(!taskPid)
    {
      printf("AddAnalysisStrangeJetsEmcal: Failed to add PID response task\n");
      return 1;
    }
//    taskPid->SetUserDataRecoPass(2);
  }

  // Event plane calculation
  if(bEventPlane)
  {
    AliEPSelectionTask* taskEP = reinterpret_cast<AliEPSelectionTask*>(gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C"))));
    if(!taskEP)
      ::Warning("AddAnalysisStrangeJetsEmcal", "AliEventplans cannot run for this train conditions - EXCLUDED");
  }

  if(bMultiplicity)
  {
    AliMultSelectionTask* task = reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ProcessLine(Form(".x %s(%d)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"), kFALSE)));
    task->SetSelectedTriggerClass(iPhysicsSelectionFlag);
    task->SetAddInfo(kTRUE);
    if(sPathOADB.Length())
      task->SetAlternateOADBFullManualBypass(sPathOADB.Data());
    if(sPathOADBMC.Length())
      task->SetAlternateOADBFullManualBypassMC(sPathOADBMC.Data());
  }

  // Names of the different objects passed around; these are the default names; added here mostly for documentation purposes
  TString inputTracks = "AODFilterTracks";
  TString tracksName = "usedefault";
  TString particlesMCName = ""; //"MCParticlesSelected";
  TString clustersName = ""; //"EmcCaloClusters";
  TString clustersCorrName = ""; //"CaloClustersCorr";
  TString rhoName = "";
  TString sType = "TPC";
  Double_t dClusPtCut = 0.;
  Double_t dGhostArea = 0.005;
  AliEmcalJetTask* jetFinderTaskSignal = 0;
  std::vector<AliEmcalJetTask*> vecJetFinder;

  AliEmcalJetTask* jetFinderTaskBg = nullptr;
  if(bV0JetAnalysis)
  {
    // Add jet finders
    TString sRunPeriod = "lhc11h";
    AliTrackContainer::SetDefTrackCutsPeriod(sRunPeriod);
    for(Int_t iR = 0; iR < iNRadiusJet; iR++)
    {
      jetFinderTaskSignal = AliEmcalJetTask::AddTaskEmcalJet(tracksName, clustersCorrName, AliJetContainer::antikt_algorithm, dRadiusJet[iR], AliJetContainer::kChargedJet, dPtTrackJetMin, dClusPtCut, dGhostArea, AliJetContainer::pt_scheme, "Jet", 0., kFALSE, kFALSE);
      jetFinderTaskSignal->SelectCollisionCandidates(iPhysicsSelectionFlag);
      if(!jetFinderTaskSignal)
        return 1;
      vecJetFinder.push_back(jetFinderTaskSignal);
    }

    jetFinderTaskBg = AliEmcalJetTask::AddTaskEmcalJet(tracksName, clustersCorrName, AliJetContainer::kt_algorithm, dRadiusJetBg, AliJetContainer::kChargedJet, dPtTrackJetMin, dClusPtCut, dGhostArea, AliJetContainer::pt_scheme, "Jet", 0., kFALSE, kFALSE);
    jetFinderTaskBg->SelectCollisionCandidates(iPhysicsSelectionFlag);

    // Background density calculation
    if(bIsPbPb)
    {
      rhoName = "Rho";
      AliAnalysisTaskRho* pRhoTask = reinterpret_cast<AliAnalysisTaskRho*>(gInterpreter->ProcessLine(Form(".x %s(\"%s\", \"%s\", \"%s\", %f)", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoNew.C"), tracksName.Data(), clustersCorrName.Data(), rhoName.Data(), dRadiusJetBg)));
      pRhoTask->SetExcludeLeadJets(2);
      pRhoTask->SelectCollisionCandidates(iPhysicsSelectionFlag);

    }
  }

  // Here you can put in your AddTaskMacro for your task

  // Analysis of V0s (K0S, Lambda, anti-Lambda), inclusive, in jets
  TString sFlag = "";
  if(bV0Analysis)
  {
    printf("AddAnalysisStrangeJetsEmcal: V0 analysis\n");
    AliAnalysisTaskV0sInJetsEmcal* taskV0 = 0;
    for(Int_t iR = 0; iR < iNRadiusJet; iR++)
    {
      if(!bV0JetAnalysis && iR > 0) // don't create more than one task if no jet analysis
        break;
      if(bCutVariations)
        sFlag = "Default";
//      printf("AddAnalysisStrangeJetsEmcal: Step 1, %d\n",iR);
      if(bV0JetAnalysis)
        taskV0 = reinterpret_cast<AliAnalysisTaskV0sInJetsEmcal*>(gInterpreter->ProcessLine(Form(".x %s(\"%s\", %f, \"%s\", %f, \"%s\", %d, \"%s\", \"%s\", \"%s\", \"%s\", \"%s\")", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskV0sInJetsEmcal.C"),
        vecJetFinder[iR]->GetName(), dRadiusJet[iR], jetFinderTaskBg->GetName(), dRadiusJetBg, kResultsFileName.Data(), bMCData, tracksName.Data(), "", rhoName.Data(), sType.Data(), sFlag.Data())));
      else
        taskV0 = reinterpret_cast<AliAnalysisTaskV0sInJetsEmcal*>(gInterpreter->ProcessLine(Form(".x %s(\"%s\", %f, \"%s\", %f, \"%s\", %d, \"%s\", \"%s\", \"%s\", \"%s\", \"%s\")", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskV0sInJetsEmcal.C"),
        "", 0., "", 0., kResultsFileName.Data(), bMCData, "", "", "", "", sFlag.Data())));
      if(!taskV0)
        return 3;

//      printf("AddAnalysisStrangeJetsEmcal: Step 2\n");
      taskV0->SetDebugLevel(iCommonTaskDebugLevel);
      // data selection
      taskV0->SetIsPbPb(bIsPbPb);
      if(bMCData && sGenerator.Length())
        taskV0->SetGeneratorName(sGenerator.Data());
      // event selection
      taskV0->SelectCollisionCandidates(iPhysicsSelectionFlag);
      taskV0->SetEventCuts(dVertexWindow, dVertexR2, dCentLo, dCentUp, dDeltaZMax, iNContribMin);
      taskV0->SetUseMultiplicity(bMultiplicity);
      taskV0->SetUseIonutCut(bIonut);
      // V0 selection
      taskV0->SetCutTPCRefit(bTPCRefit);
      taskV0->SetCutRejectKinks(bRejectKinks);
      taskV0->SetCutFindableClusters(bFindableClusters);
      taskV0->SetCutDCAToPrimVtxMin(dCutDCAToPrimVtxMin);
      taskV0->SetCutDCADaughtersMax(dCutDCADaughtersMax);
      taskV0->SetCutEtaDaughterMax(dCutEtaDaughterMax);
      taskV0->SetCutNSigmadEdxMax(dCutNSigmadEdxMax);
      taskV0->SetPtProtonPIDMax(dPtProtonPIDMax);
      taskV0->SetOnFly(bOnFly);
      taskV0->SetCutCPAKMin(dCutCPAKMin);
      taskV0->SetCutCPALMin(dCutCPALMin);
      taskV0->SetCutRadiusDecayMin(dCutRadiusDecayMin);
      taskV0->SetCutRadiusDecayMax(dCutRadiusDecayMax);
      taskV0->SetCutEtaV0Max(dCutEtaV0Max);
      taskV0->SetCutNTauKMax(dCutNTauKMax);
      taskV0->SetCutNTauLMax(dCutNTauLMax);
      taskV0->SetCutNCrossedRowsTPCMin(dCutNCrossedRowsTPCMin);
      taskV0->SetCutCrossedRowsOverFindMin(dCutCrossedRowsOverFindMin);
      taskV0->SetCutCrossedRowsOverFindMax(dCutCrossedRowsOverFindMax);
      taskV0->SetCutPtDaughterMin(dCutPtDaughterMin);
      taskV0->SetCutRapV0Max(dCutRapV0Max);
      taskV0->SetCutArmPod(bCutArmPod);
      taskV0->SetCutCross(bCutCross);
      // jet selection
      if(bV0JetAnalysis)
      {
        taskV0->SetJetSelection(1);
        taskV0->SetBgSubtraction(iBgSubtraction);
        taskV0->SetPtJetMin(dCutPtJetMin);
        taskV0->SetPtTrackJetMin(dCutPtTrackJetMin);
        taskV0->SetAreaPercJetMin(dCutAreaPercJetMin);
        taskV0->SetDistanceV0JetMax(dDistanceV0JetMax[iR]);
        taskV0->SetCompareTriggerTracks(bCompareTriggers);
        // correlations
        if(bCorrelations)
        {
          taskV0->SetCorrelations(1);
          taskV0->SetPoolParam(iSizePool, iNJetsPerPool, fFractionMin, iNEventsMin);
          taskV0->SetDeltaEtaMax(dDeltaEtaMax);
        }
      }
      else
        taskV0->SetJetSelection(0);

      if(bCutVariations)
      {
        TString sCutName = "";
        Double_t dCutVal = 0;
        cut kCut = kDCAV;

        for(Int_t iCut = 0; iCut < iNArrayVarSys; iCut++)
        {
          sCutName = sArrayVarSysName[iCut];
          kCut = kArrayVarSysType[iCut];

          printf("Variations for %s\n", sCutName.Data());
          for(Int_t iVar = 0; iVar < iArrayNVarSys[iCut]; iVar++)
          {
            dCutVal = dArrayVarSysK[iCut][iVar];
            printf("Variation for %s: %g\n", sCutName.Data(), dCutVal);
            sFlag = Form("%s%d", sCutName.Data(), iVar + 1);

            if(bV0JetAnalysis)
              taskV0 = reinterpret_cast<AliAnalysisTaskV0sInJetsEmcal*>(gInterpreter->ProcessLine(Form(".x %s(\"%s\", %f, \"%s\", %f, \"%s\", %d, \"%s\", \"%s\", \"%s\", \"%s\", \"%s\")", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskV0sInJetsEmcal.C"),
              vecJetFinder[iR]->GetName(), dRadiusJet[iR], jetFinderTaskBg->GetName(), dRadiusJetBg, kResultsFileName.Data(), bMCData, tracksName.Data(), "", rhoName.Data(), sType.Data(), sFlag.Data())));
            else
              taskV0 = reinterpret_cast<AliAnalysisTaskV0sInJetsEmcal*>(gInterpreter->ProcessLine(Form(".x %s(\"%s\", %f, \"%s\", %f, \"%s\", %d, \"%s\", \"%s\", \"%s\", \"%s\", \"%s\")", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskV0sInJetsEmcal.C"),
              "", 0., "", 0., kResultsFileName.Data(), bMCData, "", "", "", "", sFlag.Data())));
            if(!taskV0)
              return 4;

            taskV0->SetDebugLevel(iCommonTaskDebugLevel);
            // data selection
            taskV0->SetIsPbPb(bIsPbPb);
            if(bMCData && sGenerator.Length())
              taskV0->SetGeneratorName(sGenerator.Data());
            // event selection
            taskV0->SelectCollisionCandidates(iPhysicsSelectionFlag);
            taskV0->SetEventCuts(dVertexWindow, dVertexR2, dCentLo, dCentUp, dDeltaZMax, iNContribMin);
            taskV0->SetUseMultiplicity(bMultiplicity);
            taskV0->SetUseIonutCut(bIonut);
            // V0 selection
            taskV0->SetCutTPCRefit(bTPCRefit);
            taskV0->SetCutRejectKinks(bRejectKinks);
            taskV0->SetCutFindableClusters(bFindableClusters);
            taskV0->SetCutDCAToPrimVtxMin(dCutDCAToPrimVtxMin);
            taskV0->SetCutDCADaughtersMax(dCutDCADaughtersMax);
            taskV0->SetCutEtaDaughterMax(dCutEtaDaughterMax);
            taskV0->SetCutNSigmadEdxMax(dCutNSigmadEdxMax);
            taskV0->SetPtProtonPIDMax(dPtProtonPIDMax);
            taskV0->SetOnFly(bOnFly);
            taskV0->SetCutCPAKMin(dCutCPAKMin);
            taskV0->SetCutCPALMin(dCutCPALMin);
            taskV0->SetCutRadiusDecayMin(dCutRadiusDecayMin);
            taskV0->SetCutRadiusDecayMax(dCutRadiusDecayMax);
            taskV0->SetCutEtaV0Max(dCutEtaV0Max);
            taskV0->SetCutNTauKMax(dCutNTauKMax);
            taskV0->SetCutNTauLMax(dCutNTauLMax);
            taskV0->SetCutNCrossedRowsTPCMin(dCutNCrossedRowsTPCMin);
            taskV0->SetCutCrossedRowsOverFindMin(dCutCrossedRowsOverFindMin);
            taskV0->SetCutCrossedRowsOverFindMax(dCutCrossedRowsOverFindMax);
            taskV0->SetCutPtDaughterMin(dCutPtDaughterMin);
            taskV0->SetCutRapV0Max(dCutRapV0Max);
            taskV0->SetCutArmPod(bCutArmPod);
            taskV0->SetCutCross(bCutCross);
            // jet selection
            if(bV0JetAnalysis)
            {
              taskV0->SetJetSelection(1);
              taskV0->SetBgSubtraction(iBgSubtraction);
              taskV0->SetPtJetMin(dCutPtJetMin);
              taskV0->SetPtTrackJetMin(dCutPtTrackJetMin);
              taskV0->SetAreaPercJetMin(dCutAreaPercJetMin);
              taskV0->SetDistanceV0JetMax(dDistanceV0JetMax[iR]);
              taskV0->SetCompareTriggerTracks(bCompareTriggers);
              // correlations
              if(bCorrelations)
              {
                taskV0->SetCorrelations(1);
                taskV0->SetPoolParam(iSizePool, iNJetsPerPool, fFractionMin, iNEventsMin);
                taskV0->SetDeltaEtaMax(dDeltaEtaMax);
              }
            }
            else
              taskV0->SetJetSelection(0);

            switch(kCut)
            {
              case kDCAV:
                taskV0->SetCutDCAToPrimVtxMin(dCutVal);
                break;
              case kDCAD:
                taskV0->SetCutDCADaughtersMax(dCutVal);
                break;
              case kCPA:
                taskV0->SetCutCPAKMin(dCutVal);
                taskV0->SetCutCPALMin(dCutVal);
                break;
              case kTau:
                taskV0->SetCutNTauKMax(dCutVal);
                taskV0->SetCutNTauLMax(dCutVal);
                break;
              case kRMin:
                taskV0->SetCutRadiusDecayMin(dCutVal);
                break;
              case kRMax:
                taskV0->SetCutRadiusDecayMax(dCutVal);
                break;
              default:
                printf("Unknown cut\n");
                return 5;
                break;
            }
          }
        }
      }
    }
  }

  return 0;
}
