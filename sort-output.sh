#!/bin/bash

mkdir -p Event
mv canCentrality* Event

mkdir -p Counter
mv canCounter* Counter

mkdir -p Cent/BC
mv canInvMassCent* Cent
mv canBinCounter*-fh1V0InvMass*Cent_* Cent/BC

mkdir -p Cut
mv canInvMassCut* Cut

mkdir -p Pt-Inclusive/all
mkdir -p Pt-Inclusive/bg
mkdir -p Pt-Inclusive/BC
mv canPtInclusiveAll* Pt-Inclusive/all
mv canPtInclusiveBg* Pt-Inclusive/bg
mv canPtInclusiveMean* Pt-Inclusive
mv canPtInclusiveSigma* Pt-Inclusive
mv canPtInclusivePurity* Pt-Inclusive
mv canPtInclusiveSpectrum* Pt-Inclusive
mv canBinCounter*-PtInclusive* Pt-Inclusive/BC

mkdir -p PtEta-Inclusive/BC
mv canBinCounter*-PtEtaInclusive*-* PtEta-Inclusive/BC

mkdir -p Pt-In-Jets/all
mkdir -p Pt-In-Jets/bg
mkdir -p Pt-In-Jets/BC
mv canPtInJetsSpectrum* Pt-In-Jets
mv canPtInJCAll* Pt-In-Jets/all
mv canPtInJCBg* Pt-In-Jets/bg
mv canPtInJCSpectrum* Pt-In-Jets
mv canPtInPCAll* Pt-In-Jets/all
mv canPtInJCPurity* Pt-In-Jets
mv canPtInPCBg* Pt-In-Jets/bg
mv canPtInPCSpectrum* Pt-In-Jets
mv canPtInRCAll* Pt-In-Jets/all
mv canPtInRCBg* Pt-In-Jets/bg
mv canPtInRCSpectrum* Pt-In-Jets
mv canPtInFCAll* Pt-In-Jets/all
mv canPtInFCBg* Pt-In-Jets/bg
mv canPtInFCSpectrum* Pt-In-Jets
mv canBinCounter*-PtInRC*_* Pt-In-Jets/BC
mv canBinCounter*-PtInMCC*_* Pt-In-Jets/BC
mv canBinCounter*-PtInJets*_* Pt-In-Jets/BC
mv canBinCounter*-PtInPC*_* Pt-In-Jets/BC
mv canPtInBulkCompare*_* Pt-In-Jets
mv canPtInBulkSpectrum*_* Pt-In-Jets
mv canBulkOverJC*_* Pt-In-Jets
mv canBulkOverSignal*_* Pt-In-Jets
mv canBinCounter*-PtInFC*_* Pt-In-Jets/BC
mv canBinCounter*-PtInNJ*_* Pt-In-Jets/BC
mv canBinCounter*-PtInOC*_* Pt-In-Jets/BC

mkdir -p Jets
mv canPtJet* Jets
mv canPhiJet* Jets
mv canEtaJet* Jets
mv canNJet* Jets
mv canEtaPtJet* Jets
mv canPtTrack* Jets
mv canPtTrigger* Jets

mkdir -p EffPt/BC
mv canEffPtInclusive* EffPt
mv canBinCounter*-MCPtInclusive* EffPt/BC

mv canEffPtInJets* EffPt
mv canBinCounter*-MCPtInJets* EffPt/BC

mv canEffPtRatio* EffPt
mv canBinCounter*-MCPtInclusiveRebin* EffPt/BC

mkdir -p EffEtaPt/BC
mv canEffEtaPtInclusive* EffEtaPt
mv canBinCounter*-MCEtaPtInclusive* EffEtaPt/BC

mkdir -p EffPtEta/BC
mv canEffPtEtaIncl* EffPtEta
mv canBinCounter*-MCPtEtaIncl* EffPtEta/BC

mv canEffPtEtaInJets* EffPtEta
mv canBinCounter*-MCPtEtaInJet* EffPtEta/BC

mv canEffPtEtaRatio* EffPtEta
mv canGenPtEtaRatio* EffPtEta

mkdir -p EffEtaDaughters/EffPt
mkdir -p EffEtaDaughters/EffEta
mkdir -p EffEtaDaughters/EtaDaughters
mkdir -p EffEtaDaughters/PtDaughters
mkdir -p EffEtaDaughters/MeanEta
mv canEffDPtIn*_C* EffEtaDaughters/EffPt
mv canEffDPtJC*_C* EffEtaDaughters/EffPt
mv canEffDPtRatio*_C* EffEtaDaughters/EffPt
mv canEtaDIncl*Rec_C* EffEtaDaughters/EtaDaughters
mv canEtaDInJets*Rec_C* EffEtaDaughters/EtaDaughters
mv canEtaV0Incl*Eff_C* EffEtaDaughters/EffEta
mv canEtaV0Incl*Rec_C* EffEtaDaughters/EffEta
mv canEtaV0Incl*Gen_C* EffEtaDaughters/EffEta
mv canEtaV0InJets*Eff_C* EffEtaDaughters/EffEta
mv canEtaV0InJets*Rec_C* EffEtaDaughters/EffEta
mv canEtaV0InJets*Gen_C* EffEtaDaughters/EffEta
mv canEtaV0Ratio*Eff_C* EffEtaDaughters/EffEta
mv canDPt*_C* EffEtaDaughters/PtDaughters
mv canPtRatio*_C* EffEtaDaughters/PtDaughters
mv canMeanEta* EffEtaDaughters/MeanEta

mkdir -p Correlations/BC
mv canBinCounter2-Correl* Correlations/BC

mkdir -p Ratios
mv canRatio* Ratios

