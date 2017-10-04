### Description

This package provides tools to READ miniAOD events and store the results in an TNtuple optimized for the HTauTau analyses.
The package workflow is the following:
Main flags/cut/configuration tools can be set in analyzer.py
Modules and plugins are called by python/HiggsTauTauProducer.py 
	Cuts on leptons/pairs are also set here
	This module creates AF/OS and AF/SS pairs, runs the SVfit on the pairs and store the useful variables in a TNtuple
The stored variables are set in plugins/HTauTauNtupleMaker.cc

### Instructions for 8_0_6 (miniAOD 2016)

```
cmsrel CMSSW_8_0_6
cd CMSSW_8_0_6/src
cmsenv
# Z-recoil corrections
git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections
git clone https://github.com/LLRCMS/LLRHiggsTauTau
cd LLRHiggsTauTau; git checkout master
cd -
git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
cd Muon/MuonAnalysisTools ; git checkout master -- interface/MuonEffectiveArea.h
cd -
git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
cd EGamma/EGammaAnalysisTools; git checkout c0db796 -- interface/ElectronEffectiveArea.h
cd -
# FSR corrections
git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
cd UFHZZAnalysisRun2 ; git checkout master FSRPhotons
cd -
# SVfit
git clone https://github.com/veelken/SVfit_standalone TauAnalysis/SVfitStandalone
cd TauAnalysis/SVfitStandalone
git checkout HIG-16-006
cd -
scram b -j 4
```

### Quick usage:
Define the files you want to run in analyzer.py and run cmsRun analyzer.py
Please note that since 7_4_7 we switched to a new eleID recipe and we are not using anymore git cms-merge-topic sregnard:Phys14ElectronMvaIdFor745

