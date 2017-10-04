### Description

This package provides tools to READ miniAOD events and store the results in an TNtuple optimized for the HTauTau analyses.
The package workflow is the following:
Main flags/cut/configuration tools can be set in analyzer.py
Modules and plugins are called by python/HiggsTauTauProducer.py 
	Cuts on leptons/pairs are also set here
	This module creates AF/OS and AF/SS pairs, runs the SVfit on the pairs and store the useful variables in a TNtuple
The stored variables are set in plugins/HTauTauNtupleMaker.cc

### Instructions for 9_2_9 (Summer 2017)

```
scram project -n CMSSW_9_2_9_patch1_ntuple CMSSW CMSSW_9_2_9_patch1
cd CMSSW_9_2_9_patch1_ntuple/src
cmsenv

ssh lxplus.cern.ch "(cd /afs/cern.ch/user/m/mbluj/public/ProductionSummer2017; tar -cf - . )" | tar -xvf -

tar -xzvf EGamma_EGammaAnalysisTools.tar.gz
tar -xzvf HTT-utilities_RecoilCorrections.tar.gz
tar -xzvf Muon_MuonAnalysisTools.tar.gz
tar -xzvf TauAnalysis_SVfitStandalone.tar.gz
tar -xzvf UFHZZAnalysisRun2_FSRPhotons.tar.gz

rm *.tar.gz
 
git clone https://github.com/akalinow/LLRHiggsTauTau.git
cd LLRHiggsTauTau; git checkout Run2017; cd -

scram b -j 4

```

### Quick usage:
Define the files you want to run in analyzer.py and run cmsRun analyzer.py
Please note that since 7_4_7 we switched to a new eleID recipe and we are not using anymore git cms-merge-topic sregnard:Phys14ElectronMvaIdFor745

