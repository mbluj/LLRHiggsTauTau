```
scram project -n CMSSW_9_4_8_ntuple CMSSW CMSSW_9_4_8
cd CMSSW_9_4_8_ntuple/src
cmsenv

git cms-init
git cms-merge-topic mbluj:CMSSW_9_4_X_DPFIso
 
git clone https://github.com/SVfit/ClassicSVfit TauAnalysis/ClassicSVfit -b release_2018Mar20
git clone https://github.com/SVfit/SVfitTF TauAnalysis/SVfitTF
git clone https://github.com/akalinow/LLRHiggsTauTau.git
git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections 

cd LLRHiggsTauTau; git checkout Run2017; cd -

scram b -j 4

# Clone DNN datafiles to cmssw externals (it must be done after 1st compilation)
mkdir -p $CMSSW_BASE/external/$SCRAM_ARCH/data
cd $CMSSW_BASE/external/$SCRAM_ARCH/data
git clone https://github.com/cms-tau-pog/RecoTauTag-TrainingFiles -b master RecoTauTag/TrainingFiles/data
cd -

cd LLRHiggsTauTau/NtupleProducer/test/
git clone https://github.com/akalinow/Production.git
cd Production; git checkout Run2017;
```

Edit production scripts and you are done




