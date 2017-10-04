import FWCore.ParameterSet.Config as cms

## TRIGGER VERSIONING: insert paths without the explicit version number, e.g. HLT_anything_v 
## do not remove any other part of path name or the check will give meaningless results!

## leg numbering: write the object type the two legs refer to
## 11: electrons, 13: muons, 15: taus, 999: others/unused (e.g. leg 2 in single muon)
##

TRIGGERLIST=[]
#list triggers and filter paths here!
# channel: kemu=0, ketau=1,kmutau=2,ktautau=3
HLTLIST = cms.VPSet(

### === Single muon triggers
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
     cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ), 
     cms.PSet (
        HLT = cms.string("HLT_IsoMu27_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ), 
### === mu tauh triggers
   cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),   
   cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20MediumChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltOverlapFilterIsoMu20MediumChargedIsoPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),   
   cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20TightChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltOverlapFilterIsoMu20TightChargedIsoPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),   
   cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20TightChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path2 = cms.vstring ("hltOverlapFilterIsoMu20LooseChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),   
   cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20MediumChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path2 = cms.vstring ("hltOverlapFilterIsoMu20MediumChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),   
   cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20TightChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path2 = cms.vstring ("hltOverlapFilterIsoMu20TightChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),   
#
   cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu24LooseChargedIsoPFTau20"),
        path2 = cms.vstring ("hltOverlapFilterIsoMu24LooseChargedIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),   
   cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu24MediumChargedIsoPFTau20"),
        path2 = cms.vstring ("hltOverlapFilterIsoMu24MediumChargedIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),   
   cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu24TightChargedIsoPFTau20"),
        path2 = cms.vstring ("hltOverlapFilterIsoMu24TightChargedIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),   
   cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu24LooseChargedIsoTightOOSCPhotonsPFTau20"),
        path2 = cms.vstring ("hltOverlapFilterIsoMu24LooseChargedIsoTightOOSCPhotonsPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),   
   cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu24MediumChargedIsoTightOOSCPhotonsPFTau20"),
        path2 = cms.vstring ("hltOverlapFilterIsoMu24MediumChargedIsoTightOOSCPhotonsPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),   
   cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu24TightChargedIsoTightOOSCPhotonsPFTau20"),
        path2 = cms.vstring ("hltOverlapFilterIsoMu24TightChargedIsoTightOOSCPhotonsPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),   
### === single tauh triggers
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v"),
        path1 = cms.vstring ("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched1Prong"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
        ),
### === tauh tauh triggers
    cms.PSet (
        HLT = cms.string("HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1LooseChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1LooseChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1TightChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1TightChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1LooseChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1LooseChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleLooseChargedIsoPFTau45_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau45TrackPt1LooseChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau45TrackPt1LooseChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau45_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau45TrackPt1MediumChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau45TrackPt1MediumChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau45_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau45TrackPt1TightChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau45TrackPt1TightChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleLooseChargedIsoPFTau45_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau45TrackPt1LooseChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau45TrackPt1LooseChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau45_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau45TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau45TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau45_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau45TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau45TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
### === VBF tauh tauh triggers (since Run2017D)
    cms.PSet (
        HLT = cms.string("HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau20TrackPt1LooseChargedIsolationReg"),
        path2 = cms.vstring ("hltDoublePFTau20TrackPt1LooseChargedIsolationReg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    )

#now I create the trigger list for HLTconfig
for i in range(len(HLTLIST)):
    tmpl =  str(HLTLIST[i].HLT).replace('cms.string(\'','') ##CMSSW Vaffanculo
    tmpl = tmpl.replace('\')','') ##CMSSW Vaffanculo x 2
    TRIGGERLIST.append(tmpl)
#print TRIGGERLIST
