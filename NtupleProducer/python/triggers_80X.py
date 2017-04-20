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
        HLT = cms.string("HLT_IsoMu22_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoTkMu22_v"),
        path1 = cms.vstring ("hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoTkMu22_eta2p1_v"),
        path1 = cms.vstring ("hltL3fL1sMu20erL1f0Tkf22QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu22_eta2p1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu20erL1f0L2f10QL3f22QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),   
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
     cms.PSet (
        HLT = cms.string("HLT_IsoTkMu24_v"),
        path1 = cms.vstring ("hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ), 
    ### ok
### === mu tauh triggers
   cms.PSet (
        HLT = cms.string("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTauJet20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu19LooseIsoPFTau20"),
        path2 = cms.vstring ("hltOverlapFilterIsoMu19LooseIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),   
    cms.PSet (
        HLT = cms.string("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09", "hltOverlapFilterSingleIsoMu19LooseIsoPFTau20"),
        path2 = cms.vstring ("hltOverlapFilterSingleIsoMu19LooseIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),
cms.PSet (
        HLT = cms.string("HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu20erIorSingleMu22erL1f0L2f10QL3f21QL3trkIsoFiltered0p09", "hltOverlapFilterSingleIsoMu21LooseIsoPFTau20"),
        path2 = cms.vstring ("hltOverlapFilterSingleIsoMu21LooseIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),


cms.PSet (
        HLT = cms.string("HLT_VLooseIsoPFTau120_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltPFTau120TrackPt50LooseAbsOrRelVLooseIso"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
  cms.PSet (
        HLT = cms.string("HLT_VLooseIsoPFTau140_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltPFTau140TrackPt50LooseAbsOrRelVLooseIso"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),  
### === tauh tauh triggers
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1MediumIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1MediumIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
 cms.PSet (
        HLT = cms.string("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1MediumCombinedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1MediumCombinedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1MediumCombinedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1MediumCombinedIsolationDz02Reg"),
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
