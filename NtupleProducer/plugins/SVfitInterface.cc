/* \class SVfitInterface
**
** This class provides an interface to the SVfit standalone algorithm
** for the computation of the SVfit mass of the lepton pair candidates.
** 
** The decay mode (e, mu, tauh) of each lepton is the pair is asserted
** from the pdgId associated and is used in the algorithm.
**
** input type is reco::CompositeCandidate for each lepton
** that is coverted to TLorentzVector to be passed to the algorithm
**
** output type is pat::CompositeCandidate, i.e. the original pairs
** plus some userfloats containing the SVfit mass and MET px, px, pt, phi. 
**  
** \date:    18 November 2014
** \author:  L. Cadamuro, O. Davignon (LLR)
*/

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CutSet.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/LeptonIsoHelper.h>
#include <DataFormats/Candidate/interface/ShallowCloneCandidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/METReco/interface/CommonMETData.h>
#include <TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
#include <TLorentzVector.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

// ------------------------------------------------------------------

class SVfitInterface : public edm::EDProducer {
 public:
  /// Constructor
  explicit SVfitInterface(const edm::ParameterSet&);
    
  /// Destructor
  ~SVfitInterface();  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  svFitStandalone::kDecayType GetDecayTypeFlag (int pdgId);
  bool Switch (svFitStandalone::kDecayType type1, double pt1, svFitStandalone::kDecayType type2, double pt2);
  double GetMass (svFitStandalone::kDecayType type, double candMass);

  edm::EDGetTokenT<View<reco::CompositeCandidate> > theCandidateTag;
  // std::vector <edm::EDGetTokenT<View<pat::MET> > > vtheMETTag; // polymorphism of view --> good for reco::PFMET and pat::MET! 
  edm::EDGetTokenT<View<pat::MET> > theMETTag;
  edm::EDGetTokenT<double> theSigTag;
  edm::EDGetTokenT<math::Error<2>::type> theCovTag;
  bool _usePairMET;
  bool _computeForUpDownTES;
  TFile* inputFile_visPtResolution_;
  
};

// ------------------------------------------------------------------


SVfitInterface::SVfitInterface(const edm::ParameterSet& iConfig):
theCandidateTag(consumes<View<reco::CompositeCandidate> >(iConfig.getParameter<InputTag>("srcPairs"))),
//vtheMETTag(consumes<View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
theMETTag(consumes<View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
theSigTag(consumes<double>(iConfig.getParameter<edm::InputTag>("srcSig"))),
theCovTag(consumes<math::Error<2>::type>(iConfig.getParameter<edm::InputTag>("srcCov")))
{
  //theCandidateTag = iConfig.getParameter<InputTag>("srcPairs");
  _usePairMET = iConfig.getParameter<bool>("usePairMET");
  _computeForUpDownTES = iConfig.getParameter<bool>("computeForUpDownTES");
  
  // const std::vector<edm::InputTag>& inMET = iConfig.getParameter<std::vector<edm::InputTag> >("srcMET");
  // for (std::vector<edm::InputTag>::const_iterator it = inMET.begin(); it != inMET.end(); ++it)
  // {      
  //   // vtheMETTag.emplace_back(consumes<edm::View<reco::MET> >(*it) );
  //   vtheMETTag.emplace_back(consumes<edm::View<pat::MET> >(*it) );
  // }

  edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
  TH1::AddDirectory(false);  
  inputFile_visPtResolution_ = new TFile(inputFileName_visPtResolution.fullPath().data());

  produces<pat::CompositeCandidateCollection>();
  
}  

SVfitInterface::~SVfitInterface()
{
    delete inputFile_visPtResolution_;
}

void SVfitInterface::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){  }


svFitStandalone::kDecayType SVfitInterface::GetDecayTypeFlag (int pdgId)
{
    if (abs(pdgId) == 11) return svFitStandalone::kTauToElecDecay;
    if (abs(pdgId) == 13) return svFitStandalone::kTauToMuDecay;
    if (abs(pdgId) == 15) return svFitStandalone::kTauToHadDecay;
    
    edm::LogWarning("WrongDecayModePdgID")
       << "(SVfitInterface): unable to identify decay type from pdgId"
       << "     ---> Decay will be treated as an hadronic decay";
    return svFitStandalone::kTauToHadDecay;
}

// decide if leptons 1 and 2 must be switched to respect SVfit conventions
bool SVfitInterface::Switch (svFitStandalone::kDecayType type1, double pt1, svFitStandalone::kDecayType type2, double pt2)
{
    // e e, mu mu, tau tau
    if (type1 == type2) {return (pt1 < pt2);}
    
    // e tau, mu tau
    if ( (type1 == svFitStandalone::kTauToElecDecay || type1 == svFitStandalone::kTauToMuDecay) &&
         type2 == svFitStandalone::kTauToHadDecay ) {return false;}
    if ( (type2 == svFitStandalone::kTauToElecDecay || type2 == svFitStandalone::kTauToMuDecay) &&
         type1 == svFitStandalone::kTauToHadDecay ) {return true;}

    // e mu
    if (type1 == svFitStandalone::kTauToElecDecay && type2 == svFitStandalone::kTauToMuDecay) {return false;}
    if (type2 == svFitStandalone::kTauToElecDecay && type1 == svFitStandalone::kTauToMuDecay) {return true;}
    
    cout << "SVfit Standalone: ordering not done (should never happen)" << endl;
    return false;
}

// set mass (pdg ele/mu or leave cand mass for tauh
double SVfitInterface::GetMass (svFitStandalone::kDecayType type, double candMass)
{
    if (type == svFitStandalone::kTauToElecDecay) return 0.51100e-3;
    if (type == svFitStandalone::kTauToMuDecay) return 0.10566;

    return candMass; // for tauh and all exceptions return cand mass
}



#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(SVfitInterface);
