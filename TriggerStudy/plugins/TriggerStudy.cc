// -*- C++ -*-
//
// Package:    HLTStudies/TriggerStudy
// Class:      TriggerStudy
// 
/**\class TriggerStudy TriggerStudy.cc HLTStudies/TriggerStudy/plugins/TriggerStudy.cc

 Description: Trigger performance studies

*/
//
// Original Author:  Nadir Daci
//         Created:  Tue, 14 May 2019 16:05:48 GMT
//
//


// C++
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>

// Framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
//#include "CLHEP/Matrix/Vector.h"

// HLT info
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

// Gen Info
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

// DataFormats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

// Level-1
#include "DataFormats/L1Trigger/interface/EGamma.h" 
#include "DataFormats/L1Trigger/interface/Tau.h"    
#include "DataFormats/L1Trigger/interface/Jet.h"    
#include "DataFormats/L1Trigger/interface/Muon.h"   
#include "DataFormats/L1Trigger/interface/EtSum.h"  
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmAlgorithm.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

// ROOT
#include "TH1F.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TPRegexp.h"

// Namespaces
using namespace edm;
using namespace reco;
using namespace std;
using namespace pat;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TriggerStudy : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

public:
  explicit TriggerStudy(const edm::ParameterSet&);
  ~TriggerStudy();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------


  // InputTags //

  // Trigger and MET filters
  const edm::InputTag triggerResultsTag;
  const edm::InputTag filterResultsTag;
  const edm::InputTag prescalesTag;
  const edm::InputTag triggerObjectsTag; 
  const edm::InputTag IT_L1_EG;  
  const edm::InputTag IT_L1_Tau; 
  const edm::InputTag IT_L1_Jet; 
  const edm::InputTag IT_L1_Mu;  
  const edm::InputTag IT_L1_Sums;
  const edm::InputTag IT_L1_Algos; 

  // Objects
  const edm::InputTag muonsTag;
  const edm::InputTag electronsTag;
  const edm::InputTag photonsTag;
  const edm::InputTag tausTag;
  const edm::InputTag jetsTag;
  const edm::InputTag t1metTag;

  
  // Tokens //

  // Trigger and MET filters
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken; 
  edm::EDGetTokenT<edm::TriggerResults> filterResultsToken;
  edm::EDGetTokenT<l1t::EGammaBxCollection> T_L1EG;
  edm::EDGetTokenT<l1t::TauBxCollection>    T_L1Tau;
  edm::EDGetTokenT<l1t::JetBxCollection>    T_L1Jet;
  edm::EDGetTokenT<l1t::MuonBxCollection>   T_L1Mu;
  edm::EDGetTokenT<l1t::EtSumBxCollection>  T_L1Sums;
  edm::EDGetTokenT<GlobalAlgBlkBxCollection> T_L1Algos;

  // Objects
  edm::EDGetTokenT<pat::MuonCollection> muonsToken;
  edm::EDGetTokenT<pat::ElectronCollection> electronsToken;
  edm::EDGetTokenT<pat::PhotonCollection> photonsToken;
  edm::EDGetTokenT<pat::TauCollection> tausToken;
  edm::EDGetTokenT<std::vector<pat::Jet> > jetsToken;
  edm::EDGetTokenT<edm::View<pat::MET> > t1metToken;


  // Variables //
  uint32_t event, run, lumi; 

  // Utilities //

  // HLT PS
  std::unique_ptr<HLTPrescaleProvider> hltPrescaleProvider_;

  // Sorting objects
  template<typename T> 
  class PatPtSorter{
  public:
    bool operator ()(const T & i, const T & j) const {
      return (i->pt() > j->pt());
    }
  };
  //
  PatPtSorter<pat::JetRef>      jetSorter;
  PatPtSorter<pat::MuonRef>     muonSorter;
  PatPtSorter<pat::ElectronRef> electronSorter;
  PatPtSorter<pat::PhotonRef>   photonSorter;
  PatPtSorter<pat::TauRef>      tauSorter;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TriggerStudy::TriggerStudy(const edm::ParameterSet& iConfig):
  
  triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerResults")),
  filterResultsTag(iConfig.getParameter<edm::InputTag>("filterResults")),
  prescalesTag(iConfig.getParameter<edm::InputTag>("prescales")),
  triggerObjectsTag(iConfig.existsAs<edm::InputTag>("triggerObjects") ? 
		    iConfig.getParameter<edm::InputTag>("triggerObjects") : edm::InputTag("")), 

  IT_L1_EG(  iConfig.existsAs<edm::InputTag>(    "triggerL1EG")   ? 
	     iConfig.getParameter<edm::InputTag>("triggerL1EG")   : edm::InputTag("")), 
  IT_L1_Tau( iConfig.existsAs<edm::InputTag>(    "triggerL1Tau")  ? 
	     iConfig.getParameter<edm::InputTag>("triggerL1Tau")  : edm::InputTag("")), 
  IT_L1_Jet( iConfig.existsAs<edm::InputTag>(    "triggerL1Jet")  ? 
	     iConfig.getParameter<edm::InputTag>("triggerL1Jet")  : edm::InputTag("")), 
  IT_L1_Mu(  iConfig.existsAs<edm::InputTag>(    "triggerL1Mu")   ? 
	     iConfig.getParameter<edm::InputTag>("triggerL1Mu")   : edm::InputTag("")), 
  IT_L1_Sums(iConfig.existsAs<edm::InputTag>(    "triggerL1Sums") ? 
	     iConfig.getParameter<edm::InputTag>("triggerL1Sums") : edm::InputTag("")), 
  IT_L1_Algos(iConfig.existsAs<edm::InputTag>("triggerL1algos") ?
	      iConfig.getParameter<edm::InputTag>("triggerL1algos") : edm::InputTag("gtStage2Digis")),

  muonsTag(iConfig.getParameter<edm::InputTag>("muons")),
  electronsTag(iConfig.getParameter<edm::InputTag>("electrons")),
  photonsTag(iConfig.getParameter<edm::InputTag>("photons")),
  tausTag(iConfig.getParameter<edm::InputTag>("taus")),
  jetsTag(iConfig.getParameter<edm::InputTag>("jets")),
  t1metTag(iConfig.getParameter<edm::InputTag>("t1met"))

{

  //usesResource("TFileService");

  // Tokens //

  // Trigger and MET filters
  triggerResultsToken   = consumes<edm::TriggerResults> (triggerResultsTag);
  triggerPrescalesToken = consumes<pat::PackedTriggerPrescales>(prescalesTag);
  triggerObjectsToken = consumes<pat::TriggerObjectStandAloneCollection>(triggerObjectsTag); 
  T_L1EG              = consumes<l1t::EGammaBxCollection>(IT_L1_EG  ); 
  T_L1Tau             = consumes<l1t::TauBxCollection   >(IT_L1_Tau ); 
  T_L1Jet             = consumes<l1t::JetBxCollection   >(IT_L1_Jet ); 
  T_L1Mu              = consumes<l1t::MuonBxCollection  >(IT_L1_Mu  ); 
  T_L1Sums            = consumes<l1t::EtSumBxCollection >(IT_L1_Sums); 
  T_L1Algos           = consumes<GlobalAlgBlkBxCollection>(IT_L1_Algos); 
  filterResultsToken    = consumes<edm::TriggerResults> (filterResultsTag);

  // Instantiate the hltPrescaleProvider
  hltPrescaleProvider_.reset(new HLTPrescaleProvider(iConfig, consumesCollector(), *this)); 

  // Objects
  muonsToken = consumes<pat::MuonCollection> (muonsTag);
  electronsToken = consumes<pat::ElectronCollection> (electronsTag);
  photonsToken = consumes<pat::PhotonCollection> (photonsTag);
  tausToken = consumes<pat::TauCollection> (tausTag);
  jetsToken = consumes<std::vector<pat::Jet> > (jetsTag);
  t1metToken = consumes<edm::View<pat::MET> > (t1metTag);

}


TriggerStudy::~TriggerStudy()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggerStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Handles //

  // Trigger
  Handle<TriggerResults> triggerResultsH;
  iEvent.getByToken(triggerResultsToken, triggerResultsH);

  Handle<pat::PackedTriggerPrescales> triggerPrescalesH;
  iEvent.getByToken(triggerPrescalesToken, triggerPrescalesH);

  Handle<pat::TriggerObjectStandAloneCollection>   triggerObjectsH; 
  iEvent.getByToken(triggerObjectsToken, triggerObjectsH); 

  const edm::TriggerNames &trignames = iEvent.triggerNames(*triggerResultsH);

  // MET filters
  Handle<TriggerResults> filterResultsH;
  iEvent.getByToken(filterResultsToken, filterResultsH);
  
  // L1 Trigger Candidates
  edm::Handle<l1t::EGammaBxCollection> H_L1EG;
  edm::Handle<l1t::TauBxCollection>    H_L1Tau;
  edm::Handle<l1t::JetBxCollection>    H_L1Jet;
  edm::Handle<l1t::MuonBxCollection>   H_L1Mu;
  edm::Handle<l1t::EtSumBxCollection>  H_L1Sums;
  edm::Handle<GlobalAlgBlkBxCollection> H_L1Algos;
  //
  iEvent.getByToken(T_L1EG  , H_L1EG);
  iEvent.getByToken(T_L1Tau , H_L1Tau);
  iEvent.getByToken(T_L1Jet , H_L1Jet);
  iEvent.getByToken(T_L1Mu  , H_L1Mu);
  iEvent.getByToken(T_L1Sums, H_L1Sums);
  iEvent.getByToken(T_L1Algos, H_L1Algos);

  // Objects
  Handle<pat::MuonCollection> muonsH;
  iEvent.getByToken(muonsToken, muonsH);
  pat::MuonCollection muons = *muonsH;
  
  Handle<pat::ElectronCollection> electronsH;
  iEvent.getByToken(electronsToken, electronsH);
  pat::ElectronCollection electrons = *electronsH;

  Handle<pat::PhotonCollection> photonsH;
  iEvent.getByToken(photonsToken, photonsH);
  pat::PhotonCollection photons = *photonsH;

  Handle<pat::TauCollection > tausH;
  iEvent.getByToken(tausToken, tausH);
  pat::TauCollection taus = *tausH;

  Handle<vector<pat::Jet> > jetsH;
  iEvent.getByToken(jetsToken, jetsH);

  Handle<View<pat::MET> > t1metH;
  iEvent.getByToken(t1metToken, t1metH);


  // Look into the event information //

  // Event
  event = iEvent.id().event();
  run   = iEvent.id().run();
  lumi = iEvent.luminosityBlock();

  // Trigger objects
  string trgColl="";
  for(pat::TriggerObjectStandAlone obj : *triggerObjectsH) { 
    obj.unpackPathNames(trignames);
    trgColl = obj.collection();
    cout << "\tTrigger object:  " << trgColl << " | pt = " << obj.pt() 
	 << " | eta = " << obj.eta() << " | phi = " << obj.phi() 
	 << endl;
  }


}


// ------------ method called once each job just before starting event loop  ------------
void 
TriggerStudy::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerStudy::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  //desc.setUnknown();
  //descriptions.addDefault(desc);

  // Trigger and MET filters
  desc.add<edm::InputTag>("triggerResults", edm::InputTag("TriggerResults", "", "HLT"));
  desc.add<edm::InputTag>("filterResults" , edm::InputTag("TriggerResults", "", "RECO"));
  desc.add<edm::InputTag>("prescales"     , edm::InputTag("patTrigger"    , "", "RECO"));
  desc.add<edm::InputTag>("triggerObjects", edm::InputTag("slimmedPatTrigger", "", "RECO")); 
  desc.add<edm::InputTag>("triggerL1EG"         , edm::InputTag("caloStage2Digis", "EGamma", "RECO"));  
  desc.add<edm::InputTag>("triggerL1Tau"        , edm::InputTag("caloStage2Digis", "Tau"   , "RECO")); 
  desc.add<edm::InputTag>("triggerL1Jet"        , edm::InputTag("caloStage2Digis", "Jet"   , "RECO")); 
  desc.add<edm::InputTag>("triggerL1Sums"       , edm::InputTag("caloStage2Digis", "EtSum" , "RECO"));
  desc.add<edm::InputTag>("triggerL1Mu"         , edm::InputTag("gmtStage2Digis" , "Muon"  , "RECO"));  
  desc.add<edm::InputTag>("triggerL1algos"      , edm::InputTag("gtStage2Digis"  , ""      , "RECO")); 

  // Objects
  desc.add<edm::InputTag>("muons"    , edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("electrons", edm::InputTag("slimmedElectrons"));
  desc.add<edm::InputTag>("photons"  , edm::InputTag("slimmedPhotons"));
  desc.add<edm::InputTag>("taus"     , edm::InputTag("slimmedTaus"));
  desc.add<edm::InputTag>("jets"     , edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("t1met"    , edm::InputTag("slimmedMETs"));

  descriptions.add("triggerStudy", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerStudy);
