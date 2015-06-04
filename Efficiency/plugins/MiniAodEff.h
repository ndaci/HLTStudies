#ifndef MINIAODEFF_H
#define MINIAODEFF_H

// -*- C++ -*-
//
// Package:    HLTStudies/MiniAodEff
// Class:      MiniAodEff
// 
/**\class MiniAodEff MiniAodEff.cc HLTStudies/MiniAodEff/plugins/MiniAodEff.cc

 Description: [one line class summary]

*/
//
// Original Author:  ndaci
//         Created:  Wed, 07 Jan 2015 13:48:51 GMT
//
//


// C++ lib
#include <memory>

// ROOT
#include "TTree.h"
#include "TVector2.h"
#include "TClonesArray.h"

// CMSSW standard lib
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/LorentzVector.h"

// CMSSW specific lib
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// others
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV ;
using namespace std;

//
// class declaration
//

class MiniAodEff : public edm::EDAnalyzer {
 public:
  explicit MiniAodEff(const edm::ParameterSet&);
  ~MiniAodEff();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void Init();  
  double computeDeltaPhi(double phi1, double phi2);

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------

  // "global" variables
  int _verbose;
  const UInt_t _nJ=3;
  const UInt_t _nV=3;

  // Input tags
  vector<string> _namePaths;
  string         _hltProcessName;
  edm::InputTag  _trigResultsLabel;

  // InputTags for Tokens
  edm::InputTag _IT_trg_bits, _IT_trg_ps, _IT_trg_obj; 
  edm::InputTag _IT_vert ;
  edm::InputTag _IT_mu   ; 
  edm::InputTag _IT_ele  ; 
  edm::InputTag _IT_tau  ; 
  edm::InputTag _IT_pho  ; 
  edm::InputTag _IT_jet  ; 
  edm::InputTag _IT_gjet ; 
  edm::InputTag _IT_met  ; 

  // Tokens
  edm::EDGetTokenT<edm::TriggerResults> trgBitsToken_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales>            trgPrescalesToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trgObjectsToken_;
  //
  edm::EDGetTokenT<reco::VertexCollection>  vtxToken_;
  edm::EDGetTokenT<pat::MuonCollection>     muonToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<pat::TauCollection>      tauToken_;
  edm::EDGetTokenT<pat::PhotonCollection>   photonToken_;
  edm::EDGetTokenT<pat::JetCollection>      jetToken_;
  edm::EDGetTokenT<reco::GenJetCollection>  genjetToken_;
  edm::EDGetTokenT<pat::METCollection>      metToken_;

  bool _usePtMHT;
  double _minPtJetHt, _maxEtaJetHt, _minPtJetMht, _maxEtaJetMht;

  // To fill the tree

  // 4-vectors
  LV _METP4, _MHTP4, _METNoMuP4, _MHTNoMuP4;

  // Tree and its branches
  TTree* _tree;

  // Global quantities
  int _nEvent, _nRun, _nLumi, _nJet;

  // Trigger info
  TString _trig_pass;
  int _trig_n;
  //std::vector< pat::TriggerObjectStandAlone > *_trig_obj;
  TClonesArray* _trig_obj ;

  // Vertices
  int _vtx_N, _vtx_N_stored;
  double _vtx_x[3], _vtx_y[3], _vtx_z[3];
  double _vtx_normalizedChi2[3], _vtx_ndof[3], _vtx_nTracks[3], _vtx_d0[3];

  // MET
  double _met,_mht,_metnomu,_mhtnomu,
    _met_eta,_mht_eta,_metnomu_eta,_mhtnomu_eta,
    _met_phi,_mht_phi,_metnomu_phi,_mhtnomu_phi,
    _met_dphi,_mht_dphi,_metnomu_dphi,_mhtnomu_dphi;

  // Jets
  int _jet_mult_ch[3], _jet_mult_mu[3], _jet_mult_ne[3]; // multiplicities
  double _jet_eta[3], _jet_phi[3], _jet_pt[3], _jet_e[3], _jet_m[3];
  double _jet_efrac_ne_Had[3], _jet_efrac_ne_EM[3]; // neutral energy fractions
  double _jet_efrac_ch_Had[3], _jet_efrac_ch_EM[3], _jet_efrac_ch_Mu[3]; // charged energy fractions
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

#endif
