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
#include <sstream>

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
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<Double_t> > LV ;
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
  Double_t computeDeltaPhi(Double_t phi1, Double_t phi2);
  string concatenate(vector<string> vstring);

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------

  // "global" variables
  Int_t _verbose;
  const UInt_t _nJ=3;
  const UInt_t _nV=100;

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

  Bool_t _usePtMHT;
  Double_t _minPtJetHt, _maxEtaJetHt, _minPtJetMht, _maxEtaJetMht;

  // To fill the tree

  // 4-vectors
  LV _METP4, _MHTP4, _METNoMuP4, _MHTNoMuP4;

  // Tree and its branches
  TTree* _tree;

  // Global quantities
  Int_t _nEvent, _nRun, _nLumi;

  // Trigger info
  TString _trig_pass;
  Int_t _trig_n;
  Int_t _trig_obj_n;
  vector< double > _trig_obj_pt, _trig_obj_eta, _trig_obj_phi;
  vector< string > _trig_obj_col, _trig_obj_lab;
  vector< string > _trig_obj_path_FF, _trig_obj_path_FT, 
    _trig_obj_path_TF, _trig_obj_path_TT ;
  vector< vector<int> > _trig_obj_ids;

  // Vertices
  Int_t _vtx_N, _vtx_N_stored;
  Double_t _vtx_x[100], _vtx_y[100], _vtx_z[100];
  Double_t _vtx_normalizedChi2[100], _vtx_ndof[100], _vtx_nTracks[100], _vtx_d0[100];

  // MET
  Double_t _met,_mht,_metnomu,_mhtnomu,
    _met_phi,_mht_phi,_metnomu_phi,_mhtnomu_phi,
    _met_dphi,_mht_dphi,_metnomu_dphi,_mhtnomu_dphi;

  // Jets
  Int_t _jet_N;
  vector<int>    _jet_mult_ch, _jet_mult_mu, _jet_mult_ne; // multiplicities
  vector<double> _jet_eta, _jet_phi, _jet_pt, _jet_e, _jet_m;
  vector<double> _jet_efrac_ne_Had, _jet_efrac_ne_EM; // neutral energy fractions
  vector<double> _jet_efrac_ch_Had, _jet_efrac_ch_EM, _jet_efrac_ch_Mu; // charged energy fractions
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

#endif
