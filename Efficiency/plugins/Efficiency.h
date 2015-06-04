#ifndef EFFICIENCY_H
#define EFFICIENCY_H

// -*- C++ -*-
//
// Package:    HLTStudies/Efficiency
// Class:      Efficiency
// 
/**\class Efficiency Efficiency.cc HLTStudies/Efficiency/plugins/Efficiency.cc

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
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// others
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV ;
using namespace std;

//
// class declaration
//

class Efficiency : public edm::EDAnalyzer {
 public:
  explicit Efficiency(const edm::ParameterSet&);
  ~Efficiency();
  
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

  vector<string> _namePaths;
  string         _hltProcessName;
  edm::InputTag  _trigResultsLabel;
  edm::InputTag  _genParticleLabel;
  edm::InputTag  _pfjetCollection;
  edm::InputTag  _pfmetCollection;
  edm::InputTag  _muCollection;
  //edm::InputTag  _vertexCollection;
  //GlobalPoint    _vertexPosition;
  bool _usePtMHT;
  double _minPtJetHt, _maxEtaJetHt, _minPtJetMht, _maxEtaJetMht;

  // 4-vectors
  LV _METP4, _MHTP4, _METNoMuP4, _MHTNoMuP4;

  // Tree and its branches
  TTree* _tree;

  // Global quantities
  int _nEvent, _nRun, _nLumi, _nJet;

  // Trigger info
  TString _trig_pass;
  int _trig_n;

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
