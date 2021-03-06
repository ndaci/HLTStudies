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
int verbose=1;
const UInt_t nJ=3;
const UInt_t nV=3;

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
  double _vtx_x[nV], _vtx_y[nV], _vtx_z[nV];
  double _vtx_normalizedChi2[nV], _vtx_ndof[nV], _vtx_nTracks[nV], _vtx_d0[nV];

  // MET
  double _met,_mht,_metnomu,_mhtnomu,
    _met_eta,_mht_eta,_metnomu_eta,_mhtnomu_eta,
    _met_phi,_mht_phi,_metnomu_phi,_mhtnomu_phi,
    _met_dphi,_mht_dphi,_metnomu_dphi,_mhtnomu_dphi;

  // Jets
  int _jet_mult_ch[nJ], _jet_mult_mu[nJ], _jet_mult_ne[nJ]; // multiplicities
  double _jet_eta[nJ], _jet_phi[nJ], _jet_pt[nJ], _jet_e[nJ], _jet_m[nJ];
  double _jet_efrac_ne_Had[nJ], _jet_efrac_ne_EM[nJ]; // neutral energy fractions
  double _jet_efrac_ch_Had[nJ], _jet_efrac_ch_EM[nJ], _jet_efrac_ch_Mu[nJ]; // charged energy fractions
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

