#include "Efficiency.h"

//
// constructors and destructor
//
Efficiency::Efficiency(const edm::ParameterSet& pset) :

  _namePaths(pset.getParameter< vector<string> >("namePaths")),
  _hltProcessName(pset.getParameter<string>("hltProcessName")),
  _trigResultsLabel("TriggerResults", "", _hltProcessName),
  _genParticleLabel(pset.getParameter<string>("genParticleLabel")),
  _pfjetCollection( pset.getParameter<edm::InputTag>("pfjetCollection")),
  _pfmetCollection( pset.getParameter<edm::InputTag>("pfmetCollection")),
  _muCollection(    pset.getParameter<edm::InputTag>("muCollection")),
  //_vertexCollection(pset.getParameter<edm::InputTag>("vertexCollection")),
  _usePtMHT( pset.getParameter<bool>("usePtMHT") ),
  _minPtJetHt( pset.getParameter<double>("minPtJetHt") ), 
  _maxEtaJetHt( pset.getParameter<double>("maxEtaJetHt") ), 
  _minPtJetMht( pset.getParameter<double>("minPtJetMht") ), 
  _maxEtaJetMht( pset.getParameter<double>("maxEtaJetMht") )
{

  //now do what ever initialization is needed
  Init();
  edm::Service<TFileService> fs ;
  _tree = fs->make <TTree>("HLTStudies","Efficiency");

  // Event
  _tree->Branch("nEvent",&_nEvent,"nEvent/I");
  _tree->Branch("nRun",&_nRun,"nRun/I");
  _tree->Branch("nLumi",&_nLumi,"nLumi/I");
  //
  // Trigger
  _tree->Branch("trig_pass",&_trig_pass);
  _tree->Branch("trig_n",&_trig_n,"trig_n/I");  
  //
  // Vertices
  _tree->Branch("vtx_N",&_vtx_N,"vtx_N/I");
  _tree->Branch("vtx_normalizedChi2",&_vtx_normalizedChi2,"vtx_normalizedChi2[3]/D");
  _tree->Branch("vtx_ndof",&_vtx_ndof,"vtx_ndof[3]/D");
  _tree->Branch("vtx_nTracks",&_vtx_nTracks,"vtx_nTracks[3]/D");
  _tree->Branch("vtx_d0",&_vtx_d0,"vtx_d0[3]/D");
  _tree->Branch("vtx_x",&_vtx_x,"vtx_x[3]/D");
  _tree->Branch("vtx_y",&_vtx_y,"vtx_y[3]/D");
  _tree->Branch("vtx_z",&_vtx_z,"vtx_z[3]/D");
  //
  // MET
  _tree->Branch("met", &_met,"met/D");  
  _tree->Branch("mht", &_mht,"mht/D");  
  _tree->Branch("metnomu", &_metnomu,"metnomu/D");  
  _tree->Branch("mhtnomu", &_mhtnomu,"mhtnomu/D");  
  _tree->Branch("met_eta", &_met_eta,"met_eta/D");  
  _tree->Branch("mht_eta", &_mht_eta,"mht_eta/D");  
  _tree->Branch("metnomu_eta", &_metnomu_eta,"metnomu_eta/D");  
  _tree->Branch("mhtnomu_eta", &_mhtnomu_eta,"mhtnomu_eta/D");  
  _tree->Branch("met_phi", &_met_phi,"met_phi/D");  
  _tree->Branch("mht_phi", &_mht_phi,"mht_phi/D");  
  _tree->Branch("metnomu_phi", &_metnomu_phi,"metnomu_phi/D");  
  _tree->Branch("mhtnomu_phi", &_mhtnomu_phi,"mhtnomu_phi/D");  
  _tree->Branch("met_dphi", &_met_dphi,"met_dphi/D");  
  _tree->Branch("mht_dphi", &_mht_dphi,"mht_dphi/D");  
  _tree->Branch("metnomu_dphi", &_metnomu_dphi,"metnomu_dphi/D");  
  _tree->Branch("mhtnomu_dphi", &_mhtnomu_dphi,"mhtnomu_dphi/D");  
  //
  // Jets
  _tree->Branch("nJet",&_nJet,"nJet/I");
  //
  _tree->Branch("jet_eta",&_jet_eta,"jet_eta[nJet]/D");
  _tree->Branch("jet_phi",&_jet_phi,"jet_phi[nJet]/D");
  _tree->Branch("jet_pt",&_jet_pt,"jet_pt[nJet]/D");
  _tree->Branch("jet_e",&_jet_e,"jet_e[nJet]/D");
  _tree->Branch("jet_m",&_jet_m,"jet_m[nJet]/D");
  //
  _tree->Branch("jet_mult_ch",&_jet_mult_ch,"jet_mult_ch[nJet]/I");
  _tree->Branch("jet_mult_mu",&_jet_mult_mu,"jet_mult_mu[nJet]/I");
  _tree->Branch("jet_mult_ne",&_jet_mult_ne,"jet_mult_ne[nJet]/I");
  //
  _tree->Branch("jet_efrac_ne_Had", &_jet_efrac_ne_Had, "jet_efrac_ne_Had[nJet]/D");
  _tree->Branch("jet_efrac_ne_EM",  &_jet_efrac_ne_EM,  "jet_efrac_ne_EM[nJet]/D" );
  _tree->Branch("jet_efrac_ch_Had", &_jet_efrac_ch_Had, "jet_efrac_ch_Had[nJet]/D");
  _tree->Branch("jet_efrac_ch_EM",  &_jet_efrac_ch_EM,  "jet_efrac_ch_EM[nJet]/D" );
  _tree->Branch("jet_efrac_ch_Mu",  &_jet_efrac_ch_Mu,  "jet_efrac_ch_Mu[nJet]/D" );
  //

}


Efficiency::~Efficiency()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Efficiency::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


  //edm::Handle<reco::GenParticleCollection> genPart;

  // Get collections
  edm::Handle<edm::TriggerResults> H_trig;
  iEvent.getByLabel(_trigResultsLabel, H_trig);

  //edm::Handle<reco::VertexCollection> H_vert;
  //iEvent.getByLabel(_vertexCollection, H_vert);

  edm::Handle<reco::PFJetCollection> H_pfjets;
  iEvent.getByLabel(_pfjetCollection , H_pfjets);

  edm::Handle<reco::PFMETCollection> H_pfmet;
  iEvent.getByLabel(_pfmetCollection , H_pfmet);

  edm::Handle<reco::MuonCollection> H_mu;
  iEvent.getByLabel(_muCollection , H_mu);

  // Check validity
  if(!H_trig.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _trigResultsLabel << " ... skip entry !" << endl;
    return;
  }

  if(!H_pfjets.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _pfjetCollection << " ... skip entry !" << endl;
    return;
  }

  if(!H_pfmet.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _pfmetCollection << " ... skip entry !" << endl;
    return;
  }

  if(!H_mu.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _muCollection << " ... skip entry !" << endl;
    return;
  }

  // GLOBAL EVENT INFORMATIONS //
  _nRun   = iEvent.id().run();
  _nLumi  = iEvent.luminosityBlock();
  _nEvent = iEvent.id().event();

  // TRIGGER RESULTS //
  _trig_pass = "";
  _trig_n    = 0;
  TString path="";
  // loop over H_trig
  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*H_trig);
  for (int iHLT = 0 ; iHLT<static_cast<int>(H_trig->size()); ++iHLT) {	
    if (H_trig->accept (iHLT)) {
      path = TString(triggerNames.triggerName(iHLT));
      _trig_pass += "_%_"+path ;
      _trig_n++ ;
    }
  }

  /////////////////////

  // VERTICES //
  /*
  int vtx_counter=0;
  _vtx_N = H_vert->size();
  _vtx_N_stored = nV;
	
  // select the primary vertex as the one with higest sum of (pt)^2 of tracks                                                                               
  PrimaryVertexSorter PVSorter;
  std::vector<reco::Vertex> sortedVertices = PVSorter.sortedList( *(H_vert.product()) );

  if(_vtx_N > 0) {
    GlobalPoint local_vertexPosition(sortedVertices.front().position().x(),
				     sortedVertices.front().position().y(),
				     sortedVertices.front().position().z());
    _vertexPosition = local_vertexPosition;
  }
  else {

  //GlobalPoint local_vertexPosition(bs.position().x(),
  //				     bs.position().y(),
  //				     bs.position().z());

    GlobalPoint local_vertexPosition(0.,0.,0.);
    _vertexPosition = local_vertexPosition;
  }

  for( std::vector<reco::Vertex>::const_iterator PV = sortedVertices.begin(); PV != sortedVertices.end(); ++PV){
    if(vtx_counter > int(nV)) break;
		
    _vtx_normalizedChi2[vtx_counter] = PV->normalizedChi2();
    _vtx_ndof[vtx_counter] = PV->ndof();
    _vtx_nTracks[vtx_counter] = PV->tracksSize();
    _vtx_d0[vtx_counter] = PV->position().Rho();
    _vtx_x[vtx_counter] = PV->x();
    _vtx_y[vtx_counter] = PV->y();
    _vtx_z[vtx_counter] = PV->z();
		
    vtx_counter++;
  } // for loop on primary vertices
*/

  // STORE JET INFORMATION //
  // Loop over PFJets where theJet is a pointer to a PFJet
  // loop only over 3 highest-pt jets
  //
  UInt_t iJ=0;
  //
  for (reco::PFJetCollection::const_iterator theJet = H_pfjets->begin(); theJet != H_pfjets->end(); ++theJet){

    // Kinematics
    _jet_pt[iJ]  = theJet->pt();
    _jet_eta[iJ] = theJet->eta();
    _jet_phi[iJ] = theJet->phi();
    _jet_e[iJ]   = theJet->energy();
    _jet_m[iJ]   = theJet->mass();

    // Energy fractions
    _jet_efrac_ne_Had[iJ] = theJet->neutralHadronEnergyFraction();
    _jet_efrac_ne_EM[ iJ] = theJet->neutralEmEnergyFraction();
    _jet_efrac_ch_Had[iJ] = theJet->chargedHadronEnergyFraction();
    _jet_efrac_ch_EM[ iJ] = theJet->chargedEmEnergyFraction();
    _jet_efrac_ch_Mu[ iJ] = theJet->chargedMuEnergyFraction();

    // Multiplicities
    _jet_mult_ch[iJ] = theJet->chargedMultiplicity();
    _jet_mult_mu[iJ] = theJet->muonMultiplicity();
    _jet_mult_ne[iJ] = theJet->neutralMultiplicity();

    iJ++ ;
    if(iJ>=nJ) break;
  }

  _nJet = nJ;

  // MET INFORMATION //
  const reco::PFMETCollection *C_pfmet = H_pfmet.product();
  _METP4 = (*C_pfmet)[0].p4();

  // MHT
  int nj_ht = 0, nj_mht = 0;
  double ht = 0., mhx = 0., mhy = 0.;
  double pt, eta, phi, px, py;
  //
  for (reco::PFJetCollection::const_iterator theJet = H_pfjets->begin(); theJet != H_pfjets->end(); ++theJet){
    //
    pt = _usePtMHT ? theJet->pt() : theJet->et();
    eta = theJet->eta();
    phi = theJet->phi();
    px = _usePtMHT ? theJet->px() : theJet->et() * cos(phi);
    py = _usePtMHT ? theJet->py() : theJet->et() * sin(phi);
    //
    if (pt > _minPtJetHt && std::abs(eta) < _maxEtaJetHt) {
      ht += pt;
      ++nj_ht;
    }
    //
    if (pt > _minPtJetMht && std::abs(eta) < _maxEtaJetMht) {
      mhx -= px;
      mhy -= py;
      ++nj_mht;
    }
  }

  // NoMu quantities
  double mex_nomu=_METP4.Px();
  double mey_nomu=_METP4.Py();
  double mhx_nomu=mhx;
  double mhy_nomu=mhy;
  //
  for (reco::MuonCollection::const_iterator theMu = H_mu->begin(); theMu != H_mu->end(); ++theMu) {
    mex_nomu += theMu->px();
    mey_nomu += theMu->py();
    mhx_nomu += theMu->px();
    mhy_nomu += theMu->py();
  }

  //if (nj_ht < minNJetHt_ ) { ht = 0; }
  //if (nj_mht < minNJetMht_) { mhx = 0; mhy = 0; }

  // 4-vectors
  _METNoMuP4 = LV(mex_nomu, mey_nomu, 0, sqrt(mex_nomu*mex_nomu + mey_nomu*mey_nomu));
  _MHTNoMuP4 = LV(mhx_nomu, mhy_nomu, 0, sqrt(mhx_nomu*mhx_nomu + mhy_nomu*mhy_nomu));
  _MHTP4     = LV(mhx, mhy, 0, sqrt(mhx*mhx + mhy*mhy));

  // Flat values
  _met     = _METP4.Et();
  _met_eta = _METP4.Eta();
  _met_phi = _METP4.Phi();
  _met_dphi= computeDeltaPhi( _met_phi , _jet_phi[0] );
  //
  _metnomu     = _METNoMuP4.Et();
  _metnomu_eta = _METNoMuP4.Eta();
  _metnomu_phi = _METNoMuP4.Phi();
  _metnomu_dphi= computeDeltaPhi( _metnomu_phi , _jet_phi[0] );
  //
  _mht     = _MHTP4.Et();
  _mht_eta = _MHTP4.Eta();
  _mht_phi = _MHTP4.Phi();
  _mht_dphi= computeDeltaPhi( _mht_phi , _jet_phi[0] );
  //
  _mhtnomu     = _MHTNoMuP4.Et();
  _mhtnomu_eta = _MHTNoMuP4.Eta();
  _mhtnomu_phi = _MHTNoMuP4.Phi();
  _mhtnomu_dphi= computeDeltaPhi( _mhtnomu_phi , _jet_phi[0] );

  
  // FILL TREE //
  _tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
Efficiency::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Efficiency::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
Efficiency::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
Efficiency::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
Efficiency::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
Efficiency::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Efficiency::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
Efficiency::Init()
{

  _nEvent = _nRun = _nLumi = 0;
  _trig_pass = "";
  _trig_n = 0;

  // MET
  _met = _mht = _metnomu = _mhtnomu = 
    _met_eta = _mht_eta = _metnomu_eta = _mhtnomu_eta = 
    _met_phi = _mht_phi = _metnomu_phi = _mhtnomu_phi = 
    _met_dphi = _mht_dphi = _metnomu_dphi = _mhtnomu_dphi = 0;
  
  // Vertices
  _vtx_N = 0; 
  for(UInt_t iv=0;iv<nV;iv++) {
    _vtx_normalizedChi2[iv] = 0.;
    _vtx_ndof[iv] = 0.;
    _vtx_nTracks[iv] = 0.;
    _vtx_d0[iv] = 0.;
    _vtx_x[iv] = 0.;
    _vtx_y[iv] = 0.;
    _vtx_z[iv] = 0.;
  }

  // Jets
  for(UInt_t i=0 ; i<nJ ; i++) {
    _jet_eta[i] = 0;
    _jet_phi[i] = 0;
    _jet_pt[i]  = 0;
    _jet_e[i]   = 0;
    _jet_m[i]   = 0;
    _jet_mult_ch[i] = 0;
    _jet_mult_mu[i] = 0;
    _jet_mult_ne[i] = 0;
    _jet_efrac_ne_Had[i] = 0;
    _jet_efrac_ne_EM[i] = 0;
    _jet_efrac_ch_Had[i] = 0; 
    _jet_efrac_ch_EM[i] = 0; 
    _jet_efrac_ch_Mu[i] = 0;
  }

}

double Efficiency::computeDeltaPhi(double phi1, double phi2)
{
  // Return value in [0;pi] 
  /*
  float dphi0 = TMath::Abs( phi1 - phi2 );
  if(dphi0 > TMath::Pi()) return  TMath::TwoPi() - dphi0;
  else                    return dphi0;
  */
  return TMath::Abs( TVector2::Phi_mpi_pi(phi1 - phi2) );
}

//define this as a plug-in
DEFINE_FWK_MODULE(Efficiency);
