#include "MiniAodEff.h"

//
// constructors and destructor
//
MiniAodEff::MiniAodEff(const edm::ParameterSet& pset)
{

  // Get Input parameters
  _namePaths = pset.getParameter< vector<string> >("namePaths");
  _hltProcessName = pset.getParameter<string>("hltProcessName");
  _trigResultsLabel = edm::InputTag("TriggerResults", "", _hltProcessName);

  // InputTags for tokens
  _IT_trg_bits = pset.getParameter<edm::InputTag>("bits");
  _IT_trg_ps   = pset.getParameter<edm::InputTag>("prescales");
  _IT_trg_obj  = pset.getParameter<edm::InputTag>("objects");
  //
  _IT_vert = pset.getParameter<edm::InputTag>("vertices");
  _IT_mu   = pset.getParameter<edm::InputTag>("muons");
  _IT_ele  = pset.getParameter<edm::InputTag>("electrons");
  _IT_tau  = pset.getParameter<edm::InputTag>("taus");
  _IT_pho  = pset.getParameter<edm::InputTag>("photons");
  _IT_jet  = pset.getParameter<edm::InputTag>("jets");
  _IT_gjet = pset.getParameter<edm::InputTag>("genjets");
  _IT_met  = pset.getParameter<edm::InputTag>("met");

  // Tokens
  trgBitsToken_      = consumes<edm::TriggerResults>(_IT_trg_bits);
  trgPrescalesToken_ = consumes<pat::PackedTriggerPrescales>(_IT_trg_ps);
  trgObjectsToken_   = consumes<pat::TriggerObjectStandAloneCollection>(_IT_trg_obj);
  //
  vtxToken_     = consumes<reco::VertexCollection>(_IT_vert);
  muonToken_    = consumes<pat::MuonCollection>(_IT_mu);
  electronToken_= consumes<pat::ElectronCollection>(_IT_ele);
  tauToken_     = consumes<pat::TauCollection>(_IT_tau);
  photonToken_  = consumes<pat::PhotonCollection>(_IT_pho);
  jetToken_     = consumes<pat::JetCollection>(_IT_jet);
  genjetToken_  = consumes<reco::GenJetCollection>(_IT_gjet);
  metToken_     = consumes<pat::METCollection>(_IT_met);
  //
  _usePtMHT     = pset.getParameter<bool>("usePtMHT") ;
  _minPtJetHt   = pset.getParameter<double>("minPtJetHt") ; 
  _maxEtaJetHt  = pset.getParameter<double>("maxEtaJetHt") ; 
  _minPtJetMht  = pset.getParameter<double>("minPtJetMht") ; 
  _maxEtaJetMht = pset.getParameter<double>("maxEtaJetMht") ;

  //now do what ever initialization is needed
  Init();
  edm::Service<TFileService> fs ;
  _tree = fs->make <TTree>("HLTStudies","MiniAodEff");

  // Event
  _tree->Branch("nEvent",&_nEvent,"nEvent/I");
  _tree->Branch("nRun",&_nRun,"nRun/I");
  _tree->Branch("nLumi",&_nLumi,"nLumi/I");
  //
  // Trigger
  _tree->Branch("trig_pass",&_trig_pass);
  _tree->Branch("trig_n",&_trig_n,"trig_n/I");  

  //_trig_obj = new std::vector< pat::TriggerObjectStandAlone >;
  //_tree->Branch("trig_obj", "std::vector< pat::TriggerObjectStandAlone >", &_trig_obj); // trigger objects
  _trig_obj = new TClonesArray("pat::TriggerObjectStandAlone",10000);
  _tree->Branch("trig_obj", &_trig_obj); // trigger objects
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


MiniAodEff::~MiniAodEff()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  
  _trig_obj->Delete();

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MiniAodEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Get collections //
  // Trigger
  edm::Handle<edm::TriggerResults> H_trig;
  iEvent.getByLabel(_trigResultsLabel, H_trig);
  
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(trgBitsToken_, triggerBits);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(trgObjectsToken_, triggerObjects);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(trgPrescalesToken_, triggerPrescales);

  edm::Handle<reco::VertexCollection> H_vert;
  iEvent.getByToken(vtxToken_, H_vert);
  
  edm::Handle<pat::MuonCollection> H_mu;
  iEvent.getByToken(muonToken_, H_mu);

  edm::Handle<pat::ElectronCollection> H_ele;
  iEvent.getByToken(electronToken_, H_ele);

  edm::Handle<pat::PhotonCollection> H_pho;
  iEvent.getByToken(photonToken_, H_pho);

  edm::Handle<pat::TauCollection> H_tau;
  iEvent.getByToken(tauToken_, H_tau);

  edm::Handle<pat::JetCollection> H_jets;
  iEvent.getByToken(jetToken_, H_jets);

  edm::Handle<reco::GenJetCollection> H_genjets;
  iEvent.getByToken(genjetToken_, H_genjets);

  edm::Handle<pat::METCollection> H_met;
  iEvent.getByToken(metToken_, H_met);

  // Check validity
  if(!H_trig.isValid()) {
    if(_verbose>0) cout << "Missing collection : " << _trigResultsLabel << " ... skip entry !" << endl;
    return;
  }

  if(!H_vert.isValid()) {
    if(_verbose>0) cout << "Missing collection : " << _IT_vert << " ... skip entry !" << endl;
    return;
  }

  if(!H_mu.isValid()) {
    if(_verbose>0) cout << "Missing collection : " << _IT_mu << " ... skip entry !" << endl;
    return;
  }

  if(!H_ele.isValid()) {
    if(_verbose>0) cout << "Missing collection : " << _IT_ele << " ... skip entry !" << endl;
    return;
  }

  if(!H_pho.isValid()) {
    if(_verbose>0) cout << "Missing collection : " << _IT_pho << " ... skip entry !" << endl;
    return;
  }

  if(!H_tau.isValid()) {
    if(_verbose>0) cout << "Missing collection : " << _IT_tau << " ... skip entry !" << endl;
    return;
  }

  if(!H_jets.isValid()) {
    if(_verbose>0) cout << "Missing collection : " << _IT_jet << " ... skip entry !" << endl;
    return;
  }

  if(!H_genjets.isValid()) {
    if(_verbose>0) cout << "Missing collection : " << _IT_gjet << " ... skip entry !" << endl;
    return;
  }

  if(!H_met.isValid()) {
    if(_verbose>0) cout << "Missing collection : " << _IT_met << " ... skip entry !" << endl;
    return;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////

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

  // TRIGGER OBJECTS //
  _trig_obj->Clear();
  UInt_t iTrgObj=0;
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames

    //_trig_obj->push_back(obj);
    new( (*_trig_obj)[iTrgObj] ) pat::TriggerObjectStandAlone(obj);
    iTrgObj++ ;

    obj.unpackPathNames(triggerNames);

    cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << endl;
    // Print trigger object collection and type
    cout << "\t   Collection: " << obj.collection() << endl;
    cout << "\t   Type IDs:   ";
    for (unsigned h = 0; h < obj.filterIds().size(); ++h) cout << " " << obj.filterIds()[h] ;
    cout << endl;
    // Print associated trigger filters
    cout << "\t   Filters:    ";
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h) cout << " " << obj.filterLabels()[h];
    cout << endl;
    vector<string> pathNamesAll  = obj.pathNames(false);
    vector<string> pathNamesLast = obj.pathNames(true);
    // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
    // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
    // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
    cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
    for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
      bool isBoth = obj.hasPathName( pathNamesAll[h], true, true ); 
      bool isL3   = obj.hasPathName( pathNamesAll[h], false, true ); 
      bool isLF   = obj.hasPathName( pathNamesAll[h], true, false ); 
      bool isNone = obj.hasPathName( pathNamesAll[h], false, false ); 
      cout << "   " << pathNamesAll[h];
      if (isBoth) cout << "(L,3)";
      if (isL3 && !isBoth) cout << "(*,3)";
      if (isLF && !isBoth) cout << "(L,*)";
      if (isNone && !isBoth && !isL3 && !isLF) cout << "(*,*)";
    }
    cout << endl;
  }  

  /////////////////////

  // VERTICES //
  /*
  int vtx_counter=0;
  _vtx_N = H_vert->size();
  _vtx_N_stored = _nV;
	
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
    if(vtx_counter > int(_nV)) break;
		
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
  for (pat::JetCollection::const_iterator theJet = H_jets->begin(); theJet != H_jets->end(); ++theJet){

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
    if(iJ>=_nJ) break;
  }

  _nJet = _nJ;

  // MET INFORMATION //
  const pat::METCollection *C_pfmet = H_met.product();
  _METP4 = (*C_pfmet)[0].p4();

  // MHT
  int nj_ht = 0, nj_mht = 0;
  double ht = 0., mhx = 0., mhy = 0.;
  double pt, eta, phi, px, py;
  //
  for (pat::JetCollection::const_iterator theJet = H_jets->begin(); theJet != H_jets->end(); ++theJet){
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
  for (pat::MuonCollection::const_iterator theMu = H_mu->begin(); theMu != H_mu->end(); ++theMu) {
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
MiniAodEff::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAodEff::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
MiniAodEff::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MiniAodEff::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MiniAodEff::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MiniAodEff::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAodEff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
MiniAodEff::Init()
{

  _verbose = 1;
  //_nJ=3;
  //_nV=3;

  _nEvent = _nRun = _nLumi = 0;
  _trig_pass = "";
  _trig_n = 0;
  //_trig_obj = new vector< pat::TriggerObjectStandAlone >;

  // MET
  _met = _mht = _metnomu = _mhtnomu = 
    _met_eta = _mht_eta = _metnomu_eta = _mhtnomu_eta = 
    _met_phi = _mht_phi = _metnomu_phi = _mhtnomu_phi = 
    _met_dphi = _mht_dphi = _metnomu_dphi = _mhtnomu_dphi = 0;
  
  // Vertices
  _vtx_N = 0; 
  for(UInt_t iv=0;iv<_nV;iv++) {
    _vtx_normalizedChi2[iv] = 0.;
    _vtx_ndof[iv] = 0.;
    _vtx_nTracks[iv] = 0.;
    _vtx_d0[iv] = 0.;
    _vtx_x[iv] = 0.;
    _vtx_y[iv] = 0.;
    _vtx_z[iv] = 0.;
  }

  // Jets
  for(UInt_t i=0 ; i<_nJ ; i++) {
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

double MiniAodEff::computeDeltaPhi(double phi1, double phi2)
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
DEFINE_FWK_MODULE(MiniAodEff);
