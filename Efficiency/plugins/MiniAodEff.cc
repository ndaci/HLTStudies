#include "MiniAodEff.h"

//
// constructors and destructor
//
MiniAodEff::MiniAodEff(const edm::ParameterSet& pset)
{

  // Get Input parameters
  _namePaths = pset.getParameter< vector<string> >("namePaths");
  _verbose   = pset.getParameter< int >("verbose");

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
  _IT_met_filt = pset.getParameter<edm::InputTag>("metfilter");

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
  metFiltToken_ = consumes<edm::TriggerResults>(_IT_met_filt);
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
  Int_t buffersize = 32000; // trig_obj vectors

  // Event
  _tree->Branch("nEvent",&_nEvent,"nEvent/I");
  _tree->Branch("nRun",&_nRun,"nRun/I");
  _tree->Branch("nLumi",&_nLumi,"nLumi/I");
  //
  // Trigger
  _tree->Branch("trig_pass",&_trig_pass);
  _tree->Branch("trig_n",&_trig_n,"trig_n/I");  
  //
  _tree->Branch("trig_obj_n",&_trig_obj_n,"trig_obj_n/I");
  _tree->Branch("trig_obj_pt","std::vector<double>",&_trig_obj_pt,buffersize);
  _tree->Branch("trig_obj_eta","std::vector<double>",&_trig_obj_eta,buffersize);
  _tree->Branch("trig_obj_phi","std::vector<double>",&_trig_obj_phi,buffersize);
  _tree->Branch("trig_obj_col","std::vector<std::string>",&_trig_obj_col,buffersize);
  _tree->Branch("trig_obj_ids","std::vector<std::vector<std::int>>",&_trig_obj_ids,buffersize);
  _tree->Branch("trig_obj_lab", "std::vector<std::string>",&_trig_obj_lab, buffersize);
  _tree->Branch("trig_obj_path_FF","std::vector<std::string>",&_trig_obj_path_FF,buffersize);
  _tree->Branch("trig_obj_path_FT","std::vector<std::string>",&_trig_obj_path_FT,buffersize);
  _tree->Branch("trig_obj_path_TF","std::vector<std::string>",&_trig_obj_path_TF,buffersize);
  _tree->Branch("trig_obj_path_TT","std::vector<std::string>",&_trig_obj_path_TT,buffersize);

  //
  // Vertices
  _tree->Branch("vtx_N",&_vtx_N,"vtx_N/I");
  _tree->Branch("vtx_normalizedChi2",&_vtx_normalizedChi2,"vtx_normalizedChi2[vtx_N]/D");
  _tree->Branch("vtx_ndof",&_vtx_ndof,"vtx_ndof[vtx_N]/D");
  _tree->Branch("vtx_nTracks",&_vtx_nTracks,"vtx_nTracks[vtx_N]/D");
  _tree->Branch("vtx_d0",&_vtx_d0,"vtx_d0[vtx_N]/D");
  _tree->Branch("vtx_x",&_vtx_x,"vtx_x[vtx_N]/D");
  _tree->Branch("vtx_y",&_vtx_y,"vtx_y[vtx_N]/D");
  _tree->Branch("vtx_z",&_vtx_z,"vtx_z[vtx_N]/D");
  //
  // MET
  _tree->Branch("met_filters","std::vector<std::string>",&_met_filters,buffersize);
  _tree->Branch("met", &_met,"met/D");  
  _tree->Branch("mht", &_mht,"mht/D");  
  _tree->Branch("metnomu", &_metnomu,"metnomu/D");  
  _tree->Branch("mhtnomu", &_mhtnomu,"mhtnomu/D");  
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
  _tree->Branch("jet_N",&_jet_N,"jet_N/I");
  //
  _tree->Branch("jet_eta","std::vector<double>",&_jet_eta,buffersize);
  _tree->Branch("jet_phi","std::vector<double>",&_jet_phi,buffersize);
  _tree->Branch("jet_pt","std::vector<double>",&_jet_pt,buffersize);
  _tree->Branch("jet_e","std::vector<double>",&_jet_e,buffersize);
  _tree->Branch("jet_m","std::vector<double>",&_jet_m,buffersize);
  //
  _tree->Branch("jet_mult_ch","std::vector<int>",&_jet_mult_ch,buffersize);
  _tree->Branch("jet_mult_mu","std::vector<int>",&_jet_mult_mu,buffersize);
  _tree->Branch("jet_mult_ne","std::vector<int>",&_jet_mult_ne,buffersize);
  //
  _tree->Branch("jet_efrac_ne_Had","std::vector<double>", &_jet_efrac_ne_Had,buffersize);
  _tree->Branch("jet_efrac_ne_EM","std::vector<double>",  &_jet_efrac_ne_EM,buffersize);
  _tree->Branch("jet_efrac_ch_Had","std::vector<double>", &_jet_efrac_ch_Had,buffersize);
  _tree->Branch("jet_efrac_ch_EM","std::vector<double>",  &_jet_efrac_ch_EM,buffersize);
  _tree->Branch("jet_efrac_ch_Mu","std::vector<double>",  &_jet_efrac_ch_Mu,buffersize);
  //
  // Leptons/Photons
  _tree->Branch("nphotons",&_nphotons,"nphotons/I");
  _tree->Branch("nelectrons",&_nelectrons,"nelectrons/I");
  _tree->Branch("nmuons",&_nmuons,"nmuons/I");
  _tree->Branch("ntaus",&_ntaus,"ntaus/I");
  //
  // Control regions
  _tree->Branch("mu1pid", &_mu1pid,"mu1pid/I");  
  _tree->Branch("mu1pt", &_mu1pt,"mu1pt/D");  
  _tree->Branch("mu1p", &_mu1p,"mu1p/D");  
  _tree->Branch("mu1eta", &_mu1eta,"mu1eta/D");  
  _tree->Branch("mu1phi", &_mu1phi,"mu1phi/D");  
  _tree->Branch("mu2pid", &_mu2pid,"mu2pid/I");  
  _tree->Branch("mu2pt", &_mu2pt,"mu2pt/D");  
  _tree->Branch("mu2p", &_mu2p,"mu2p/D");  
  _tree->Branch("mu2eta", &_mu2eta,"mu2eta/D");  
  _tree->Branch("mu2phi", &_mu2phi,"mu2phi/D");  
  //
  _tree->Branch("wmt", &_wmt,"wmt/D");  
  _tree->Branch("zmass", &_zmass,"zmass/D");  
  _tree->Branch("zpt", &_zpt,"zpt/D");  
  _tree->Branch("zeta", &_zeta,"zeta/D");  
  _tree->Branch("zphi", &_zphi,"zphi/D");  

}


MiniAodEff::~MiniAodEff()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MiniAodEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Initialize //
  Init();

  // Get collections //

  // Trigger
  edm::Handle<edm::TriggerResults> H_trg_bits;
  iEvent.getByToken(trgBitsToken_, H_trg_bits);

  edm::Handle<pat::TriggerObjectStandAloneCollection> H_trg_obj;
  iEvent.getByToken(trgObjectsToken_, H_trg_obj);

  edm::Handle<pat::PackedTriggerPrescales> H_trg_ps;
  iEvent.getByToken(trgPrescalesToken_, H_trg_ps);

  // Objects
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

  edm::Handle<edm::TriggerResults> H_metfilt;
  iEvent.getByToken(metFiltToken_, H_metfilt);

  // Check validity
  if(!H_trg_bits.isValid()) {
    if(_verbose>0) cout << "Missing collection : " << _IT_trg_bits << " ... skip entry !" << endl;
    return;
  }

  if(!H_trg_obj.isValid()) {
    if(_verbose>0) cout << "Missing collection : " << _IT_trg_obj << " ... skip entry !" << endl;
    return;
  }

  if(!H_trg_ps.isValid()) {
    if(_verbose>0) cout << "Missing collection : " << _IT_trg_ps << " ... skip entry !" << endl;
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

  //  if(!H_metfilt.isValid()) {

    //return;
  //}
  ////////////////////////////////////////////////////////////////////////////////////////////

  // GLOBAL EVENT INFORMATIONS //
  _nRun   = iEvent.id().run();
  _nLumi  = iEvent.luminosityBlock();
  _nEvent = iEvent.id().event();

  // TRIGGER RESULTS //
  _trig_pass = "";
  _trig_n    = 0;
  TString path="";
  //
  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*H_trg_bits);
  for (int iHLT = 0 ; iHLT<static_cast<int>(H_trg_bits->size()); ++iHLT) {	
    if (H_trg_bits->accept (iHLT)) {
      path = TString(triggerNames.triggerName(iHLT));
      _trig_pass += "_%_"+path ;
      _trig_n++ ;
    }
  }

  // TRIGGER OBJECTS //
  string trgColl,trgFiltStr,trgPathsFFStr;
  vector<int> trgIds;
  vector<string> trgFilt, trgPathsFF, trgPathsFT, trgPathsTF, trgPathsTT;
  bool isNone, isL3, isLF, isBoth;
  isNone = isL3 = isLF = isBoth = false;

  if(_verbose>1) {
    cout << "$$$ SIZE OF H_trg_obj = " 
	 << H_trg_obj->size()
	 << " $$$" << endl;
  }

  for (pat::TriggerObjectStandAlone obj : *H_trg_obj) { // note: not "const &" since we want to call unpackPathNames

    obj.unpackPathNames(triggerNames);

    // pt,eta,phi
    if(_verbose>1) cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << endl;
    //
    _trig_obj_pt.push_back( obj.pt());
    _trig_obj_eta.push_back(obj.eta());
    _trig_obj_phi.push_back(obj.phi());

    // Collection
    trgColl = obj.collection();
    if(_verbose>1) cout << "\t   Collection: " << trgColl << endl;
    //
    _trig_obj_col.push_back(trgColl);

    // Trigger Filter Ids
    if(_verbose>1) cout << "\t   Type IDs:   ";
    trgIds = obj.filterIds();
    //
    _trig_obj_ids.push_back(trgIds);
    //
    for (unsigned h = 0; h < trgIds.size(); ++h) {
      if(_verbose>1) cout << " " << trgIds[h] ;
    }
    if(_verbose>1) cout << endl;

    // Trigger filters
    if(_verbose>1) cout << "\t   Filters:    ";
    trgFilt = obj.filterLabels();
    for (unsigned h = 0; h < trgFilt.size(); ++h) {
      if(_verbose>1) cout << " " << trgFilt[h];
      trgFiltStr += trgFilt[h]+"_%_" ;
    }
    if(_verbose>1) cout << endl;
    _trig_obj_lab.push_back(trgFiltStr);

    // Trigger paths
    trgPathsFF = obj.pathNames(false,false);
    trgPathsFT = obj.pathNames(false,true);
    trgPathsTF = obj.pathNames(true,false);
    trgPathsTT = obj.pathNames(true,true);
    //
    _trig_obj_path_FT.push_back(concatenate(trgPathsFT));
    _trig_obj_path_TF.push_back(concatenate(trgPathsTF));
    _trig_obj_path_TT.push_back(concatenate(trgPathsTT));
    //
    if(_verbose>1) cout << "\t   Paths (" << trgPathsFF.size()<<"/"<<trgPathsTT.size()<<"):    ";
    //

    // Loop over all associated paths
    for (unsigned h = 0, n = trgPathsFF.size(); h < n; ++h) {

      trgPathsFFStr += trgPathsFF[h]+"_%_";      

      isNone = obj.hasPathName( trgPathsFF[h], false, false ); 
      isL3   = obj.hasPathName( trgPathsFF[h], false, true ); 
      isLF   = obj.hasPathName( trgPathsFF[h], true, false ); 
      isBoth = obj.hasPathName( trgPathsFF[h], true, true ); 

      if(_verbose>1) {
	cout << "   " << trgPathsFF[h];
	if (isBoth)  cout << "(L,3)";
	if (isL3 && !isBoth)  cout << "(*,3)";
	if (isLF && !isBoth)  cout << "(L,*)";
	if (isNone && !isBoth && !isL3 && !isLF)  cout << "(*,*)";
      }
    }
    if(_verbose>1) cout << endl;
    //
    _trig_obj_path_FF.push_back(trgPathsFFStr);

    // Clear vectors
    trgIds.clear();
    trgFilt.clear(); 
    trgFiltStr="";
    trgColl="";
    trgPathsFFStr="";

    // Increment trigger object index
    _trig_obj_n++ ;
  }  

  /////////////////////

  // VERTICES //

  int vtx_counter=0;
  _vtx_N = H_vert->size();
  _vtx_N_stored = _nV;
	
  // select the primary vertex as the one with higest sum of (pt)^2 of tracks                                                                               
  PrimaryVertexSorter PVSorter;
  std::vector<reco::Vertex> sortedVertices = PVSorter.sortedList( *(H_vert.product()) );

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

  /////////////////////

  // JET/MET/MHT //

  // MET Filters
  if(_verbose>1) cout << "MET Filters : " ;
  edm::TriggerNames metFilters;
  if(H_metfilt.isValid()) {
    metFilters = iEvent.triggerNames(*H_metfilt);
    for (int iHLT = 0 ; iHLT<static_cast<int>(H_metfilt->size()); ++iHLT) {	
      if (H_metfilt->accept (iHLT)) {
	_met_filters.push_back(metFilters.triggerName(iHLT));
	if(_verbose>1) cout << "   " << metFilters.triggerName(iHLT);
      }
    }
  }
  else if(_verbose>0) {
    cout << "Missing collection : " << _IT_met_filt << " ... skip entry !" << endl;  
  }
  
  // MET 4-vector
  const pat::METCollection *C_pfmet = H_met.product();
  _METP4 = (*C_pfmet)[0].p4();

  // JETS AND MHT
  int nj_ht = 0, nj_mht = 0;
  double ht = 0., mhx = 0., mhy = 0.;
  double pt, eta, phi, px, py;
  //
  _jet_N = H_jets->size();
  //
  for (pat::JetCollection::const_iterator theJet = H_jets->begin(); theJet != H_jets->end(); ++theJet){

    // Jets
    //
    // Kinematics
    _jet_pt.push_back  ( theJet->pt() );
    _jet_eta.push_back ( theJet->eta() );
    _jet_phi.push_back ( theJet->phi() );
    _jet_e.push_back   ( theJet->energy() );
    _jet_m.push_back   ( theJet->mass() );
    //
    // Energy fractions
    _jet_efrac_ne_Had.push_back( theJet->neutralHadronEnergyFraction() );
    _jet_efrac_ne_EM.push_back ( theJet->neutralEmEnergyFraction() );
    _jet_efrac_ch_Had.push_back( theJet->chargedHadronEnergyFraction() );
    _jet_efrac_ch_EM.push_back ( theJet->chargedEmEnergyFraction() );
    _jet_efrac_ch_Mu.push_back ( theJet->chargedMuEnergyFraction() );
    //
    // Multiplicities
    _jet_mult_ch.push_back ( theJet->chargedMultiplicity() );
    _jet_mult_mu.push_back ( theJet->muonMultiplicity() );
    _jet_mult_ne.push_back ( theJet->neutralMultiplicity() );

    // MHT //
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

  // 4-vectors
  _METNoMuP4 = LV(mex_nomu, mey_nomu, 0, sqrt(mex_nomu*mex_nomu + mey_nomu*mey_nomu));
  _MHTNoMuP4 = LV(mhx_nomu, mhy_nomu, 0, sqrt(mhx_nomu*mhx_nomu + mhy_nomu*mhy_nomu));
  _MHTP4     = LV(mhx, mhy, 0, sqrt(mhx*mhx + mhy*mhy));

  // Flat values
  _met     = _METP4.Et();
  _met_phi = _METP4.Phi();
  _met_dphi= computeDeltaPhi( _met_phi , _jet_phi[0] );
  //
  _metnomu     = _METNoMuP4.Et();
  _metnomu_phi = _METNoMuP4.Phi();
  _metnomu_dphi= computeDeltaPhi( _metnomu_phi , _jet_phi[0] );
  //
  _mht     = _MHTP4.Et();
  _mht_phi = _MHTP4.Phi();
  _mht_dphi= computeDeltaPhi( _mht_phi , _jet_phi[0] );
  //
  _mhtnomu     = _MHTNoMuP4.Et();
  _mhtnomu_phi = _MHTNoMuP4.Phi();
  _mhtnomu_dphi= computeDeltaPhi( _mhtnomu_phi , _jet_phi[0] );

  /////////////////////

  // LEPTONS //
  _nphotons   = H_pho->size();
  _nelectrons = H_ele->size();
  _nmuons     = H_mu->size();

  //
  _ntaus = 0;
  for (pat::TauCollection::const_iterator taus_iter = H_tau->begin(); taus_iter != H_tau->end(); ++taus_iter) {
    if( taus_iter->pt() <= 20) continue; 
    if( fabs(taus_iter->eta()) >= 2.3) continue;
    if( taus_iter->tauID("decayModeFinding") <= 0.5) continue;
    if( taus_iter->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") <= 0.5) continue;
    _ntaus++;
  }  

  // CTRL REGIONS (just in case) //
  //
  TLorentzVector mu1vec, mu2vec, zvec; 
  //
  if (_nmuons >= 1) {
    pat::MuonRef muon = H_mu[0];
    _mu1pid = muon->pdgId(); 
    _mu1pt  = muon->pt(); 
    _mu1p   = muon->p(); 
    _mu1eta = muon->eta(); 
    _mu1phi = muon->phi();

    /*
    for (std::size_t i = 0; i < tightmuons.size(); i++) {
      if (muon == tightmuons[i]) _mu1id = 1;
    }
    */
    
    _wmt = sqrt(2.0 * _mu1pt * _met * (1.0 - cos(computeDeltaPhi(_mu1phi, _met_phi))));

    if (_nmuons >= 2) {        
      pat::MuonRef muon = H_mu[1];
      _mu2pid = muon->pdgId(); 
      _mu2pt  = muon->pt(); 
      _mu2p   = muon->p(); 
      _mu2eta = muon->eta(); 
      _mu2phi = muon->phi();
      
      /*
	for (std::size_t i = 0; i < tightmuons.size(); i++) {
	if (muon == tightmuons[i]) mu2id = 1;
	}
      */
            
      mu1vec.SetPtEtaPhiE(_mu1pt, _mu1eta, _mu1phi, _mu1p);
      mu2vec.SetPtEtaPhiE(_mu2pt, _mu2eta, _mu2phi, _mu2p);
      zvec = mu1vec + mu2vec;
      
      _zmass = zvec.M();
      _zpt   = zvec.Pt();
      _zeta  = zvec.Eta();            
      _zphi  = zvec.Phi();
    }
  }

  /////////////////////

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

  _nEvent = _nRun = _nLumi = 0;

  // Trigger
  _trig_pass = "";
  _trig_n = 0;
  _trig_obj_n = 0;
  _trig_obj_pt.clear();
  _trig_obj_eta.clear();
  _trig_obj_phi.clear();
  _trig_obj_col.clear();
  _trig_obj_lab.clear();
  _trig_obj_ids.clear();
  _trig_obj_path_FF.clear();
  _trig_obj_path_FT.clear();
  _trig_obj_path_TF.clear();
  _trig_obj_path_TT.clear();

  // Vertices
  _vtx_N = _vtx_N_stored = 0; 
  for(UInt_t iv=0;iv<_nV;iv++) {
    _vtx_normalizedChi2[iv] = 0.;
    _vtx_ndof[iv] = 0.;
    _vtx_nTracks[iv] = 0.;
    _vtx_d0[iv] = 0.;
    _vtx_x[iv] = 0.;
    _vtx_y[iv] = 0.;
    _vtx_z[iv] = 0.;
  }

  // MET
  _met = _mht = _metnomu = _mhtnomu = 
    _met_phi = _mht_phi = _metnomu_phi = _mhtnomu_phi = 
    _met_dphi = _mht_dphi = _metnomu_dphi = _mhtnomu_dphi = 0;
  _met_filters.clear();
  
  // Jets
  _jet_N=0;
  _jet_eta.clear() ;
  _jet_phi.clear() ;
  _jet_pt.clear()  ;
  _jet_e.clear()   ;
  _jet_m.clear()   ;
  _jet_mult_ch.clear() ;
  _jet_mult_mu.clear() ;
  _jet_mult_ne.clear() ;
  _jet_efrac_ne_Had.clear() ;
  _jet_efrac_ne_EM.clear() ;
  _jet_efrac_ch_Had.clear() ; 
  _jet_efrac_ch_EM.clear() ; 
  _jet_efrac_ch_Mu.clear() ;

  // Leptons/Photons
  _nphotons = _nelectrons = _nmuons = _ntaus = 0;


}

Double_t MiniAodEff::computeDeltaPhi(double phi1, double phi2)
{
  // Return value in [0;pi] 
  /*
  float dphi0 = TMath::Abs( phi1 - phi2 );
  if(dphi0 > TMath::Pi()) return  TMath::TwoPi() - dphi0;
  else                    return dphi0;
  */
  return TMath::Abs( TVector2::Phi_mpi_pi(phi1 - phi2) );
}

string MiniAodEff::concatenate(vector<string> vstring)
{
  string result;
  for(UInt_t i=0 ; i<vstring.size() ; i++) {
    result += vstring[i]+"_%_";
  }
  return result;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAodEff);
