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
  _vertexCollection(pset.getParameter<edm::InputTag>("vertexCollection"))
{

  //now do what ever initialization is needed
  edm::Service<TFileService> fs ;
  _tree = fs->make <TTree>("HLTStudies","Efficiency");

  // Event
  _tree->Branch("nEvent",&_nEvent,"nEvent/I");
  _tree->Branch("nRun",&_nRun,"nRun/I");
  _tree->Branch("nLumi",&_nLumi,"nLumi/I");
  //
  // Vertices
  _tree->Branch("vtx_N",&_vtx_N,"vtx_N/I");
  _tree->Branch("vtx_normalizedChi2",&_vtx_normalizedChi2,"vtx_normalizedChi2[15]/D");
  _tree->Branch("vtx_ndof",&_vtx_ndof,"vtx_ndof[15]/D");
  _tree->Branch("vtx_nTracks",&_vtx_nTracks,"vtx_nTracks[15]/D");
  _tree->Branch("vtx_d0",&_vtx_d0,"vtx_d0[15]/D");
  _tree->Branch("vtx_x",&_vtx_x,"vtx_x[15]/D");
  _tree->Branch("vtx_y",&_vtx_y,"vtx_y[15]/D");
  _tree->Branch("vtx_z",&_vtx_z,"vtx_z[15]/D");
  //
  // Jets
  _tree->Branch("nJet",&_nJet,"nJet/I");
  //
  _tree->Branch("jet_eta",&_jet_eta,"jet_eta/D");
  _tree->Branch("jet_phi",&_jet_phi,"jet_phi/D");
  _tree->Branch("jet_pt",&_jet_pt,"jet_pt/D");
  _tree->Branch("jet_e",&_jet_e,"jet_e/D");
  _tree->Branch("jet_m",&_jet_m,"jet_m/D");
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

  edm::Handle<reco::VertexCollection> H_vert;
  iEvent.getByLabel(_vertexCollection, H_vert);

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

  // VERTICES //
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
    /*
    GlobalPoint local_vertexPosition(bs.position().x(),
				     bs.position().y(),
				     bs.position().z());
    */
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
  _trig_pass.clear();

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

//define this as a plug-in
DEFINE_FWK_MODULE(Efficiency);
