#include "Efficiency.h"

//
// constructors and destructor
//
Efficiency::Efficiency(const edm::ParameterSet& pset) :

  _hltProcessName(pset.getParameter<string>("hltProcessName")),
  _trigResultsLabel("TriggerResults", "", _hltProcessName),
  _genParticleLabel(pset.getParameter<string>("genParticleLabel")),
  _pfjetCollection(pset.getParameter<edm::InputTag>("pfjetCollection")),
  _pfmetCollection(pset.getParameter<edm::InputTag>("pfmetCollection")),
  _muCollection(pset.getParameter<edm::InputTag>("muCollection"))

{

  //now do what ever initialization is needed
  edm::Service<TFileService> fs ;
  _tree = fs->make <TTree>("HLTStudies","Efficiency");

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

//define this as a plug-in
DEFINE_FWK_MODULE(Efficiency);
