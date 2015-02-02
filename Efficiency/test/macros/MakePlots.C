#define MakePlots_cxx
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("MakePlots.C")
// Root > T->Process("MakePlots.C","some options")
// Root > T->Process("MakePlots.C+")
//

#include "MakePlots.h"

int verbose=1;

void MakePlots::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   // Prepare output file
   TString resultName = "results";
   _outfile = new TFile("results/"+resultName+".root","recreate");
   _outfile->cd();

   // Initialize tree variables
   InitVar();
   
   // Prepare counters
   _nProcessed = 0;
   
   // Prepare list of paths
   DefinePaths();

   // Prepare histograms
   if(verbose>0) cout << "- Initialize Histograms : do 1D, " ;
   InitHistos();
   if(verbose>0) cout << "done." << endl;

}

void MakePlots::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t MakePlots::Process(Long64_t entry)
{

  //if(_nProcessed>=5) return kTRUE;

  fChain->GetEntry(entry);

  FillHistos();

  return kTRUE;
}

void MakePlots::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void MakePlots::Terminate()
{
  
  // Rescale histograms to get efficiencies
  //Float_t scaleFactor = _nProcessed!=0 ? 1/_nProcessed : 1;
  EndHistos();

  // Write histograms on output file
  if(verbose>1) cout << "- Write outfile" << endl;
  _outfile->Write();
  //
  if(verbose>1) cout << "- Close outfile" << endl;
  _outfile->Close();

  cout << "Processed " << _nProcessed << " events." << endl;
}

Int_t MakePlots::InitVar()
{

  return 0;
}

Int_t MakePlots::DefinePaths()
{
  _paths.clear();
  _paths.push_back("all");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu100_PFMHTNoMu140");
  return 0;
}

Int_t MakePlots::InitHistos()
{

  const UInt_t nV = 3;
  TString label[nV] = {"PFJet","PFMETNoMu","PFMHTNoMu"};
  TString name,title,xtitle,ytitle,cut;

  for(UInt_t iV=0 ; iV<nV ; iV++) {
    for(UInt_t iP=0 ; iP<_paths.size() ; iP++) {
      //
      cut   = (iV==0) ? " p_{T}" : "" ;
      name  = "hEff_"  +label[iV]+"_"+_paths[iP];
      title = "Eff vs "+label[iV]+cut+" : "+_paths[iP];
      xtitle= label[iV]+cut+" [GeV]";
      ytitle= "Efficiency";
      //
      _mHltVarH1D[_paths[iP]][label[iV]] = new TH1D(name,title,50,0,500);
      _mHltVarH1D[_paths[iP]][label[iV]]->SetXTitle(xtitle);
      _mHltVarH1D[_paths[iP]][label[iV]]->SetYTitle(ytitle);
      _mHltVarH1D[_paths[iP]][label[iV]]->GetYaxis()->SetRangeUser(0,1.1);
      _mHltVarH1D[_paths[iP]][label[iV]]->Sumw2();
    }
  }

  return 0;
}

Int_t MakePlots::FillHistos() {

  _nProcessed++ ;

  //cout << *trig_pass << endl;
  const UInt_t nV = 3;
  TString  label[nV] = {"PFJet","PFMETNoMu","PFMHTNoMu"};
  Double_t value[nV] = {jet_pt[0],metnomu,mhtnomu};

  for(UInt_t iP=0 ; iP<_paths.size() ; iP++) {
    if( _paths[iP]=="all" || trig_pass->Contains(_paths[iP]) ) {
      for(UInt_t iV=0 ; iV<nV ; iV++) {
	_mHltVarH1D[_paths[iP]][label[iV]]->Fill(value[iV]);
      }
    }
  }


  return 0;
}

Int_t MakePlots::EndHistos() 
{

  TString thePath,theVar;

  for(_it_mHltVarH1D=_mHltVarH1D.begin() ; 
      _it_mHltVarH1D!=_mHltVarH1D.end()  ; 
      _it_mHltVarH1D++) {

    thePath  = _it_mHltVarH1D->first;
    _mVarH1D = _it_mHltVarH1D->second;

    for(_it_mVarH1D = _mVarH1D.begin(); 
	_it_mVarH1D !=_mVarH1D.end()  ; 
	_it_mVarH1D++) {

      theVar = _it_mVarH1D->first;
      if( _mHltVarH1D["all"][theVar] ) {
	_it_mVarH1D->second->Divide( _mHltVarH1D["all"][theVar] );
      }
    }
  }

  /*
  for(_it_hTH2D=_hTH2D.begin() ; _it_hTH2D!=_hTH2D.end() ; _it_hTH2D++) {

  }

  for(_it_hTH3D=_hTH3D.begin() ; _it_hTH3D!=_hTH3D.end() ; _it_hTH3D++) {

  }
  */
  
  return 0;
}

