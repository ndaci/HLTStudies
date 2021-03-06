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
   _resultName = "res_"+TString(fOption);
   _outfile = new TFile("results/"+_resultName+"/"+_resultName+".root","recreate");
   _outfile->cd();

   // Initialize tree variables
   InitVar();
   
   // Prepare counters
   _nProcessed = 0;
   _nPass = 0;
   _nPassOR = 0;
   
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
  if(trig_pass->Contains("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120"))
    _nPass++ ;

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

  cout << "Processed " << _nProcessed << " events." << endl
       << "Passed " << _nPass << " events." << endl
       << "Passed OR " << _nPassOR << " events." << endl;
}

Int_t MakePlots::InitVar()
{

  return 0;
}


Int_t MakePlots::InitHistos()
{

  const UInt_t nV = 3;
  const UInt_t nC = 5;
  TString label[nV] = {"PFJet","PFMETNoMu","PFMHTNoMu"};
  TString ancut[nC] = {"nocut","jpt100","jpt150","jpt200","analysis"};
  TString name,title,xtitle,ytitle,cutname,maptag;

  for(UInt_t iV=0 ; iV<nV ; iV++) {
    for(UInt_t iC=0 ; iC<nC ; iC++) {
      for(UInt_t iP=0 ; iP<_paths.size() ; iP++) {
	//
	maptag = label[iV]+"_"+ancut[iC];
	cutname   = (iV==0) ? " p_{T}" : "" ;
	name  = "h_"+maptag+"_"+_paths[iP];
	title = "Eff vs "+label[iV]+cutname+" : "+_paths[iP]+" ("+ancut[iC]+")";
	xtitle= label[iV]+cutname+" [GeV]";
	ytitle= "Efficiency";
	//
	_mHltVarH1D[_paths[iP]][maptag] = new TH1D(name,title,50,0,500);
	_mHltVarH1D[_paths[iP]][maptag]->SetXTitle(xtitle);
	_mHltVarH1D[_paths[iP]][maptag]->SetYTitle(ytitle);
	//_mHltVarH1D[_paths[iP]][maptag]->GetYaxis()->SetRangeUser(0,1.1); //ND FIXME
	//_mHltVarH1D[_paths[iP]][maptag]->Sumw2();
      }
    }
  }

  return 0;
}

Int_t MakePlots::FillHistos() {

  _nProcessed++ ;

  //cout << *trig_pass << endl;
  const UInt_t nV = 3;
  const UInt_t nC = 5;
  TString  ancut[nC] = {"nocut","jpt100","jpt150","jpt200","analysis"};
  bool     check[nC] = {true, jet_pt[0]>100, jet_pt[0]>150 , jet_pt[0]>200 , true};

  bool     leadjet    = jet_pt[0]>100 && TMath::Abs(jet_eta[0])<2.4 ;
  bool     subleadjet = jet_pt[1]>30 && TMath::Abs(jet_eta[1])<4.5 && TMath::Abs(jet_phi[0]-jet_phi[1]);
  bool     monojet    = nJet==1 || ( nJet>1 && !subleadjet);
  bool     dijet      = nJet>1 && subleadjet;
  bool     jetveto    = nJet<3 || !(jet_pt[2]>30 && TMath::Abs(jet_eta[2])<4.5);
  check[4] = leadjet && (monojet || dijet) && jetveto;

  TString  label[nV] = {"PFJet","PFMETNoMu","PFMHTNoMu"};
  Double_t value[nV] = {jet_pt[0],metnomu,mhtnomu};
  TString maptag;

  for(UInt_t iP=0 ; iP<_paths.size() ; iP++) {
    //if( _paths[iP]=="all" || trig_pass->Contains(_paths[iP]) ) {
    if( FillNumerator( _paths[iP] ) ) {
      if(_paths[iP].Contains("14e33")) _nPassOR++ ;
      for(UInt_t iV=0 ; iV<nV ; iV++) {
	for(UInt_t iC=0 ; iC<nC ; iC++) {
	  maptag = label[iV]+"_"+ancut[iC];
	  if(check[iC]) _mHltVarH1D[_paths[iP]][maptag]->Fill(value[iV]);
	}
      }
    }
  }


  return 0;
}

Int_t MakePlots::EndHistos() 
{

  TString thePath,theVar;
  TH1D *hPass, *hTot;
  TEfficiency *pEff;
  Double_t nPass, nTot,globalEff;
  nPass=nTot=globalEff=0;

  // Loop over M_HLT_VAR_H1D : 1 path/entry
  for(_it_mHltVarH1D=_mHltVarH1D.begin() ; 
      _it_mHltVarH1D!=_mHltVarH1D.end()  ; 
      _it_mHltVarH1D++) {

    thePath  = _it_mHltVarH1D->first;
    _mVarH1D = _it_mHltVarH1D->second;

    if(thePath=="all") continue;

    // Loop over M_VAR_H1D : 1 x-axis variable / entry
    for(_it_mVarH1D = _mVarH1D.begin(); 
	_it_mVarH1D !=_mVarH1D.end()  ; 
	_it_mVarH1D++) {

      theVar = _it_mVarH1D->first;  // x-axis variable
      hPass  = _it_mVarH1D->second; // events passing thePath
      hTot   = _mHltVarH1D["all"][theVar]; // all events

      nPass = (Double_t)hPass->Integral();
      nTot  = (Double_t)hTot ->Integral();
      globalEff = nTot!=0 ? nPass/nTot : -1.0;
      TString s_globalEff = "#epsilon = "+TString(Form("%.1f",100*globalEff))+" %";

      if(theVar=="PFMETNoMu") cout << thePath << " : " << 100*globalEff << " %" << endl;

      if( hTot && TEfficiency::CheckConsistency(*hPass,*hTot) ) {
	pEff = new TEfficiency(*hPass,*hTot);
	pEff->SetNameTitle( "t"+TString(hPass->GetName()) , 
			    hPass->GetTitle() );
  
  TF1 *f1 = new TF1("fit",evaluate,50,400,5);
  f1->SetParName(0, "m0");
  f1->SetParName(1, "sigma");
  f1->SetParName(2, "alpha");
  f1->SetParName(3, "n");
  f1->SetParName(4, "norm");
  f1->SetParameter(0, 120);
  f1->SetParameter(1, 1);
  f1->SetParameter(2, 1);
  f1->SetParameter(3, 5);
  f1->SetParameter(4, 1);
//   f1->SetParLimits(1, 0.01, 50);
//   f1->SetParLimits(2, 0.01, 8);
//   f1->SetParLimits(3, 1.1, 35);
//   f1->SetParLimits(4, 0.6, 1);
  
  TF1 *f2 = new TF1("fit2",evaluate2,0,500,3);
  f2->SetParName(0, "midpoint");
  f2->SetParName(1, "steepness");
  f2->SetParName(2, "max");
  f2->SetParameter(0, 120);
  f2->SetParameter(1, 0.1);
  f2->SetParameter(2, 1);
  f2->SetParLimits(2, 0.8, 1);

  pEff->Fit(f2);
	pEff->Write();
	TCanvas c("c","c",0,0,600,600);
	pEff->Draw("AP");
  
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.4);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.1);

	TPaveText *pt2 = new TPaveText(0.58,0.15,0.85,0.22,"brNDC"); 
	pt2->SetLineColor(1);
	pt2->SetTextColor(1);
	pt2->SetTextFont(42);
	pt2->SetTextSize(0.03);
	pt2->SetFillColor(kWhite);
	pt2->SetShadowColor(kWhite);
	pt2->AddText(s_globalEff);
	pt2->Draw();

	c.Print("results/"+_resultName+"/"+TString(hPass->GetName())+".png","png");
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

Double_t MakePlots::evaluate(double *x, double *par)
{ 
  double m = x[0];
  double m0 = par[0];
  double sigma = par[1];
  double alpha = par[2];
  double n = par[3];
  double norm = par[4];
  
  const double sqrtPiOver2 = 1.2533141373; // sqrt(pi/2)
  const double sqrt2 = 1.4142135624;

  Double_t sig = fabs((Double_t) sigma);
  
  Double_t t = (m - m0)/sig ;
  
  if (alpha < 0)
    t = -t;

  Double_t absAlpha = fabs(alpha / sig);
  Double_t a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
  Double_t b = absAlpha - n/absAlpha;

  //cout << a << " " << b << endl;

  ////// Pour la crystal ball
  // if (t <= absAlpha){
  //   return norm * exp(-0.5*t*t);
  // }
  // else
  //   {
  //     return norm * a * TMath::Power(t-b,-n) ;
  //   }

  Double_t aireGauche = (1 + ApproxErf( absAlpha / sqrt2 )) * sqrtPiOver2 ;
  Double_t aireDroite = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
  Double_t aire = aireGauche + aireDroite;

  if ( t <= absAlpha ){
    return norm * (1 + ApproxErf( t / sqrt2 )) * sqrtPiOver2 / aire ;
  }
  else{
    return norm * (aireGauche +  a * (1/TMath::Power(t-b,n-1) - 1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / aire ;
  }
  
 } 

Double_t MakePlots::ApproxErf(Double_t arg)
{
  static const double erflim = 5.0;
  if( arg > erflim )
    return 1.0;
  if( arg < -erflim )
    return -1.0;
  
  return TMath::Erf(arg);
}

Double_t MakePlots::evaluate2(double *x, double *par)
{ 
  return par[2] / (1 + TMath::Exp(-par[1]*(x[0] - par[0])));
 } 

Bool_t MakePlots::FillNumerator( TString path )
{
  // Denominator histogram
  if(path.Contains("all")) return true;

  // OR phase space
  TString myOption = TString(fOption);

  if(path.Contains("14e33") && myOption.Contains("14e33")) {
    //if(trig_pass->Contains("HLT_CaloJet500_NoJetID")) return true;
    //if(trig_pass->Contains("HLT_PFJet500")) return true;
    if(trig_pass->Contains("HLT_CaloMET200_NoiseCleaned")) return true;
    //if(trig_pass->Contains("HLT_DiCentralPFJet70_PFMET120_NoiseCleaned")) return true;
    if(trig_pass->Contains("HLT_PFMET170_NoiseCleaned")) return true;
    if(trig_pass->Contains("HLT_PFMET120_PFMHT120_IDTight")) return true;
    if(trig_pass->Contains("HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight")) return true;
    if(trig_pass->Contains("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_NoiseCleaned")) return true;
  }

  else if(path.Contains("7e33") && myOption.Contains("7e33")) {
    //if(trig_pass->Contains("HLT_CaloJet500_NoJetID")) return true;
    //if(trig_pass->Contains("HLT_PFJet500")) return true;
    if(trig_pass->Contains("HLT_CaloMET200_NoiseCleaned")) return true;
    //if(trig_pass->Contains("HLT_DiCentralPFJet70_PFMET120_NoiseCleaned")) return true;
    if(trig_pass->Contains("HLT_PFMET170_NoiseCleaned")) return true;
    if(trig_pass->Contains("HLT_PFMET90_PFMHT90_IDTight")) return true;
    if(trig_pass->Contains("HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight")) return true;
    if(trig_pass->Contains("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_NoiseCleaned")) return true;
  }

  else if(path.Contains("5e33") && myOption.Contains("5e33")) {
    //if(trig_pass->Contains("HLT_CaloJet500_NoJetID")) return true;
    //if(trig_pass->Contains("HLT_PFJet500")) return true;
    if(trig_pass->Contains("HLT_CaloMET200_NoiseCleaned")) return true;
    //if(trig_pass->Contains("HLT_DiCentralPFJet70_PFMET120_NoiseCleaned")) return true;
    if(trig_pass->Contains("HLT_PFMET170_NoiseCleaned")) return true;
    if(trig_pass->Contains("HLT_PFMET90_PFMHT90_IDTight")) return true;
    if(trig_pass->Contains("HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight")) return true;
    if(trig_pass->Contains("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_NoiseCleaned")) return true;
  }

  /*
  if(path.Contains("7e33") && myOption.Contains("7e33")) {
    if(trig_pass->Contains("")) return true;
  }

  if(path.Contains("5e33") && myOption.Contains("5e33")) {
    if(trig_pass->Contains("")) return true;
  }
  */

  // Standard HLT path
  if(trig_pass->Contains(path)) return true;

  // default value
  return false;
}

Int_t MakePlots::DefinePaths()
{
  _paths.clear();
  _paths.push_back("all");
  _paths.push_back("OR_14e33");
  _paths.push_back("OR_7e33");
  _paths.push_back("OR_5e33");
  
  _paths.push_back("HLT_PFMET170_NoiseCleaned_v1");
  _paths.push_back("HLT_Overlap_PFMET170_CaloMET170_v1");
  _paths.push_back("HLT_Overlap_PFMET170_CaloMET180_v1");
  _paths.push_back("HLT_Overlap_PFMET170_CaloMET190_v1");
  _paths.push_back("HLT_Overlap_PFMET170_CaloMET200_v1");
  _paths.push_back("HLT_CaloMET170_NoiseCleaned_v1");
  _paths.push_back("HLT_CaloMET180_NoiseCleaned_v1");
  _paths.push_back("HLT_CaloMET190_NoiseCleaned_v1");
  _paths.push_back("HLT_CaloMET200_NoiseCleaned_v1");
  _paths.push_back("HLT_PFMET90_PFMHT90_IDTight_v1");
  _paths.push_back("HLT_PFMET120_PFMHT120_IDTight_v1");
  _paths.push_back("HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v1");
  _paths.push_back("HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v1");
  _paths.push_back("HLT_Overlap_MET90_MHT90_v1");
  _paths.push_back("HLT_Overlap_MET120_MHT120_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_NoiseCleaned_NoCalo_v1");
/*
  _paths.push_back("HLT_CaloJet500_NoJetID_v1");
  _paths.push_back("HLT_PFJet500_v1");
  _paths.push_back("HLT_CaloMET170_NoiseCleaned_v1");
  _paths.push_back("HLT_DiCentralPFJet70_PFMET120_NoiseCleaned_v1");
  _paths.push_back("HLT_PFMET90_PFMHT90_IDLoose_v1");
  _paths.push_back("HLT_PFMET100_PFMHT100_IDLoose_v1");
  _paths.push_back("HLT_PFMET110_PFMHT110_IDLoose_v1");
  _paths.push_back("HLT_PFMET120_PFMHT120_IDLoose_v1");
  _paths.push_back("HLT_PFMET170_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu100_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu140_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet150_PFMETNoMu150_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu90_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu80_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu80_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu90_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu150_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu150_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu160_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMHTNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet80_PFMETNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu100_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu110_PFMHTNoMu90_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu90_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu80_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu120_PFMHTNoMu80_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu110_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu120_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu100_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu130_PFMHTNoMu90_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu90_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu120_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu110_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu130_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu100_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu140_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu130_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu120_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu140_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu110_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu150_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu140_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu130_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu150_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu120_PFMHTNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu160_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMHTNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet90_PFMETNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu100_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu110_PFMHTNoMu90_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu90_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu80_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu120_PFMHTNoMu80_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu110_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu120_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu100_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu130_PFMHTNoMu90_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu90_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu120_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu110_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu130_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu100_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu140_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu130_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu120_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu140_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu110_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu150_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu140_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu130_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu150_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu120_PFMHTNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu160_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMHTNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet100_PFMETNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu100_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu110_PFMHTNoMu90_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu90_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu80_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu120_PFMHTNoMu80_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu110_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu120_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu100_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu130_PFMHTNoMu90_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu90_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu120_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu110_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu130_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu100_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu140_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu130_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu120_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu140_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu110_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu150_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu140_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu130_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu150_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu120_PFMHTNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu160_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMHTNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet110_PFMETNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu100_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu110_PFMHTNoMu90_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu90_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu80_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu120_PFMHTNoMu80_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu110_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu120_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu100_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu130_PFMHTNoMu90_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu90_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu120_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu110_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu130_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu100_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu140_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu130_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu120_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu140_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu110_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu150_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu140_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu130_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu150_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu120_PFMHTNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu160_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMHTNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet120_PFMETNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu100_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu110_PFMHTNoMu90_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu90_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu80_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu120_PFMHTNoMu80_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu110_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu120_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu100_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu130_PFMHTNoMu90_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu90_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu120_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu110_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu130_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu100_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu140_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu130_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu120_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu140_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu110_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu150_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu140_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu130_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu150_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu120_PFMHTNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu160_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMHTNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet130_PFMETNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu100_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu110_PFMHTNoMu90_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu90_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu80_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu120_PFMHTNoMu80_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu110_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu120_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu100_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu130_PFMHTNoMu90_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu90_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu120_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu110_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu130_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu100_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu140_PFMHTNoMu100_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu130_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu120_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu140_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu110_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu150_PFMHTNoMu110_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu140_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu130_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu150_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu120_PFMHTNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu160_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMHTNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMHTNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMHTNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMHTNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMHTNoMu160_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu120_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu130_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu140_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu150_NoiseCleaned_v1");
  _paths.push_back("HLT_MonoCentralPFJet140_PFMETNoMu160_NoiseCleaned_v1");
  */

  return 0;
}
