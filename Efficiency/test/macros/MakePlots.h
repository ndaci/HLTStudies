#ifndef MakePlots_h
#define MakePlots_h

// Headers
#include <iostream>
#include <algorithm>
#include <map>
#include <vector>
//
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TEfficiency.h>
#include <TCanvas.h>
//
// Header file for the classes stored in the TTree if any.
#include <TString.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
typedef map<TString,TH1D*> M_VAR_H1D; // map var and histo
typedef map<TString,TH2D*> M_VAR_H2D;
typedef map<TString,TH3D*> M_VAR_H3D;
typedef map< TString , M_VAR_H1D > M_HLT_VAR_H1D; // map path and previous map
typedef map< TString , M_VAR_H2D > M_HLT_VAR_H2D;
typedef map< TString , M_VAR_H3D > M_HLT_VAR_H3D;

using namespace std;

class MakePlots : public TSelector {
 public :
  ClassDef(MakePlots,0);
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  
  // Declaration of leaf types
  Int_t           nEvent;
  Int_t           nRun;
  Int_t           nLumi;
  TString         *trig_pass;
  Int_t           trig_n;
  Int_t           vtx_N;
  Double_t        vtx_normalizedChi2[15];
  Double_t        vtx_ndof[15];
  Double_t        vtx_nTracks[15];
  Double_t        vtx_d0[15];
  Double_t        vtx_x[15];
  Double_t        vtx_y[15];
  Double_t        vtx_z[15];
  Double_t        met;
  Double_t        mht;
  Double_t        metnomu;
  Double_t        mhtnomu;
  Double_t        met_eta;
  Double_t        mht_eta;
  Double_t        metnomu_eta;
  Double_t        mhtnomu_eta;
  Double_t        met_phi;
  Double_t        mht_phi;
  Double_t        metnomu_phi;
  Double_t        mhtnomu_phi;
  Double_t        met_dphi;
  Double_t        mht_dphi;
  Double_t        metnomu_dphi;
  Double_t        mhtnomu_dphi;
  Int_t           nJet;
  Double_t        jet_eta[3];
  Double_t        jet_phi[3];
  Double_t        jet_pt[3];
  Double_t        jet_e[3];
  Double_t        jet_m[3];
  Int_t           jet_mult_ch[3];   //[nJet]
  Int_t           jet_mult_mu[3];   //[nJet]
  Int_t           jet_mult_ne[3];   //[nJet]
  Double_t        jet_efrac_ne_Had[3];   //[nJet]
  Double_t        jet_efrac_ne_EM[3];   //[nJet]
  Double_t        jet_efrac_ch_Had[3];   //[nJet]
  Double_t        jet_efrac_ch_EM[3];   //[nJet]
  Double_t        jet_efrac_ch_Mu[3];   //[nJet]

  // List of branches
  TBranch        *b_nEvent;   //!
  TBranch        *b_nRun;   //!
  TBranch        *b_nLumi;   //!
  TBranch        *b_trig_pass;   //!
  TBranch        *b_trig_n;   //!
  TBranch        *b_vtx_N;   //!
  TBranch        *b_vtx_normalizedChi2;   //!
  TBranch        *b_vtx_ndof;   //!
  TBranch        *b_vtx_nTracks;   //!
  TBranch        *b_vtx_d0;   //!
  TBranch        *b_vtx_x;   //!
  TBranch        *b_vtx_y;   //!
  TBranch        *b_vtx_z;   //!
  TBranch        *b_met;   //!
  TBranch        *b_mht;   //!
  TBranch        *b_metnomu;   //!
  TBranch        *b_mhtnomu;   //!
  TBranch        *b_met_eta;   //!
  TBranch        *b_mht_eta;   //!
  TBranch        *b_metnomu_eta;   //!
  TBranch        *b_mhtnomu_eta;   //!
  TBranch        *b_met_phi;   //!
  TBranch        *b_mht_phi;   //!
  TBranch        *b_metnomu_phi;   //!
  TBranch        *b_mhtnomu_phi;   //!
  TBranch        *b_met_dphi;   //!
  TBranch        *b_mht_dphi;   //!
  TBranch        *b_metnomu_dphi;   //!
  TBranch        *b_mhtnomu_dphi;   //!
  TBranch        *b_nJet;   //!
  TBranch        *b_jet_eta;   //!
  TBranch        *b_jet_phi;   //!
  TBranch        *b_jet_pt;   //!
  TBranch        *b_jet_e;   //!
  TBranch        *b_jet_m;   //!
  TBranch        *b_jet_mult_ch;   //!
  TBranch        *b_jet_mult_mu;   //!
  TBranch        *b_jet_mult_ne;   //!
  TBranch        *b_jet_efrac_ne_Had;   //!
  TBranch        *b_jet_efrac_ne_EM;   //!
  TBranch        *b_jet_efrac_ch_Had;   //!
  TBranch        *b_jet_efrac_ch_EM;   //!
  TBranch        *b_jet_efrac_ch_Mu;   //!

 MakePlots(TTree * /*tree*/ =0) : fChain(0) { }
  virtual ~MakePlots() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();

  // Initialize members
  Int_t InitVar();
  // Define paths to analyze
  Int_t DefinePaths();
  Bool_t FillNumerator(TString path);
  // Initialize ; fill ; rescale and write histograms
  Int_t InitHistos();
  Int_t FillHistos();
  Int_t EndHistos();

 private:

  vector<TString> _paths;
  Double_t _nProcessed, _nPass, _nPassOR;
  TFile* _outfile;
  TString _resultName;

  // Map of histograms and iterators
  M_HLT_VAR_H1D _mHltVarH1D;
  M_HLT_VAR_H2D _mHltVarH2D;
  M_HLT_VAR_H3D _mHltVarH3D;
  //
  M_VAR_H1D _mVarH1D;
  M_VAR_H2D _mVarH2D;
  M_VAR_H3D _mVarH3D;
  //
  M_HLT_VAR_H1D::iterator _it_mHltVarH1D;
  M_HLT_VAR_H2D::iterator _it_mHltVarH2D;
  M_HLT_VAR_H3D::iterator _it_mHltVarH3D;
  //
  M_VAR_H1D::iterator _it_mVarH1D;
  M_VAR_H2D::iterator _it_mVarH2D;
  M_VAR_H3D::iterator _it_mVarH3D;
  //
};

void MakePlots::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trig_pass = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nEvent", &nEvent, &b_nEvent);
   fChain->SetBranchAddress("nRun", &nRun, &b_nRun);
   fChain->SetBranchAddress("nLumi", &nLumi, &b_nLumi);
   fChain->SetBranchAddress("trig_pass", &trig_pass, &b_trig_pass);
   fChain->SetBranchAddress("trig_n", &trig_n, &b_trig_n);
   fChain->SetBranchAddress("vtx_N", &vtx_N, &b_vtx_N);
   fChain->SetBranchAddress("vtx_normalizedChi2", vtx_normalizedChi2, &b_vtx_normalizedChi2);
   fChain->SetBranchAddress("vtx_ndof", vtx_ndof, &b_vtx_ndof);
   fChain->SetBranchAddress("vtx_nTracks", vtx_nTracks, &b_vtx_nTracks);
   fChain->SetBranchAddress("vtx_d0", vtx_d0, &b_vtx_d0);
   fChain->SetBranchAddress("vtx_x", vtx_x, &b_vtx_x);
   fChain->SetBranchAddress("vtx_y", vtx_y, &b_vtx_y);
   fChain->SetBranchAddress("vtx_z", vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("mht", &mht, &b_mht);
   fChain->SetBranchAddress("metnomu", &metnomu, &b_metnomu);
   fChain->SetBranchAddress("mhtnomu", &mhtnomu, &b_mhtnomu);
   fChain->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
   fChain->SetBranchAddress("mht_eta", &mht_eta, &b_mht_eta);
   fChain->SetBranchAddress("metnomu_eta", &metnomu_eta, &b_metnomu_eta);
   fChain->SetBranchAddress("mhtnomu_eta", &mhtnomu_eta, &b_mhtnomu_eta);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("mht_phi", &mht_phi, &b_mht_phi);
   fChain->SetBranchAddress("metnomu_phi", &metnomu_phi, &b_metnomu_phi);
   fChain->SetBranchAddress("mhtnomu_phi", &mhtnomu_phi, &b_mhtnomu_phi);
   fChain->SetBranchAddress("met_dphi", &met_dphi, &b_met_dphi);
   fChain->SetBranchAddress("mht_dphi", &mht_dphi, &b_mht_dphi);
   fChain->SetBranchAddress("metnomu_dphi", &metnomu_dphi, &b_metnomu_dphi);
   fChain->SetBranchAddress("mhtnomu_dphi", &mhtnomu_dphi, &b_mhtnomu_dphi);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_e", &jet_e, &b_jet_e);
   fChain->SetBranchAddress("jet_m", &jet_m, &b_jet_m);
   fChain->SetBranchAddress("jet_mult_ch", jet_mult_ch, &b_jet_mult_ch);
   fChain->SetBranchAddress("jet_mult_mu", jet_mult_mu, &b_jet_mult_mu);
   fChain->SetBranchAddress("jet_mult_ne", jet_mult_ne, &b_jet_mult_ne);
   fChain->SetBranchAddress("jet_efrac_ne_Had", jet_efrac_ne_Had, &b_jet_efrac_ne_Had);
   fChain->SetBranchAddress("jet_efrac_ne_EM", jet_efrac_ne_EM, &b_jet_efrac_ne_EM);
   fChain->SetBranchAddress("jet_efrac_ch_Had", jet_efrac_ch_Had, &b_jet_efrac_ch_Had);
   fChain->SetBranchAddress("jet_efrac_ch_EM", jet_efrac_ch_EM, &b_jet_efrac_ch_EM);
   fChain->SetBranchAddress("jet_efrac_ch_Mu", jet_efrac_ch_Mu, &b_jet_efrac_ch_Mu);
}

Bool_t MakePlots::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif

// The class definition in MakePlots.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//

   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either MakePlots::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
