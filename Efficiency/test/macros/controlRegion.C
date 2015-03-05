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
#include <TPaveText.h>
#include <TString.h>
#include <TLegend.h>

Int_t controlRegion(TString tag)
{

  TChain* ch = new TChain("eff/HLTStudies");
  ch->Add("../outeff_"+tag);

  TH1F* h_met = new TH1F("h_met_"+tag,"h_met_"+tag,100,0,1000);
  TH1F* h_metnomu = new TH1F("h_metnomu_"+tag,"h_metnomu_"+tag,100,0,1000);
  h_met->SetLineColor(kBlue);
  h_metnomu->SetLineColor(kRed);

  TCanvas* c = new TCanvas("c","c",10,10,600,600);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetLogy();
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);

  ch->Draw("met >> h_met","");
  ch->Draw("metnomu >> h_metnomu","","sames");

  TLegend* leg = new TLegend(0.45,0.74,0.68,0.88,"MET","brNDC");
  leg->AddEntry(h_met,"PFMET");
  leg->AddEntry(h_metnomu,"PFMETNoMu");
  leg->Draw();

  c->Print("plot_met_"+tag+".pdf","pdf");

  return 0;
}
