#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <stdio.h>
#include <TEfficiency.h>

void superimpose(TString cut){
  TFile *normal = new TFile("results/res_M1_74X_V8_14e33_v7/res_M1_74X_V8_14e33_v7.root", "READ");
  TFile *HCAL0 = new TFile("results/res_M1_73X_V13d_14e33_v4_HCAL0/res_M1_73X_V13d_14e33_v4_HCAL0.root", "READ");
  TFile *HCAL2 = new TFile("results/res_M1_73X_V13d_14e33_v4_HCAL2/res_M1_73X_V13d_14e33_v4_HCAL2.root", "READ");
  
  TString plot = "th_PFMETNoMu_" + cut + "_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_NoiseCleaned_v1";
  
  TEfficiency *eff_normal = (TEfficiency*)normal->Get(plot);  
  TEfficiency *eff_HCAL0 = (TEfficiency*)HCAL0->Get(plot);  
  TEfficiency *eff_HCAL2 = (TEfficiency*)HCAL2->Get(plot);

  TCanvas *c = new TCanvas();
  
  eff_normal->SetLineColor(kBlue);
  eff_HCAL0->SetLineColor(kRed);
  eff_HCAL2->SetLineColor(kGreen);
  
  TF1 *eff_normal_fit = eff_normal->GetListOfFunctions()->FindObject("fit2");
  TF1 *eff_HCAL0_fit = eff_HCAL0->GetListOfFunctions()->FindObject("fit2");
  TF1 *eff_HCAL2_fit = eff_HCAL2->GetListOfFunctions()->FindObject("fit2");
  
  eff_normal_fit->SetLineColor(kBlue);
  eff_HCAL0_fit->SetLineColor(kRed);
  eff_HCAL2_fit->SetLineColor(kGreen);
  
  eff_normal->Draw();
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.25);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.1);
  gStyle->SetStatTextColor(4);
  c->Update();
  eff_HCAL0->Draw("SAME");
  gStyle->SetStatY(0.4);
  gStyle->SetStatTextColor(2);
  c->Update();
  eff_HCAL2->Draw("SAME");
  gStyle->SetStatY(0.55);
  gStyle->SetStatTextColor(3);
  c->Update();
  
  TLegend *leg = new TLegend(0.58,0.58,0.73,0.68);
  leg->SetFillColor(kWhite);
  leg->AddEntry(eff_normal_fit,"14e33","l");
  leg->AddEntry(eff_HCAL0_fit,"14e33_HCAL0","l");
  leg->AddEntry(eff_HCAL2_fit,"14e33_HCAL2","l");
  leg->Draw();
  c->Update();
  
  TString output = "superimposed_" + cut + ".png";
  c->Print(output,"png");
}
