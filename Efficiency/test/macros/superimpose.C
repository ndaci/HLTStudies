#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <stdio.h>
#include <TEfficiency.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>

void superimpose(TString cut="analysis"){
  //TFile *normal = new TFile("results/res_M1_74X_V8_14e33_v7/res_M1_74X_V8_14e33_v7.root", "READ");
  TFile *HCAL0 = new TFile("results/res_r740p9_GRun740V33_L1P_14e33_HCAL0/res_r740p9_GRun740V33_L1P_14e33_HCAL0.root", "READ");
  TFile *HCAL2 = new TFile("results/res_r740p9_GRun740V33_L1P_14e33_HCAL2/res_r740p9_GRun740V33_L1P_14e33_HCAL2.root", "READ");
  TFile *HCAL3 = new TFile("results/res_r740p9_GRun740V33_L1P_14e33_HCAL3/res_r740p9_GRun740V33_L1P_14e33_HCAL3.root", "READ");

  //TString plot = "th_PFMETNoMu_" + cut + "_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_NoiseCleaned_v1";
  TString plot = "th_PFMETNoMu_" + cut + "_HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v1";
  
  TEfficiency *eff_HCAL0 = (TEfficiency*)HCAL0->Get(plot);  
  TEfficiency *eff_HCAL2 = (TEfficiency*)HCAL2->Get(plot);
  TEfficiency *eff_HCAL3 = (TEfficiency*)HCAL3->Get(plot);  

  TCanvas *c = new TCanvas();
  
  eff_HCAL0->SetLineColor(kRed);
  eff_HCAL2->SetLineColor(kBlue);
  eff_HCAL3->SetLineColor(kGreen+2);
  
  TF1 *eff_HCAL0_fit = (TF1*) eff_HCAL0->GetListOfFunctions()->FindObject("fit2");
  TF1 *eff_HCAL2_fit = (TF1*) eff_HCAL2->GetListOfFunctions()->FindObject("fit2");
  TF1 *eff_HCAL3_fit = (TF1*) eff_HCAL3->GetListOfFunctions()->FindObject("fit2");
  
  eff_HCAL0_fit->SetLineColor(kRed);
  eff_HCAL2_fit->SetLineColor(kBlue);
  eff_HCAL3_fit->SetLineColor(kGreen+2);
  
  eff_HCAL0->Draw();
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.25);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.1);
  //gStyle->SetStatTextColor(4);
  gStyle->SetStatTextColor(kRed);
  c->Update();
  eff_HCAL2->Draw("SAME");
  gStyle->SetStatY(0.4);
  //gStyle->SetStatTextColor(3);
  gStyle->SetStatTextColor(kBlue);
  c->Update();
  eff_HCAL3->Draw("SAME");
  gStyle->SetStatY(0.55);
  //gStyle->SetStatTextColor(2);
  gStyle->SetStatTextColor(kGreen+2);
  c->Update();
  
  TLegend *leg = new TLegend(0.58,0.58,0.73,0.68);
  leg->SetFillColor(kWhite);
  leg->AddEntry(eff_HCAL3_fit,"14e33_HCAL3","l");
  leg->AddEntry(eff_HCAL2_fit,"14e33_HCAL2","l");
  leg->AddEntry(eff_HCAL0_fit,"14e33_HCAL0","l");
  leg->Draw();
  c->Update();
  
  TString output = "results/superimposed_" + cut + ".pdf";
  c->Print(output,"pdf");
}
