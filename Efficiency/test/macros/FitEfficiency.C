#include "FitEfficiency.h"

FitEfficiency::FitEfficiency(TString fileIn, TString dirResults, TString base)
{
  DefinePaths();
  _tree = new TChain("eff/HLTStudies");
  _tree->Add(fileIn);
  _dirResults = dirResults;
  _base = base;
}

Int_t FitEfficiency::PerformFit(TString path, TString cutdef)
{

  // STYLE //
  gROOT->Reset();
  loadPresentationStyle();  
  gROOT->ForceStyle();

  // OUTPUT //
  TString name_image = _dirResults+"/"+_base+"_"+path;
  ofstream outlog(name_image+".txt", ios::out);

  // BINNING //
  //const int nbins = 14;
  //Double_t bins[nbins] = {50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250};
  //RooBinning binning = RooBinning(nbins-1, bins, "binning");
  RooBinning binning = RooBinning(50, 0, 500, "binning");

  const UInt_t nBinsEta=3;
  const UInt_t nBinsPt=4;
  TString etacut[nBinsEta]={"fabs(jet_eta[0])<2.4",  // all
			    "fabs(jet_eta[0])<1.3",  // Barrel
			    "fabs(jet_eta[0])>1.3 && fabs(jet_eta[0])<2.4"}; // Endcap

  TString ptcut[nBinsPt]={"jet_pt[0]>=0",
			  "jet_pt[0]<100",
			  "jet_pt[0]>=100 && jet_pt[0]<200",
			  "jet_pt[0]>=200"};

  vector<TString> cut;
  for(UInt_t iE=0 ; iE<nBinsEta ; iE++) {
    for(UInt_t iP=0 ; iP<nBinsPt ; iP++) {
      if(cutdef=="") cutdef=="nEvent>=0";
      cut.push_back(cutdef+" && "+etacut[iE]+" && "+ptcut[iP]);
    }
  }

  // INPUT DATA //

  // Input variables
  RooRealVar xaxis("metnomu","PFMETNoMu [GeV]",0,500);
  RooRealVar pfMhtNoMu("mhtnomu","PFMHTNoMu [GeV]",0,50000);

  // Trigger matching
  //RooCategory *category = new RooCategory("trig_pass.Contains('caca')","");
  RooCategory *category = new RooCategory("","");
  category->defineType("accept",1) ;
  category->defineType("reject",0) ;

  // Input dataset
  RooDataSet dataSet("data","data from tree",
		     RooArgSet(xaxis, *category, pfMhtNoMu),
		     Import(*_tree), Cut("mhtnomu>=0") );
  dataSet.Print();

  // PLOT EFFICIENCY //
  TString lumi ="20 fb";
  Double_t thres=100;
  Int_t color=kRed;
  Int_t style=kOpenCircle;
  RooPlot* frame = xaxis.frame(Bins(18000),Title("Fitted efficiency")) ;
  dataSet.plotOn(frame, Binning(binning), Efficiency(*category), MarkerColor(color), LineColor(color), MarkerStyle(style) );


  // My modifs below //
  /*
  // FIT //

  // Fit parameters
  RooRealVar norm("norm","N",0.95,0.6,1);
  RooRealVar alpha("alpha","#alpha",0.2,0.01,8);
  RooRealVar n("n","n",2,1.1,35);
  RooRealVar mean("mean","mean",10,5,30);
  mean.setVal(100);
  RooRealVar sigma("sigma","#sigma",0.23,0.01,5);

  // Fit function
  FuncCB cb("cb","Fit function (cb)",xaxis,mean,sigma,alpha,n,norm) ;
  RooEfficiency eff("eff","efficiency", cb, *category, "accept");
  */

  /*
  for(UInt_t iC=0 ; iC<cut.size() ; iC++) {
    RooDataSet dataSet("data","data from tree",
		       RooArgSet(xaxis, *category),
		       Import(*_tree), Cut(cut[iC]) );
    dataSet.Print();
  }
  */

  // Original code below //
  /*  
  // FIT //
  double fit_cuts_min = thres-1.5 ;
  double fit_cuts_max = 150;
  xaxis.setRange("interesting",fit_cuts_min,fit_cuts_max);

  outlog << "Fit characteristics :" << endl
	  << "Threshold : "          << thres << endl 
	  << "Fit Range : ["         << fit_cuts_min << "," << fit_cuts_max << "]" << endl 
	  << endl << endl;

  RooFitResult* roofitres_EB = new RooFitResult("roofitres_EB","roofitres_EB");
  RooFitResult* roofitres_EE = new RooFitResult("roofitres_EE","roofitres_EE");

  // Fit #1 //
  roofitres_EB = eff_EB.fitTo(dataSetEB,ConditionalObservables(xaxis),Range("interesting"),Minos(kTRUE),Warnings(kFALSE),NumCPU(nCPU),Save(kTRUE),SumW2Error(kTRUE));
  cb_EB.plotOn(frame,LineColor(color1),LineWidth(2));

  outlog << "<----------------- EB ----------------->" << endl
	  << "double res_mean="  << mean.getVal()   << "; "
	  << "double res_sigma=" << sigma.getVal()  << "; "
          << "double res_alpha=" << alpha.getVal()  << "; "
          << "double res_n="     << n.getVal()      << "; "
          << "double res_norm="  << norm.getVal()   << "; "
	  << endl
	  << "double err_mean="  << mean.getError()  << "; "
	  << "double err_sigma=" << sigma.getError() << "; "
          << "double err_alpha=" << alpha.getError() << "; "
          << "double err_n="     << n.getError()     << "; "
          << "double err_norm="  << norm.getErrorLo()<< "; "
	  << endl;

  // Fit #2 //
  roofitres_EE = eff_EE.fitTo(dataSetEE,ConditionalObservables(xaxis),Range("interesting"),Minos(kTRUE),Warnings(kFALSE),NumCPU(nCPU),Save(kTRUE),SumW2Error(kTRUE));
  cb_EE.plotOn(frame,LineColor(color2),LineWidth(2));

  outlog << "<----------------- EE ----------------->" << endl
	  << "double res_mean="  << mean.getVal()   << "; "
	  << "double res_sigma=" << sigma.getVal()  << "; "
          << "double res_alpha=" << alpha.getVal()  << "; "
          << "double res_n="     << n.getVal()      << "; "
          << "double res_norm="  << norm.getVal()   << "; "
	  << endl
	  << "double err_mean="  << mean.getError()  << "; "
	  << "double err_sigma=" << sigma.getError() << "; "
          << "double err_alpha=" << alpha.getError() << "; "
          << "double err_n="     << n.getError()     << "; "
          << "double err_norm="  << norm.getErrorLo()<< "; "
	  << endl;
  */

  ////////////////////////////  DRAWING PLOTS AND LEGENDS /////////////////////////////////
  TCanvas* ca = new TCanvas("ca","Trigger Efficiency") ;

  ca->SetGridx();
  ca->SetGridy();
  ca->cd();
  
  //gPad->SetLogx();
  gPad->SetObjectStat(1);

  frame->GetYaxis()->SetRangeUser(0,1.05);
  frame->GetXaxis()->SetRangeUser(1,150.);
  frame->GetYaxis()->SetTitle("Efficiency");
  frame->GetXaxis()->SetTitle("PFMETNoMu [GeV]");
  frame->Draw() ;

  TH1F *SCeta = new TH1F("SCeta","SCeta",50,-2.5,2.5);

  SCeta->SetLineColor(color) ;
  SCeta->SetMarkerColor(color);
  SCeta->SetMarkerStyle(style);

  TLegend *leg = new TLegend(0.40, 0.14, 0.63, 0.34, NULL, "brNDC");
  leg->SetLineColor(1);
  leg->SetTextColor(1);
  leg->SetTextFont(42);
  leg->SetTextSize(0.0244755);
  leg->SetShadowColor(kWhite);
  leg->SetFillColor(kWhite);  
  leg->AddEntry("NULL","monojet efficiency","h");
  //    entry->SetLineColor(1);
  //    entry->SetLineStyle(1);
  //    entry->SetLineWidth(1);
  //    entry->SetMarkerColor(1);
  //    entry->SetMarkerStyle(21);
  //    entry->SetMarkerSize(1);
  //    entry->SetTextFont(62);
  leg->AddEntry(SCeta,"Barrel","p");
  //leg->Draw(); // ND FIXME
  
  ostringstream ossi("");
  ossi << thres;
  TString tossi = ossi.str();

  leg = new TLegend(0.40, 0.30, 0.50, 0.50, NULL, "brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.0297203);
  leg->SetLineColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry("NULL","CMS Preliminary 2012 pp  #sqrt{s}=8 TeV","h");
  leg->AddEntry("NULL","#int L dt = "+lumi+"^{-1}","h");
  leg->AddEntry("NULL","Threshold : "+tossi+" GeV","h");
  leg->Draw();
  
  //ca->Print(name_image+".C","C");
  //ca->Print(name_image+".cxx","cxx");
  ca->Print(name_image+".png","png");
  //ca->Print(name_image+".gif","gif");
  //ca->Print(name_image+".pdf","pdf");
  //ca->Print(name_image+".ps","ps");

  /////////////////////////////
  // SAVE THE ROO FIT RESULT //
  /////////////////////////////
  /*
  RooWorkspace *w = new RooWorkspace("workspace","workspace") ;

  w->import(dataSetEB);
  
  w->import(*roofitres,"roofitres");

  cout << "CREATES WORKSPACE : " << endl;
  w->Print();
  
  w->writeToFile(name_image+"_fitres.root") ;
  */
  //gDirectory->Add(w) ;

  // DELETE POINTERS
  // int a=0;
//   cin >> a;
//   delete treeTnP; delete cut; delete frame; delete roofitres; delete roofitres_EE; delete ca; delete SCeta; delete SCeta_EE; delete leg; delete w;

  return 0;
}

Double_t FitEfficiency::ExtractThreshold(TString hltVar, TString path)
{

  return 100;
}

Int_t FitEfficiency::DefinePaths()
{

  _paths.push_back("");
  return 0;
}


