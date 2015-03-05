#ifndef DEF_FITEFFHLT_H
#define DEF_FITEFFHLT_H

// General C++
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

// RooFit headers
#include "RooGlobalFunc.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooEfficiency.h"
#include "RooDataSet.h"
#include "RooBinning.h"
#include "RooHist.h"
#include "RooWorkspace.h"

// Root headers
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFrame.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TString.h"

// Personal headers
#include "FuncCB.h"
#include "tdrstyle.h"

using namespace RooFit ;
using namespace std;

class FitEfficiency {

 public:

  FitEfficiency(TString, TString, TString);
  Int_t DefinePaths();
  Double_t ExtractThreshold(TString hltVar, TString path);
  Int_t PerformFit(TString path, TString cutdef);

 private:

  TChain* _tree;
  TString _dirResults;
  TString _base;
  vector<TString> _paths;
  
};

#endif
