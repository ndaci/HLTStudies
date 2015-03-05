{
  TChain* ch = new TChain("eff/HLTStudies");
  ch->Add("../outeff.root");
  ch->Process("MakePlots.C++");
}
