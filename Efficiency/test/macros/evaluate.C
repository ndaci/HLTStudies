Double_t evaluate(Double_t x, Double_t mid, Double_t steep, Double_t max) {
  return max / (1 + TMath::Exp(-steep*(x - mid)));                                                                             
}

