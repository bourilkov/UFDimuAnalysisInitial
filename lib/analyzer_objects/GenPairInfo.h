
#ifndef GENPAIR_INFO
#define GENPAIR_INFO

#include <vector>
#include "TMath.h"
#include "ParticleInfo.h"

struct GenPairInfo : public ParticleInfo{

  Int_t iMu1;
  Int_t iMu2;
  Int_t mother_ID;
  Int_t postFSR;

  Double_t mass;
  Double_t pt  ;
  Double_t eta ;
  Double_t y   ;
  Double_t phi ;
  Double_t angle;

  void init();
  Float_t getMass();
  TLorentzVector get4vec();
  TString outputInfo();
  Double_t iso();

  ClassDef(GenPairInfo, 2);

};

#endif  // #ifndef GENPAIR_INFO
