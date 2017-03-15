
#ifndef GENMUPAIR_INFO
#define GENMUPAIR_INFO

#include <vector>
#include "TMath.h"
#include "ParticleInfo.h"

struct GenMuPairInfo : public ParticleInfo{

  Int_t iMu1;
  Int_t iMu2;
  Int_t mother_ID;
  Int_t postFSR;
  Int_t charge;

  Double_t mass ;
  Double_t pt   ;
  Double_t eta  ;
  Double_t rapid;
  Double_t phi  ;
  Double_t dR   ;
  Double_t dEta ;
  Double_t dPhi ;

  void init();
  Float_t getMass();
  TLorentzVector get4vec();
  TString outputInfo();
  Double_t iso();

  ClassDef(GenMuPairInfo, 2);

};

#endif  // #ifndef GENMUPAIR_INFO
