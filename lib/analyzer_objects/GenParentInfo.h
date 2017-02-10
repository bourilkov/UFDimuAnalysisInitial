
#ifndef GENPARENT_INFO
#define GENPARENT_INFO

#include <vector>
#include "TMath.h"
#include "ParticleInfo.h"

struct GenParentInfo : public ParticleInfo{

  Int_t ID;
  Int_t status;
  Int_t nDaughters;

  Float_t pt;
  Float_t eta;
  Float_t phi;
  Float_t mass;
  Int_t charge;

  Int_t daughter_1_ID;
  Int_t daughter_2_ID;
  Int_t daughter_1_status;
  Int_t daughter_2_status;
  Int_t daughter_1_idx;
  Int_t daughter_2_idx;

  void init();
  Float_t getMass();
  TLorentzVector get4vec();
  TString outputInfo();
  Double_t iso();

  ClassDef(GenParentInfo, 2)
};

#endif  // #ifndef GENPARENT_INFO
