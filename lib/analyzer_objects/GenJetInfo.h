
#ifndef GENJET_INFO
#define GENJET_INFO

#include <vector>
#include "TMath.h"
#include "ParticleInfo.h"

struct GenJetInfo : public ParticleInfo{

  Float_t px    ;
  Float_t py    ;
  Float_t pz    ;
  Float_t pt    ;
  Float_t eta   ;
  Float_t phi   ;
  Float_t mass  ;
  Float_t charge;

  void init();
  Float_t getMass();
  TLorentzVector get4vec();
  TString outputInfo();
  Double_t iso();


  ClassDef(GenJetInfo, 2)
};


#endif  // #ifndef GENJET_INFO
