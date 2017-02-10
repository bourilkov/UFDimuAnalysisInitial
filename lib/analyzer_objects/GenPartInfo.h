
#ifndef GENPART_INFO
#define GENPART_INFO

#include <vector>
#include "TMath.h"
#include "ParticleInfo.h"

struct GenPartInfo : public ParticleInfo{

  Float_t charge;
  Float_t mass;
  Float_t pt;
  Float_t ptErr;
  Float_t eta;  // pseudo rapidity
  Float_t y;    // rapidity
  Float_t phi;  // phi

  void init();
  Float_t getMass();
  TLorentzVector get4vec();
  TString outputInfo();
  Double_t iso();

  ClassDef(GenPartInfo, 2);

};

#endif  // #ifndef GENPART_INFO
