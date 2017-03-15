#ifndef JET_PAIR_INFO
#define JET_PAIR_INFO

#include <vector>
#include "TMath.h"
#include "ParticleInfo.h"

struct JetPairInfo : public ParticleInfo{

  Int_t iJet1;
  Int_t iJet2;

  Float_t mass;
  Float_t pt;
  Float_t eta;
  Float_t phi;
  Float_t dR;
  Float_t dEta;
  Float_t dPhi;

  void init();
  Float_t getMass();
  TLorentzVector get4vec();
  TString outputInfo();
  Double_t iso();

  ClassDef(JetPairInfo, 2)
};

#endif  // #ifndef JET_PAIR_INFO
