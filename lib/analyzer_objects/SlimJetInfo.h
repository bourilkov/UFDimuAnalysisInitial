
#ifndef SLIMJET_INFO
#define SLIMJET_INFO

#include <vector>
#include "TMath.h"
#include "ParticleInfo.h"

struct SlimJetInfo : public ParticleInfo{
  Float_t pt      ;
  Float_t eta     ;
  Float_t phi     ;
  Float_t mass    ;
  Int_t   partonID;

  // Jet Correction Factor--Above momentum is already corrected!!
  // This factor will return momentum to uncorrected value!!
  Float_t jecFactor;
  Float_t jecUnc   ;  // Jet Energy Correction Uncertainty

  Float_t CSV ;  // Btag CSV_v2
  Float_t puID;  // PUID

  void init();
  Float_t getMass();
  TLorentzVector get4vec();
  TString outputInfo();
  Double_t iso();
};


#endif  // #ifndef SLIMJET_INFO
