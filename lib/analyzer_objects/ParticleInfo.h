
#ifndef PARTICLE_INFO
#define PARTICLE_INFO

#include <vector>
#include "TMath.h"
#include "TString.h"
#include "TLorentzVector.h"

struct ParticleInfo {

  virtual void init() = 0;
  virtual TLorentzVector get4vec() = 0;
  virtual TString outputInfo() = 0;
  virtual Double_t iso() = 0;
  ClassDef(ParticleInfo, 1)
};

#endif                                                                                                               
