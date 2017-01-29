
#ifndef ELE_INFO
#define ELE_INFO

#include <vector>
#include "TMath.h"
#include "ParticleInfo.h"

struct EleInfo : public ParticleInfo{

  Bool_t isTightID         ;
  Bool_t isMediumID        ;
  Bool_t isLooseID         ;
  Bool_t isVetoID          ;
  Bool_t passConversionVeto;

  Int_t   charge;
  Float_t pt    ;
  Float_t eta   ;
  Float_t phi   ;

  Float_t d0_PV           ;
  Float_t dz_PV           ;
  Float_t missingInnerHits;
 
  Bool_t isPF; 

  Float_t relIso               ;
  Float_t sumChargedHadronPtR03;  // sum-pt of charged Hadron 
  Float_t sumNeutralHadronEtR03;  // sum pt of neutral hadrons
  Float_t sumPhotonEtR03       ;  // sum pt of PF photons
  Float_t sumPUPtR03           ;  // sum pt of charged Particles not from PV (for Pu corrections)

  void init();
  Float_t getMass();
  TLorentzVector get4vec();
  TString outputInfo();
  Double_t iso();

};


#endif  // #ifndef ELE_INFO                                                                                                                      
