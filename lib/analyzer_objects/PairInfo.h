
#ifndef PAIR_INFO
#define PAIR_INFO

#include <vector>
#include "TMath.h"

struct PairInfo {

  Int_t iMu1;
  Int_t iMu2;

  Double_t mass;
  Double_t pt  ;
  Double_t eta ;
  Double_t y   ;
  Double_t phi ;
  Double_t angle;

  Double_t mass_PF;
  Double_t pt_PF  ;

  // need eta_PF, y_PF, phi_PF

  Double_t mass_trk;
  Double_t pt_trk  ;

  Double_t mass_KaMu;
  Double_t pt_KaMu  ;
  Double_t mass_KaMu_clos_up;
  Double_t pt_KaMu_clos_up  ;
  Double_t mass_KaMu_clos_down;
  Double_t pt_KaMu_clos_down  ;
  Double_t mass_KaMu_sys_up;
  Double_t pt_KaMu_sys_up  ;
  Double_t mass_KaMu_sys_down;
  Double_t pt_KaMu_sys_down  ;

  Double_t mass_Roch;
  Double_t pt_Roch  ;
  Double_t mass_Roch_sys_up;
  Double_t pt_Roch_sys_up  ;
  Double_t mass_Roch_sys_down;
  Double_t pt_Roch_sys_down  ;

  void init();

};


#endif  // #ifndef PAIR_INFO
