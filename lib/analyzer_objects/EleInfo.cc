
#include "EleInfo.h"

void EleInfo::init() {

  isPF               = -999;
  isTightID          = -999;
  isMediumID         = -999;
  isLooseID          = -999;
  isVetoID           = -999;
  passConversionVeto = -999;
  
  charge = -999;
  pt     = -999;
  eta    = -999;
  phi    = -999;
  
  d0_PV = -999;
  dz_PV = -999;
  
  missingInnerHits = -999;
  
  sumChargedHadronPtR03 = -999;
  sumNeutralHadronEtR03 = -999;
  sumPhotonEtR03        = -999;
  sumPUPtR03            = -999;

} // End void EleInfo::init()

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

Float_t EleInfo::getMass()
{
  Float_t mass = 0.0005109989461; //GeV/c2s
  return mass;
}

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

TLorentzVector EleInfo::get4vec()
{
    TLorentzVector v;
    v.SetPtEtaPhiM(pt, eta, phi, getMass());
    return v;
}

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////
    
TString EleInfo::outputInfo()
{
    TString s = Form("pt: %7.3f, eta: %7.3f, phi: %7.3f, isTight: %d, isMedium: %d", 
                      pt, eta, phi, isTightID, isMediumID);

    return s;
}

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

Double_t EleInfo::iso()
{
    return (sumChargedHadronPtR03 + TMath::Max(0.0,sumNeutralHadronEtR03 + sumPhotonEtR03 - 0.5*sumPUPtR03))/pt;
}

