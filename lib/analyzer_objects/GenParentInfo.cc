#include "GenParentInfo.h"

void GenParentInfo::init() {

  ID         = -999;
  status     = -999;
  nDaughters = -999;

  pt     = -999;
  eta    = -999;
  phi    = -999;
  mass   = -999;
  charge = -999;

  daughter_1_ID     = -999;
  daughter_2_ID     = -999;
  daughter_1_status = -999;
  daughter_2_status = -999;
  daughter_1_idx    = -999;
  daughter_2_idx    = -999;
}
///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

Float_t GenParentInfo::getMass()
{
    return mass;
}

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

TLorentzVector GenParentInfo::get4vec()
{
    TLorentzVector v;
    v.SetPtEtaPhiM(pt, eta, phi, getMass());
    return v;
}

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

TString GenParentInfo::outputInfo()
{
    TString s = Form("id: %d, status: %d, pt: %7.3f, eta: %7.3f, phi: %7.3f, mass: %7.3f, d1_idx: %d, d2_idx: %d", 
                      ID, status, pt, eta, phi, mass, daughter_1_idx, daughter_2_idx);
    return s;
}

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

Double_t GenParentInfo::iso()
{
    return 0.0;
}
