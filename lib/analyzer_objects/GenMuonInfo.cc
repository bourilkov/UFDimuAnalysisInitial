
#include "GenMuonInfo.h"

void GenMuonInfo::init() {

  status     = -999;
  nMothers   = -999;
  postFSR    = -999;

  pt     = -999;
  eta    = -999;
  phi    = -999;
  mass   = -999;
  charge = -999;

  FSR_pt     = -999;
  FSR_eta    = -999;
  FSR_phi    = -999;
  FSR_mass   = -999;

  mother_ID     = -999;
  mother_status = -999;
  mother_idx    = -999;

}
///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

Float_t GenMuonInfo::getMass()
{
    return mass;
}

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

TLorentzVector GenMuonInfo::get4vec()
{
    TLorentzVector v;
    v.SetPtEtaPhiM(pt, eta, phi, getMass());
    return v;
}

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

TString GenMuonInfo::outputInfo()
{
    TString s = Form("status: %d, pt: %7.3f, ptf: %7.3f, eta: %7.3f, phi: %7.3f, mass: %7.3f, postFSR: %d, m_id: %d, m_idx %d", 
                      status, pt, FSR_pt, eta, phi, mass, postFSR, mother_ID, mother_idx);
    return s;
}

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

Double_t GenMuonInfo::iso()
{
    return 0.0;
}

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

bool GenMuonInfo::operator%(GenMuonInfo& other)
{
    TLorentzVector left = this->get4vec();
    TLorentzVector right = other.get4vec();
    if(left.DeltaR(right) < 0.1 && (this->postFSR != other.postFSR)) return true;
    else return false;
}

