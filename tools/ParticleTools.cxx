///////////////////////////////////////////////////////////////////////////
//                           ParticleTools.cxx                           //
//=======================================================================//
//                                                                       //
//        Miscellaneous tools for particles: output 4vec info to terminal//
//        get the correct gen muons from a gen event, some other tools   //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <vector>
#include <iostream>
#include <cstdlib>

#include "ParticleTools.h"
#include "TLorentzVector.h"
#include "VarSet.h"

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////


TString ParticleTools::output4vecInfo(TLorentzVector& v)
{
    TString s = Form("pt: %7.3f, eta: %7.3f, phi: %7.3f, mass: %11.5f", v.Pt(), v.Eta(), v.Phi(), v.M());
    return s;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

TLorentzVector ParticleTools::getMotherPtEtaPhiM(float pt0, float eta0, float phi0, float m0, float pt1, float eta1, float phi1, float m1)
{
// get a 4vec for the mother particle based upon the daughters
    TLorentzVector p0, p1, mother;
    p0.SetPtEtaPhiM(pt0, eta0, phi0, m0);
    p1.SetPtEtaPhiM(pt1, eta1, phi1, m1);
    mother = p0 + p1;
    return mother;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

bool ParticleTools::isValid4Vec(TLorentzVector& v)
{
    if(v.Pt() == 0 && v.Eta() == 0 && v.Phi() == 0 && v.M() == 0) return false;
    return true;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

float ParticleTools::dR(float eta1, float phi1, float eta2, float phi2)
{
    // it's easiest to use the TLorentzVector class to account for the dPhi > pi problem
    TLorentzVector v1; 
    TLorentzVector v2; 
    v1.SetPtEtaPhiM(10, eta1, phi1, 0); // doesn't matter what pt and mass are we only care about the dR value
    v2.SetPtEtaPhiM(10, eta2, phi2, 0); // same thing
    return v1.DeltaR(v2);
}
