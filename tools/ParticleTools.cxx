// ParticleTools.cxx

#include <vector>
#include <iostream>
#include <cstdlib>

#include "ParticleTools.h"
#include "TLorentzVector.h"
#include "VarSet.h"

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

TLorentzVector ParticleTools::getMotherPtEtaPhiM(float pt0, float eta0, float phi0, float m0, float pt1, float eta1, float phi1, float m1)
{
    TLorentzVector p0, p1, mother;
    p0.SetPtEtaPhiM(pt0, eta0, phi0, m0);
    p1.SetPtEtaPhiM(pt1, eta1, phi1, m1);
    mother = p0 + p1;
    return mother;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

_TrackInfo ParticleTools::getGenMuDY(bool m, bool postFSR, VarSet& vars)
{
// m tells us which muon to get
    _TrackInfo genZmu = vars.genM1ZpostFSR;
    _TrackInfo genGmu = vars.genM1GpostFSR;

    if(!m && !postFSR)
    {
        genZmu = vars.genM1ZpreFSR;
        genGmu = vars.genM1GpreFSR;
    }
    if(m && !postFSR)
    {
        genZmu = vars.genM2ZpreFSR;
        genGmu = vars.genM2GpreFSR;
    }
    if(m && postFSR)
    {
        genZmu = vars.genM2ZpostFSR;
        genGmu = vars.genM2GpostFSR;
    }

    if(genZmu.pt > 0)
        return genZmu;

    else if(genGmu.pt > 0)
        return genGmu;
    
    return genZmu;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

TLorentzVector ParticleTools::convertTrackMuTo4Vec(_TrackInfo& track)
{
    TLorentzVector mu;
    // if the track is valid then set up the 4vec, if the track is not valid just return a blank 4vec with mass = 0
    if(track.pt > 0) mu.SetPtEtaPhiM(track.pt, track.eta, track.phi, MASS_MUON);
    return mu;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

bool ParticleTools::isValid4Vec(TLorentzVector& v)
{
    if(v.Pt() == 0 && v.Eta() == 0 && v.Phi() == 0 && v.M() == 0) return false;
    return true;
}
