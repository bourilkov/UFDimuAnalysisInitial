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

_TrackInfo ParticleTools::getGenMuDY(bool m, bool postFSR, VarSet& vars)
{
// m = muon 0 or muon 1, postFSR = want postFSR gen muon or preFSR gen muon?
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

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

bool ParticleTools::isMuPairGenMatchedDY(float dRcut, float ptCut, VarSet& vars)
{
// See if the candidate reco muon pair matches the pair 
    std::vector<TLorentzVector> recoMuons;
    std::vector<TLorentzVector> genMuons;

    TLorentzVector reco0, reco1, gen0, gen1;

    reco0.SetPtEtaPhiM(vars.recoMuons->at(0).pt,vars.recoMuons->at(0).eta, vars.recoMuons->at(0).phi, MASS_MUON);
    reco1.SetPtEtaPhiM(vars.recoMuons->at(1).pt,vars.recoMuons->at(1).eta, vars.recoMuons->at(1).phi, MASS_MUON);
    recoMuons.push_back(reco0);
    recoMuons.push_back(reco1);

   _TrackInfo tgen0 = getGenMuDY(0, 1, vars);
   _TrackInfo tgen1 = getGenMuDY(1, 1, vars);

   // no dy gen muon in event, obviously no match
   if(tgen0.pt < 0 || tgen1.pt < 0) return false;

   gen0.SetPtEtaPhiM(tgen0.pt, tgen0.eta, tgen0.phi, MASS_MUON);
   gen1.SetPtEtaPhiM(tgen1.pt, tgen1.eta, tgen1.phi, MASS_MUON);
   genMuons.push_back(gen0);
   genMuons.push_back(gen1);

   // binary number to keep track of the matches
   unsigned int matchResults = 0;
   for(unsigned int i=0; i<recoMuons.size(); i++)
   {
       for(unsigned int j=0; j<genMuons.size(); j++)
       {
           if(recoMuons[i].DeltaR(genMuons[j]) < dRcut && abs(recoMuons[i].Pt()-genMuons[i].Pt())/genMuons[i].Pt() < ptCut)
           {
               // 11, 10, 01, 00 -> bit 3, bit 2, bit 1, bit 0
               // reco1 matches gen1, reco1 matches gen0, reco0 matches gen1, reco0 matches gen0
               unsigned int bit = (i << 1) | (j << 0);   // figure out which bit should turn on 
               matchResults = matchResults | (1 << bit); // turn on the bit
           }
       }
   }
   // match0 = bits 0 and 1, match 1 = bits 2 and 3
   unsigned int match0 = matchResults & 0x3;
   unsigned int match1 = (matchResults & 0xC) >> 2;

   // reco0 has a bit on aka matches a gen, so does match1
   if(match0 && match1)
   {
       // if they both match only one and the same gen particle then only one bit is on 
       // and 10 | 10 = 10, 01 | 01 = 01, everything else gives 11 == 0x3
       if((match0 | match1) == 0x3) return true; 
   }
   return false;
}
