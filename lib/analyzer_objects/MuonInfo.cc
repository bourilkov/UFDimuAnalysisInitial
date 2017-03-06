
#include "MuonInfo.h"

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

void MuonInfo::init() {


    isTracker    = -999;
    isStandAlone = -999;
    isGlobal     = -999;
    
    isTightID        = -999;
    isMediumID       = -999;
    isMediumID2016   = -999;
    isLooseID        = -999;
    
    charge = -999;
    pt     = -999;
    eta    = -999;
    phi    = -999;

    pt_trk    = -999;
    ptErr_trk = -999;
    eta_trk   = -999;
    phi_trk   = -999;
    
    pt_KaMu           = -999;
    ptErr_KaMu        = -999;
    pt_KaMu_sys_up    = -999;
    pt_KaMu_sys_down  = -999;
    pt_KaMu_clos_up   = -999;
    pt_KaMu_clos_down = -999;

    pt_Roch          = -999;
    pt_Roch_sys_up   = -999;
    pt_Roch_sys_down = -999;
     
    d0_BS = -999;
    dz_BS = -999;
    
    d0_PV = -999;
    dz_PV = -999;

    relIso = -999;
    
    trackIsoSumPt     = -999;
    trackIsoSumPtCorr = -999;
    hcalIso           = -999;
    ecalIso           = -999;
    relCombIso        = -999;
    
    isPF = -999;
    
    pt_PF  = -999;
    eta_PF = -999;
    phi_PF = -999;
    
    sumChargedHadronPtR03   = -999;
    sumChargedParticlePtR03 = -999;
    sumNeutralHadronEtR03   = -999;
    sumPhotonEtR03          = -999;
    sumPUPtR03              = -999;
    
    sumChargedHadronPtR04   = -999;
    sumChargedParticlePtR04 = -999;
    sumNeutralHadronEtR04   = -999;
    sumPhotonEtR04          = -999;
    sumPUPtR04              = -999;
    
    for (unsigned int iTrig = 0; iTrig < nTrig; iTrig++) {
      isHltMatched[iTrig] = -999;
      hltPt[iTrig]        = -999;
      hltEta[iTrig]       = -999;
      hltPhi[iTrig]       = -999;
    }

} // End void MuonInfo::init()

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

Float_t MuonInfo::getMass()
{
    Float_t mass = 0.105658367;     //GeV/c2
    return mass;
}

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

TLorentzVector MuonInfo::get4vec()
{
    TLorentzVector v;
    v.SetPtEtaPhiM(pt, eta, phi, getMass());
    return v;
}

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

TString MuonInfo::outputInfo()
{
    TString s = Form("pt: %7.3f, pt_PF: %7.3f, eta: %7.3f, phi: %7.3f, isoMu24: %d, isoTkMu24: %d, isTight: %d, isMedium: %d", 
                      pt, pt_PF, eta, phi, isHltMatched[2], isHltMatched[3], isTightID, isMediumID);
    return s;
}

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

Double_t MuonInfo::iso()
{
    return iso04();
}

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

Double_t MuonInfo::iso03()
{
    return (sumChargedHadronPtR03 + TMath::Max(0.0,sumNeutralHadronEtR03 + sumPhotonEtR03 - 0.5*sumPUPtR03))/pt;
}

///////////////////////////////////////////////////////////
//--------------------------------------------------------
///////////////////////////////////////////////////////////

Double_t MuonInfo::iso04()
{
    return (sumChargedHadronPtR04 + TMath::Max(0.0,sumNeutralHadronEtR04 + sumPhotonEtR04 - 0.5*sumPUPtR04))/pt;
}
