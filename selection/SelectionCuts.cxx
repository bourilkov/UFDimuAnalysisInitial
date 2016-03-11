///////////////////////////////////////////////////////////////////////////
//                             SelectionCuts.cxx                         //
//=======================================================================//
//                                                                       //
//        Define the Event and Muon Selection cuts for the analysis.     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "SelectionCuts.h"
#include "TMath.h"

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________TightMuonID ___________________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

TightMuonIdCuts::TightMuonIdCuts()
{
// Standard values for the tight muon id cuts

    cNumTrackerLayers = 5;     // >
    cNumValidMuonHits = 0;     // >
    cNumValidPixelHits = 0;    // >
    cNumOfMatchedStations = 1; // >
    cNormChiSquare = 10;       // <
    cd0_PV = 0.2;              // <
    cdz_PV = 0.5;              // <
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TightMuonIdCuts::TightMuonIdCuts(int numTrackerLayers, int numValidMuonHits, int numValidPixelHits, int numOfMatchedStations, 
int normChiSquare, float d0_PV, float dz_PV)
{
// Custom values for the tight muon id cuts

    cNumTrackerLayers = numTrackerLayers;
    cNumValidMuonHits = numValidMuonHits;
    cNumValidPixelHits = numValidPixelHits;
    cNumOfMatchedStations = numOfMatchedStations;
    cNormChiSquare = normChiSquare;
    cd0_PV = d0_PV;
    cdz_PV = dz_PV;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool TightMuonIdCuts::evaluate(VarSet& vars)
{
    // Cuts on reco1
    // if the event fails a single cut return false
    if(!vars.reco1.isGlobal) return false;
    if(!vars.reco1.isPFMuon) return false;
    if(!(vars.reco1.numTrackerLayers > cNumTrackerLayers)) return false;
    if(!(vars.reco1.numValidMuonHits > cNumValidMuonHits)) return false;
    if(!(vars.reco1.numValidPixelHits > cNumValidPixelHits)) return false;
    if(!(vars.reco1.numOfMatchedStations > cNumOfMatchedStations)) return false;
    if(!(vars.reco1.normChiSquare < cNormChiSquare)) return false;
    if(!(TMath::Abs(vars.reco1.d0_PV) < cd0_PV)) return false;
    if(!(TMath::Abs(vars.reco1.dz_PV) < cdz_PV)) return false;

    // Cuts on reco2
    if(!vars.reco2.isGlobal) return false;
    if(!vars.reco2.isPFMuon) return false;
    if(!(vars.reco2.numTrackerLayers > cNumTrackerLayers)) return false;
    if(!(vars.reco2.numValidMuonHits > cNumValidMuonHits)) return false;
    if(!(vars.reco2.numValidPixelHits > cNumValidPixelHits)) return false;
    if(!(vars.reco2.numOfMatchedStations > cNumOfMatchedStations)) return false;
    if(!(vars.reco2.normChiSquare < cNormChiSquare)) return false;
    if(!(TMath::Abs(vars.reco2.d0_PV) < cd0_PV)) return false;
    if(!(TMath::Abs(vars.reco2.dz_PV) < cdz_PV)) return false;

    // The event passed all the cuts return true
    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString TightMuonIdCuts::string()
{
    return TString("Tight Muon ID");
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________SynchEventSelection____________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

SynchEventSelectionCuts::SynchEventSelectionCuts()
{
// Default values for the event selection for the synchronization exercise

    cDimuMassMin = 100;         // >
    cDimuMassMax = 110;         // <
    cTrigMuPtMin = 20;          // >
    cTrigMuEtaMax = 2.4;        // <
    cPVzMax = 24;               // < 
    cNDFpv = 4;                 // >
    cNPV = 0;                   // > 
    cNJets = 2;                 // <= 
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

SynchEventSelectionCuts::SynchEventSelectionCuts(float dimuMassMin, float dimuMassMax, float trigMuPtMin, float trigMuEtaMax, float PVzMax, int NDFpv, int NPV, int nJets)
{
// Default values for the event selection for the synchronization exercise

    cDimuMassMin = dimuMassMin;       // >
    cDimuMassMax = dimuMassMax;       // <
    cTrigMuPtMin = trigMuPtMin;       // >
    cTrigMuEtaMax = trigMuEtaMax;     // <
    cPVzMax = PVzMax;                 // < 
    cNDFpv = NDFpv;                   // >
    cNPV = NPV;                       // > 
    cNJets = nJets;                   // <=
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool SynchEventSelectionCuts::evaluate(VarSet& vars)
{
    // if the event fails a single cut return false
    // Require oppositely charged muons in the event
    if(!(vars.reco1.charge != vars.reco2.charge)) return false;

    // Look at a certain mass range for synchronization purposes
    if(!(vars.recoCandMass > cDimuMassMin && vars.recoCandMass < cDimuMassMax)) return false;

    // One muon in the pair must pass one of the HLT triggers. This muon have the appropriate pt and eta.
    // Should probably make this into a function so that we can look at a larger number of triggers without cluttering this too much.
    if(vars.reco1.isHltMatched[0] && vars.reco1.pt > cTrigMuPtMin && TMath::Abs(vars.reco1.eta) < cTrigMuEtaMax) ;
    else if(vars.reco1.isHltMatched[1] && vars.reco1.pt > cTrigMuPtMin && TMath::Abs(vars.reco1.eta) < cTrigMuEtaMax) ;
    else if(vars.reco2.isHltMatched[0] && vars.reco2.pt > cTrigMuPtMin && TMath::Abs(vars.reco2.eta) < cTrigMuEtaMax) ;
    else if(vars.reco2.isHltMatched[1] && vars.reco2.pt > cTrigMuPtMin && TMath::Abs(vars.reco2.eta) < cTrigMuEtaMax) ;
    else return false;

    // Vertex selection
    if(!passesVertexSelection(vars.vertices)) return false;

    // nJets selection
    if(!(vars.validJets.size() <= cNJets)) return false;

    // Passed all the cuts
    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool SynchEventSelectionCuts::passesVertexSelection(_VertexInfo& vertices) 
{
// For synchronization we have to check that at least one vertex meets certain criteria.

    if(!(vertices.nVertices > cNPV)) return false;

    bool passesZreq = false;
    bool passesNDFreq = false;
    for(int v=0; v < vertices.nVertices && v < 20; ++v)
    {   
        if(TMath::Abs(vertices.z[v]) < cPVzMax) passesZreq = true;
        if(vertices.ndf[v] > cNDFpv) passesNDFreq = true;
 
        // There is at least one vertex that passes the Z req and at least one that passes the NDF req 
        // Maybe the condition should be that there is at least one that passes the Z req AND NDF req
        if(passesZreq && passesNDFreq) return true;
    }   
    return false;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString SynchEventSelectionCuts::string()
{
    return TString("Synch Event Selection Cuts");
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________SynchMuonSelection_____________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

SynchMuonSelectionCuts::SynchMuonSelectionCuts()
{
// Default muon selection values for the synchronization exercise

    cMinPt = 10;       // >
    cMaxEta = 2.4;     // <
    cMaxRelIso = 0.25; // <
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

SynchMuonSelectionCuts::SynchMuonSelectionCuts(float minPt, float maxEta, float maxRelIso)
{
// Custom muon selection values for the synchronization exercise

    cMinPt = minPt;         // >
    cMaxEta = maxEta;       // <
    cMaxRelIso = maxRelIso; // <
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool SynchMuonSelectionCuts::evaluate(VarSet& vars)
{
    // if either muon fails the selection return false
    if(!(evaluate(vars.reco1, vars.rho))) return false;
    if(!(evaluate(vars.reco2, vars.rho))) return false;
    return true;
    // Some references for the iso cuts (from old code)
    //synchCuts.push_back("(reco1.trackIsoSumPt + max(0.0,reco1.ecalIso + reco1.hcalIso - rho*pi*0.09))/reco1.pt < 0.25");
    //synchCuts.push_back("(reco2.trackIsoSumPt + max(0.0,reco2.ecalIso + reco2.hcalIso - rho*pi*0.09))/reco2.pt < 0.25");
    //synchCuts.push_back("(reco1.sumChargedHadronPtR03 + max(0.0,reco1.sumNeutralHadronEtR03+reco1.sumPhotonEtR03 - 0.5*reco1.sumPUPtR03))/reco1.pt <= 0.12");
    //synchCuts.push_back("(reco2.sumChargedHadronPtR03 + max(0.0,reco2.sumNeutralHadronEtR03+reco2.sumPhotonEtR03 - 0.5*reco2.sumPUPtR03))/reco2.pt <= 0.12");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool SynchMuonSelectionCuts::evaluate(_MuonInfo& recoMu, float rho)
{
    // if the event fails a single cut return false
    // pt cuts
    if(!(recoMu.pt > cMinPt)) return false;

    // eta cuts
    if(!(TMath::Abs(recoMu.eta) < cMaxEta)) return false;

    // isolation cuts
    if(!((recoMu.trackIsoSumPt + TMath::Max(0.0,recoMu.ecalIso + recoMu.hcalIso - rho*TMath::Pi()*0.09))/recoMu.pt < cMaxRelIso)) return false;

    return true;

    // Some references for the iso cuts (from old code)
    //synchCuts.push_back("(reco1.trackIsoSumPt + max(0.0,reco1.ecalIso + reco1.hcalIso - rho*pi*0.09))/reco1.pt < 0.25");
    //synchCuts.push_back("(reco2.trackIsoSumPt + max(0.0,reco2.ecalIso + reco2.hcalIso - rho*pi*0.09))/reco2.pt < 0.25");
    //synchCuts.push_back("(reco1.sumChargedHadronPtR03 + max(0.0,reco1.sumNeutralHadronEtR03+reco1.sumPhotonEtR03 - 0.5*reco1.sumPUPtR03))/reco1.pt <= 0.12");
    //synchCuts.push_back("(reco2.sumChargedHadronPtR03 + max(0.0,reco2.sumNeutralHadronEtR03+reco2.sumPhotonEtR03 - 0.5*reco2.sumPUPtR03))/reco2.pt <= 0.12");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString SynchMuonSelectionCuts::string()
{
    return TString("Synch Muon Selection Cuts");
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________Run1EventSelection_____________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

Run1EventSelectionCuts::Run1EventSelectionCuts()
{
// Default values for the modified run 1 event selection

    cTrigMuPtMin = 20;          // >, originally 24 GeV in accordance with the IsoMu24_Eta2p1 trigger
                                // we are using IsoMu20 triggers so this has been reduced to 20
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Run1EventSelectionCuts::Run1EventSelectionCuts(float trigMuPtMin)
{
// Default values for the modified run 1 event selection

    cTrigMuPtMin = trigMuPtMin;      // >
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool Run1EventSelectionCuts::evaluate(VarSet& vars)
{
    // if the event fails a single cut return false
    // Require oppositely charged muons in the event
    if(!(vars.reco1.charge != vars.reco2.charge)) return false;

    // One muon in the pair must pass one of the HLT triggers. This muon have the appropriate pt and eta.
    // Should probably make this into a function so that we can look at a larger number of triggers without cluttering this too much.
    if(vars.reco1.isHltMatched[0] && vars.reco1.pt > cTrigMuPtMin) ;
    else if(vars.reco1.isHltMatched[1] && vars.reco1.pt > cTrigMuPtMin) ;
    else if(vars.reco2.isHltMatched[0] && vars.reco2.pt > cTrigMuPtMin) ;
    else if(vars.reco2.isHltMatched[1] && vars.reco2.pt > cTrigMuPtMin) ;
    else return false;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString Run1EventSelectionCuts::string()
{
    return TString("Run1 Event Selection Cuts");
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________Run1MuonSelection______________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

Run1MuonSelectionCuts::Run1MuonSelectionCuts()
{
// Default muon selection values for the synchronization exercise

    cMinPt = 10;       // >, originally 15
    cMaxEta = 2.4;     // <, originally 2.1
    cMaxRelIso = 0.12; // <
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Run1MuonSelectionCuts::Run1MuonSelectionCuts(float minPt, float maxEta, float maxRelIso)
{
// Custom muon selection values for the synchronization exercise

    cMinPt = minPt;         // >
    cMaxEta = maxEta;       // <
    cMaxRelIso = maxRelIso; // <
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool Run1MuonSelectionCuts::evaluate(VarSet& vars)
{
    // if either muon fails the selection return false
    if(!evaluate(vars.reco1)) return false;
    if(!evaluate(vars.reco2)) return false;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool Run1MuonSelectionCuts::evaluate(_MuonInfo& recoMu)
{
// Test a single muon 

    // if the event fails a single cut return false
    // pt cuts
    if(!(recoMu.pt > cMinPt)) return false;

    // eta cuts
    if(!(TMath::Abs(recoMu.eta) < cMaxEta)) return false;

    // isolation cuts
    if(!((recoMu.sumChargedHadronPtR03 + 
          TMath::Max(0.0,recoMu.sumNeutralHadronEtR03+recoMu.sumPhotonEtR03 - 0.5*recoMu.sumPUPtR03))/recoMu.pt <= cMaxRelIso)) return false;
    
    return true;

    // Some references for the iso cuts (from old code)
    //synchCuts.push_back("(reco1.trackIsoSumPt + max(0.0,reco1.ecalIso + reco1.hcalIso - rho*pi*0.09))/reco1.pt < 0.25");
    //synchCuts.push_back("(reco2.trackIsoSumPt + max(0.0,reco2.ecalIso + reco2.hcalIso - rho*pi*0.09))/reco2.pt < 0.25");
    //synchCuts.push_back("(reco1.sumChargedHadronPtR03 + max(0.0,reco1.sumNeutralHadronEtR03+reco1.sumPhotonEtR03 - 0.5*reco1.sumPUPtR03))/reco1.pt <= 0.12");
    //synchCuts.push_back("(reco2.sumChargedHadronPtR03 + max(0.0,reco2.sumNeutralHadronEtR03+reco2.sumPhotonEtR03 - 0.5*reco2.sumPUPtR03))/reco2.pt <= 0.12");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString Run1MuonSelectionCuts::string()
{
    return TString("Run1 Muon Selection Cuts");
}
