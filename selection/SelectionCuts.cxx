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
#include <iostream>

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
    cutset.cuts = std::vector<CutInfo>(18, CutInfo());
    makeCutSet();
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
    cutset.cuts = std::vector<CutInfo>(18, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool TightMuonIdCuts::evaluate(VarSet& vars)
{
    // Make sure passed is true at first and the value is reset
    cutset.resetCuts();

    // Keep track of the value that was cut on
    cutset.cuts[0].value = vars.reco1.isGlobal;
    cutset.cuts[1].value = vars.reco1.isPFMuon;
    cutset.cuts[2].value = vars.reco1.numTrackerLayers;
    cutset.cuts[3].value = vars.reco1.numValidMuonHits;
    cutset.cuts[4].value = vars.reco1.numValidPixelHits;
    cutset.cuts[5].value = vars.reco1.numOfMatchedStations;
    cutset.cuts[6].value = vars.reco1.normChiSquare;
    cutset.cuts[7].value = TMath::Abs(vars.reco1.d0_PV);
    cutset.cuts[8].value = TMath::Abs(vars.reco1.dz_PV);

    cutset.cuts[9].value = vars.reco2.isGlobal;
    cutset.cuts[10].value = vars.reco2.isPFMuon;
    cutset.cuts[11].value = vars.reco2.numTrackerLayers;
    cutset.cuts[12].value = vars.reco2.numValidMuonHits;
    cutset.cuts[13].value = vars.reco2.numValidPixelHits;
    cutset.cuts[14].value = vars.reco2.numOfMatchedStations;
    cutset.cuts[15].value = vars.reco2.normChiSquare;
    cutset.cuts[16].value = TMath::Abs(vars.reco2.d0_PV);
    cutset.cuts[17].value = TMath::Abs(vars.reco2.dz_PV);


    // Cuts on reco1
    // if the event fails a single cut return false
    if(!vars.reco1.isGlobal && cutset.cuts[0].on)
    { 
        cutset.cuts[0].passed = false; 
        return false; 
    }
    if(!vars.reco1.isPFMuon && cutset.cuts[1].on)
    { 
        cutset.cuts[1].passed = false; 
        return false; 
    }
    if(!(vars.reco1.numTrackerLayers > cNumTrackerLayers) && cutset.cuts[2].on)
    { 
        cutset.cuts[2].passed = false; 
        return false; 
    }
    if(!(vars.reco1.numValidMuonHits > cNumValidMuonHits) && cutset.cuts[3].on)
    { 
        cutset.cuts[3].passed = false; 
        return false; 
    }
    if(!(vars.reco1.numValidPixelHits > cNumValidPixelHits) && cutset.cuts[4].on)
    { 
        cutset.cuts[4].passed = false; 
        return false; 
    }
    if(!(vars.reco1.numOfMatchedStations > cNumOfMatchedStations) && cutset.cuts[5].on)
    { 
        cutset.cuts[5].passed = false; 
        return false; 
    }
    if(!(vars.reco1.normChiSquare < cNormChiSquare) && cutset.cuts[6].on)
    { 
        cutset.cuts[6].passed = false; 
        return false; 
    }
    if(!(TMath::Abs(vars.reco1.d0_PV) < cd0_PV) && cutset.cuts[7].on)
    { 
        cutset.cuts[7].passed = false; 
        return false; 
    }
    if(!(TMath::Abs(vars.reco1.dz_PV) < cdz_PV) && cutset.cuts[8].on)
    { 
        cutset.cuts[8].passed = false; 
        return false; 
    }

    // Cuts on reco2
    if(!vars.reco2.isGlobal && cutset.cuts[9].on)
    { 
        cutset.cuts[9].passed = false; 
        return false;
    }
    if(!vars.reco2.isPFMuon && cutset.cuts[10].on)
    { 
        cutset.cuts[10].passed = false; 
        return false;
    }
    if(!(vars.reco2.numTrackerLayers > cNumTrackerLayers) && cutset.cuts[11].on)
    { 
        cutset.cuts[11].passed = false; 
        return false;
    }
    if(!(vars.reco2.numValidMuonHits > cNumValidMuonHits) && cutset.cuts[12].on)
    { 
        cutset.cuts[12].passed = false; 
        return false;
    }
    if(!(vars.reco2.numValidPixelHits > cNumValidPixelHits) && cutset.cuts[13].on)
    { 
        cutset.cuts[13].passed = false; 
        return false;
    }
    if(!(vars.reco2.numOfMatchedStations > cNumOfMatchedStations) && cutset.cuts[14].on)
    { 
        cutset.cuts[14].passed = false; 
        return false;
    }
    if(!(vars.reco2.normChiSquare < cNormChiSquare) && cutset.cuts[15].on)
    { 
        cutset.cuts[15].passed = false; 
        return false;
    }
    if(!(TMath::Abs(vars.reco2.d0_PV) < cd0_PV) && cutset.cuts[16].on)
    { 
        cutset.cuts[16].passed = false; 
        return false;
    }
    if(!(TMath::Abs(vars.reco2.dz_PV) < cdz_PV) && cutset.cuts[17].on)
    { 
        cutset.cuts[17].passed = false; 
        return false;
    }

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

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void TightMuonIdCuts::makeCutSet()
{

    //cNumTrackerLayers      
    //cNumValidMuonHits      
    //cNumValidPixelHits     
    //cNumOfMatchedStations 
    //cNormChiSquare       
    //cd0_PV               
    //cdz_PV               

    // Cuts on reco1
    cutset.cuts[0].name = "reco1.isGlobal";
    cutset.cuts[0].tstring = "reco1.isGlobal";
    cutset.cuts[0].bins = 2;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 2;

    cutset.cuts[1].name = "reco1.isPFMuon";
    cutset.cuts[1].tstring = "reco1.isPFMuon";
    cutset.cuts[1].bins = 2;
    cutset.cuts[1].min = 0;
    cutset.cuts[1].max = 2;

    cutset.cuts[2].name = "reco1.numTrackerLayers";
    cutset.cuts[2].tstring.Form("reco1.numTrackerLayers > %d", cNumTrackerLayers);
    cutset.cuts[2].bins = 20;
    cutset.cuts[2].min = 0;
    cutset.cuts[2].max = 20;

    cutset.cuts[3].name = "reco1.numValidMuonHits";
    cutset.cuts[3].tstring.Form("reco1.numValidMuonHits > %d", cNumValidMuonHits);
    cutset.cuts[3].bins = 20;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 20;

    cutset.cuts[4].name = "reco1.numValidPixelHits";
    cutset.cuts[4].tstring.Form("reco1.numValidPixelHits > %d", cNumValidPixelHits);
    cutset.cuts[4].bins = 20;
    cutset.cuts[4].min = 0;
    cutset.cuts[4].max = 20;

    cutset.cuts[5].name = "reco1.numOfMatchedStations";
    cutset.cuts[5].tstring.Form("reco1.numOfMatchedStations > %d", cNumOfMatchedStations);
    cutset.cuts[5].bins = 20;
    cutset.cuts[5].min = 0;
    cutset.cuts[5].max = 20;

    cutset.cuts[6].name = "reco1.normChiSquare";
    cutset.cuts[6].tstring.Form("reco1.normChiSquare < %d", cNormChiSquare);
    cutset.cuts[6].bins = 50;
    cutset.cuts[6].min = 0;
    cutset.cuts[6].max = 50;

    cutset.cuts[7].name = "reco1.d0_PV";
    cutset.cuts[7].tstring.Form("abs(reco1.d0_PV) < %4.2f", cd0_PV);
    cutset.cuts[7].bins = 100;
    cutset.cuts[7].min = 0;
    cutset.cuts[7].max = 5;

    cutset.cuts[8].name = "reco1.dz_PV";
    cutset.cuts[8].tstring.Form("abs(reco1.dz_PV) < %4.2f", cdz_PV);
    cutset.cuts[8].bins = 100;
    cutset.cuts[8].min = 0;
    cutset.cuts[8].max = 5;


    // Cuts on reco2
    cutset.cuts[9].name = "reco2.isGlobal";
    cutset.cuts[9].tstring = "reco2.isGlobal";
    cutset.cuts[9].bins = 2;
    cutset.cuts[9].min = 0;
    cutset.cuts[9].max = 2;

    cutset.cuts[10].name = "reco2.isPFMuon";
    cutset.cuts[10].tstring = "reco2.isPFMuon";
    cutset.cuts[10].bins = 2;
    cutset.cuts[10].min = 0;
    cutset.cuts[10].max = 2;

    cutset.cuts[11].name = "reco2.numTrackerLayers";
    cutset.cuts[11].tstring.Form("reco2.numTrackerLayers > %d", cNumTrackerLayers);
    cutset.cuts[11].bins = 20;
    cutset.cuts[11].min = 0;
    cutset.cuts[11].max = 20;

    cutset.cuts[12].name = "reco2.numValidMuonHits";
    cutset.cuts[12].tstring.Form("reco2.numValidMuonHits > %d", cNumValidMuonHits);
    cutset.cuts[12].bins = 20;
    cutset.cuts[12].min = 0;
    cutset.cuts[12].max = 20;

    cutset.cuts[13].name = "reco2.numValidPixelHits";
    cutset.cuts[13].tstring.Form("reco2.numValidPixelHits > %d", cNumValidPixelHits);
    cutset.cuts[13].bins = 20;
    cutset.cuts[13].min = 0;
    cutset.cuts[13].max = 20;

    cutset.cuts[14].name = "reco2.numOfMatchedStations";
    cutset.cuts[14].tstring.Form("reco2.numOfMatchedStations > %d", cNumOfMatchedStations);
    cutset.cuts[14].bins = 20;
    cutset.cuts[14].min = 0;
    cutset.cuts[14].max = 20;

    cutset.cuts[15].name = "reco2.normChiSquare";
    cutset.cuts[15].tstring.Form("reco2.normChiSquare < %d", cNormChiSquare);
    cutset.cuts[15].bins = 50;
    cutset.cuts[15].min = 0;
    cutset.cuts[15].max = 50;

    cutset.cuts[16].name = "reco2.d0_PV";
    cutset.cuts[16].tstring.Form("abs(reco2.d0_PV) < %4.2f", cd0_PV);
    cutset.cuts[16].bins = 100;
    cutset.cuts[16].min = 0;
    cutset.cuts[16].max = 5;

    cutset.cuts[17].name = "reco2.dz_PV";
    cutset.cuts[17].tstring.Form("abs(reco2.dz_PV) < %4.2f", cdz_PV);
    cutset.cuts[17].bins = 100;
    cutset.cuts[17].min = 0;
    cutset.cuts[17].max = 5;
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
    cutset.cuts = std::vector<CutInfo>(5, CutInfo());
    makeCutSet();
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
    cutset.cuts = std::vector<CutInfo>(5, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool SynchEventSelectionCuts::evaluate(VarSet& vars)
{
    // Make sure passed is true at first and the value is reset
    cutset.resetCuts();

    // keep track of the values that were cut on
    cutset.cuts[0].value = vars.reco1.charge != vars.reco2.charge;
    cutset.cuts[1].value = vars.recoCandMassPF;
    cutset.cuts[2].value = 1;
    cutset.cuts[3].value = 1;
    cutset.cuts[4].value = vars.validJets.size();

    // if the event fails a single cut return false
    // Require oppositely charged muons in the event
    if(!(vars.reco1.charge != vars.reco2.charge) && cutset.cuts[0].on){ cutset.cuts[0].passed = false; return false;}

    // Look at a certain mass range for synchronization purposes
    if(!(vars.recoCandMassPF > cDimuMassMin && vars.recoCandMassPF < cDimuMassMax) && cutset.cuts[1].on){ cutset.cuts[1].passed = false; return false;}

    // One muon in the pair must pass one of the HLT triggers. This muon have the appropriate pt and eta.
    // Should probably make this into a function so that we can look at a larger number of triggers without cluttering this too much.
    if(!cutset.cuts[2].on) ;
    else if(vars.reco1.isHltMatched[0] && vars.reco1.pt > cTrigMuPtMin && TMath::Abs(vars.reco1.eta) < cTrigMuEtaMax) ; // reco1 passes trigger0
    else if(vars.reco1.isHltMatched[1] && vars.reco1.pt > cTrigMuPtMin && TMath::Abs(vars.reco1.eta) < cTrigMuEtaMax) ; // reco1 passes trigger1
    else if(vars.reco2.isHltMatched[0] && vars.reco2.pt > cTrigMuPtMin && TMath::Abs(vars.reco2.eta) < cTrigMuEtaMax) ; // reco2 passes trigger0
    else if(vars.reco2.isHltMatched[1] && vars.reco2.pt > cTrigMuPtMin && TMath::Abs(vars.reco2.eta) < cTrigMuEtaMax) ; // reco2 passes trigger1
    else 
    {
        cutset.cuts[2].value = 0;
        cutset.cuts[2].passed = false;
        return false;                                                                                                   // no muon passed a trigger
    }

    // Vertex selection
    if(!passesVertexSelection(vars.vertices) && cutset.cuts[3].on)
    { 
        cutset.cuts[3].passed = false;
        cutset.cuts[3].value = 0;
        return false;
    }

    // nJets selection
    if(!(vars.validJets.size() <= cNJets) && cutset.cuts[4].on){ cutset.cuts[4].passed = false;  return false;}

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

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void SynchEventSelectionCuts::makeCutSet()
{
    cutset.cuts[0].name = "reco1.charge != reco2.charge";
    cutset.cuts[0].tstring = "reco1.charge != reco2.charge";
    cutset.cuts[0].bins = 2;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 2;

    cutset.cuts[1].name = "recoCandMassPF";
    cutset.cuts[1].tstring.Form("recoCandMassPF > %f && recoCandMassPF < %f", cDimuMassMin, cDimuMassMax);
    cutset.cuts[1].bins = 140;
    cutset.cuts[1].min = 60;
    cutset.cuts[1].max = 200;

    cutset.cuts[2].name = "passes HLT Trigger Selection";
    cutset.cuts[2].tstring.Form("recoMu.isHltMatched[0||1] && recoMu.pt > %5.2f && TMath::Abs(recoMu.eta) < %4.2f", cTrigMuPtMin, cTrigMuEtaMax);
    cutset.cuts[2].bins = 2;
    cutset.cuts[2].min = 0;
    cutset.cuts[2].max = 2;

    cutset.cuts[3].name = "passes Vertex Selection";
    cutset.cuts[3].tstring.Form("v in V s.t. abs(vertices.z[v]) < %4.2f && v in V s.t. vertices.ndf[v] > %d", cPVzMax, cNDFpv);
    cutset.cuts[3].bins = 2;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 2;

    cutset.cuts[4].name = "nValidJets";
    cutset.cuts[4].tstring.Form("validJets.size() <= %d", cNJets);
    cutset.cuts[4].bins = 10;
    cutset.cuts[4].min = 0;
    cutset.cuts[4].max = 10;
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
    cMaxRelIso = 0.10; // <
    cutset.cuts = std::vector<CutInfo>(6, CutInfo());
    makeCutSet();
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
    cutset.cuts = std::vector<CutInfo>(6, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool SynchMuonSelectionCuts::evaluate(VarSet& vars)
{
    // Make sure passed is true at first and the value is reset
    cutset.resetCuts();

    // keep track of the values that were cut on
    // reco1
    cutset.cuts[0].value = vars.reco1.pt;
    cutset.cuts[1].value = vars.reco1.eta;
    cutset.cuts[2].value = vars.reco1.trackIsoSumPt/vars.reco1.pt;

    // reco2
    cutset.cuts[3].value = vars.reco2.pt;
    cutset.cuts[4].value = vars.reco2.eta;
    cutset.cuts[5].value = vars.reco2.trackIsoSumPt/vars.reco2.pt;

    // if either muon fails the selection return false
    if(!(evaluate(vars.reco1, 1))) return false;
    if(!(evaluate(vars.reco2, 2))) return false;
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

bool SynchMuonSelectionCuts::evaluate(_MuonInfo& recoMu, int m)
{
    // if the event fails a single cut return false
    // pt cuts
    // only return false if the cut is activated for the specific muon m=1or2
    if(!(recoMu.pt > cMinPt) && ((m==1 && cutset.cuts[0].on) || (m==2 && cutset.cuts[3].on)) ) return false;

    // eta cuts
    if(!(TMath::Abs(recoMu.eta) < cMaxEta)  && ((m==1 && cutset.cuts[1].on) || (m==2 && cutset.cuts[4].on))) return false;

    // isolation cuts
    if(!(recoMu.trackIsoSumPt/recoMu.pt < cMaxRelIso)  && ((m==1 && cutset.cuts[2].on) || (m==2 && cutset.cuts[5].on))) return false;
    //if(!( (recoMu.sumChargedHadronPtR03 + max(0.0,recoMu.sumNeutralHadronEtR03+recoMu.sumPhotonEtR03 - 0.5*recoMu.sumPUPtR03))/recoMu.pt <= cMaxRelIso) );

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

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void SynchMuonSelectionCuts::makeCutSet()
{
    cutset.cuts[0].name = "reco1.pt";
    cutset.cuts[0].tstring.Form("reco1.pt > %4.1f", cMinPt);
    cutset.cuts[0].bins = 200;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 200;

    cutset.cuts[1].name = "reco1.eta";
    cutset.cuts[1].tstring.Form("TMath::Abs(reco1.eta) < %4.2f", cMaxEta);
    cutset.cuts[1].bins = 50;
    cutset.cuts[1].min = -3;
    cutset.cuts[1].max = 3;

    cutset.cuts[2].name = "reco1.iso";
    cutset.cuts[2].tstring.Form("reco1.trackIsoSumPt/reco1.pt < %4.2f", cMaxRelIso);
    cutset.cuts[2].bins = 100;
    cutset.cuts[2].min = 0;
    cutset.cuts[2].max = 5;

    cutset.cuts[3].name = "reco2.pt";
    cutset.cuts[3].tstring.Form("reco2.pt > %f", cMinPt);
    cutset.cuts[3].bins = 200;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 200;

    cutset.cuts[4].name = "reco2.eta";
    cutset.cuts[4].tstring.Form("TMath::Abs(reco2.eta) < %4.2f", cMaxEta);
    cutset.cuts[4].bins = 50;
    cutset.cuts[4].min = -3;
    cutset.cuts[4].max = 3;

    cutset.cuts[5].name = "reco2.iso";
    cutset.cuts[5].tstring.Form("reco2.trackIsoSumPt/reco2.pt < %4.2f", cMaxRelIso);
    cutset.cuts[5].bins = 100;
    cutset.cuts[5].min = 0;
    cutset.cuts[5].max = 5;
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
    cDimuMassMin = 122.5;
    cDimuMassMax = 127.5;
    cutset.cuts = std::vector<CutInfo>(4, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Run1EventSelectionCuts::Run1EventSelectionCuts(float trigMuPtMin, float dimuMassMin, float dimuMassMax)
{
// Default values for the modified run 1 event selection

    cTrigMuPtMin = trigMuPtMin;      // >
    cDimuMassMin = dimuMassMin;      // >
    cDimuMassMax = dimuMassMax;      // <
    cutset.cuts = std::vector<CutInfo>(4, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool Run1EventSelectionCuts::evaluate(VarSet& vars)
{
    // Make sure passed is true at first and the value is reset
    cutset.resetCuts();

    cutset.cuts[0].value = vars.reco1.charge != vars.reco2.charge;
    cutset.cuts[1].value = -1;
    cutset.cuts[2].value = vars.recoCandMassPF;
    cutset.cuts[3].value = vars.recoCandMassPF;

    // Set cuts[1].value correctly
    if(vars.reco1.isHltMatched[0]) cutset.cuts[1].value = TMath::Abs(vars.reco1.pt);
    else if(vars.reco1.isHltMatched[1]) cutset.cuts[1].value = TMath::Abs(vars.reco1.pt);
    else if(vars.reco2.isHltMatched[0]) cutset.cuts[1].value = TMath::Abs(vars.reco2.pt);
    else if(vars.reco2.isHltMatched[1]) cutset.cuts[1].value = TMath::Abs(vars.reco2.pt);
    else
    {
         cutset.cuts[1].value = -1;
    }

    // if the event fails a single cut return false
    // Require oppositely charged muons in the event
    if(!(vars.reco1.charge != vars.reco2.charge) && cutset.cuts[0].on) return false;

    // One muon in the pair must pass one of the HLT triggers. This muon have the appropriate pt and eta.
    // Should probably make this into a function so that we can look at a larger number of triggers without cluttering this too much.
    if(!cutset.cuts[1].on) ;
    else if(vars.reco1.isHltMatched[0] && TMath::Abs(vars.reco1.pt) > cTrigMuPtMin) ;
    else if(vars.reco1.isHltMatched[1] && TMath::Abs(vars.reco1.pt) > cTrigMuPtMin) ;
    else if(vars.reco2.isHltMatched[0] && TMath::Abs(vars.reco2.pt) > cTrigMuPtMin) ;
    else if(vars.reco2.isHltMatched[1] && TMath::Abs(vars.reco2.pt) > cTrigMuPtMin) ;
    else
    {
         return false;
    }

    if(!(vars.recoCandMassPF > cDimuMassMin) && cutset.cuts[2].on) return false;
    if(!(vars.recoCandMassPF < cDimuMassMax) && cutset.cuts[3].on) return false;

    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString Run1EventSelectionCuts::string()
{
    return TString("Run1 Event Selection Cuts");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void Run1EventSelectionCuts::makeCutSet()
{
    //cTrigMuPtMin; 
    //cDimuMassMin; 

    cutset.cuts[0].name = "reco1.charge != reco2.charge";
    cutset.cuts[0].tstring = "reco1.charge != reco2.charge";
    cutset.cuts[0].bins = 2;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 2;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = "trigMatchedRecoMu.pt";
    cutset.cuts[1].tstring.Form("reco1.isHltMatched[0||1] && reco1.pt > %5.2f", cTrigMuPtMin);
    cutset.cuts[1].bins = 201;
    cutset.cuts[1].min = -1;
    cutset.cuts[1].max = 200;
    cutset.cuts[1].cutvalue = &cTrigMuPtMin;
    cutset.cuts[1].ismin = true;

    cutset.cuts[2].name = "recoCandMassPF_Min";
    cutset.cuts[2].tstring.Form("recoCandMassPF > %5.2f", cDimuMassMin);
    cutset.cuts[2].bins = 150;
    cutset.cuts[2].min = 50;
    cutset.cuts[2].max = 200;
    cutset.cuts[2].cutvalue = &cDimuMassMin;
    cutset.cuts[2].ismin = true;

    cutset.cuts[3].name = "recoCandMassPF_Max";
    cutset.cuts[3].tstring.Form("recoCandMassPF < %5.2f", cDimuMassMax);
    cutset.cuts[3].bins = 150;
    cutset.cuts[3].min = 50;
    cutset.cuts[3].max = 200;
    cutset.cuts[3].cutvalue = &cDimuMassMin;
    cutset.cuts[3].ismin = false;

    //cutset.cuts[2].name = "passes Vertex Selection";
    //cutset.cuts[2].tstring.Form("v in V s.t. abs(vertices.z[v]) < %4.2f && v in V s.t. vertices.ndf[v] > %d", cPVzMax, cNDFpv);

    //cutset.cuts[3].name = "nValidJets";
    //cutset.cuts[3].tstring.Form("validJets.size() <= %d", cNJets);
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
    cutset.cuts = std::vector<CutInfo>(6, CutInfo());
    makeCutSet();
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
    cutset.cuts = std::vector<CutInfo>(6, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool Run1MuonSelectionCuts::evaluate(VarSet& vars)
{
    // Make sure passed is true at first and the value is reset
    cutset.resetCuts();

    // keep track of the values that were cut on
    // reco1
    cutset.cuts[0].value = TMath::Abs(vars.reco1.pt);
    cutset.cuts[1].value = TMath::Abs(vars.reco1.eta);
    cutset.cuts[2].value = (vars.reco1.sumChargedHadronPtR03 + TMath::Max(0.0,vars.reco1.sumNeutralHadronEtR03+vars.reco1.sumPhotonEtR03 - 0.5*vars.reco1.sumPUPtR03))/vars.reco1.pt;

    // reco2
    cutset.cuts[3].value = TMath::Abs(vars.reco2.pt);
    cutset.cuts[4].value = TMath::Abs(vars.reco2.eta);
    cutset.cuts[5].value = (vars.reco2.sumChargedHadronPtR03 + TMath::Max(0.0,vars.reco2.sumNeutralHadronEtR03+vars.reco2.sumPhotonEtR03 - 0.5*vars.reco2.sumPUPtR03))/vars.reco2.pt;

    // if either muon fails the selection return false
    if(!evaluate(vars.reco1, 1)) return false;
    if(!evaluate(vars.reco2, 2)) return false;
    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool Run1MuonSelectionCuts::evaluate(_MuonInfo& recoMu, int m)
{
// Test a single muon 

    // if the event fails a single cut return false
    // pt cuts
    // only return false if the cut is activated for the specific muon m=1or2
    if(!(TMath::Abs(recoMu.pt) > cMinPt) && ((m==1 && cutset.cuts[0].on) || (m==2 && cutset.cuts[3].on)) ) return false;

    // eta cuts
    if(!(TMath::Abs(recoMu.eta) < cMaxEta)  && ((m==1 && cutset.cuts[1].on) || (m==2 && cutset.cuts[4].on))) return false;

    // isolation cuts
    if(!((recoMu.sumChargedHadronPtR03 + 
          TMath::Max(0.0,recoMu.sumNeutralHadronEtR03+recoMu.sumPhotonEtR03 - 0.5*recoMu.sumPUPtR03))/recoMu.pt <= cMaxRelIso) && ((m==1 && cutset.cuts[2].on) || (m==2 && cutset.cuts[5].on))) return false;
    
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

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void Run1MuonSelectionCuts::makeCutSet()
{
    //cMinPt 
    //cMaxEta
    //cMaxRelIso

    // reco1 cuts
    cutset.cuts[0].name = "reco1.pt";
    cutset.cuts[0].tstring.Form("reco1.pt > %4.1f", cMinPt);
    cutset.cuts[0].bins = 200;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 200;
    cutset.cuts[0].cutvalue = &cMinPt;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = "reco1.eta";
    cutset.cuts[1].tstring.Form("TMath::Abs(reco1.eta) < %4.2f", cMaxEta);
    cutset.cuts[1].bins = 50;
    cutset.cuts[1].min = 0;
    cutset.cuts[1].max = 3;
    cutset.cuts[1].cutvalue = &cMaxEta;
    cutset.cuts[1].ismin = false;

    cutset.cuts[2].name = "reco1.iso";
    cutset.cuts[2].tstring.Form("reco1.iso < %4.2f", cMaxRelIso);
    cutset.cuts[2].bins = 50;
    cutset.cuts[2].min = 0;
    cutset.cuts[2].max = 1;
    cutset.cuts[2].cutvalue = &cMaxRelIso;
    cutset.cuts[2].ismin = false;

     // reco2 cuts
    cutset.cuts[3].name = "reco2.pt";
    cutset.cuts[3].tstring.Form("reco2.pt > %f", cMinPt);
    cutset.cuts[3].bins = 200;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 200;
    cutset.cuts[3].cutvalue = &cMinPt;
    cutset.cuts[3].ismin = true;

    cutset.cuts[4].name = "reco2.eta";
    cutset.cuts[4].tstring.Form("TMath::Abs(reco2.eta) < %4.2f", cMaxEta);
    cutset.cuts[4].bins = 50;
    cutset.cuts[4].min = 0;
    cutset.cuts[4].max = 3;
    cutset.cuts[4].cutvalue = &cMaxEta;
    cutset.cuts[4].ismin = false;

    cutset.cuts[5].name = "reco2.iso";
    cutset.cuts[5].tstring.Form("reco2.iso < %4.2f", cMaxRelIso);
    cutset.cuts[5].bins = 50;
    cutset.cuts[5].min = 0;
    cutset.cuts[5].max = 1;
    cutset.cuts[5].cutvalue = &cMaxRelIso;
    cutset.cuts[5].ismin = false;
}

