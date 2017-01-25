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
#include "ParticleTools.h"
#include "TMath.h"
#include <iostream>

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
    cTrigMuPtMin = 24;          // >
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
    cutset.cuts[0].value = vars.recoMuons.charge[0] != vars.recoMuons.charge[1];
    cutset.cuts[1].value = vars.dimuCand.recoCandMassPF;
    cutset.cuts[2].value = 1;
    cutset.cuts[3].value = 1;
    cutset.cuts[4].value = vars.validJets.size();

    // if the event fails a single cut return false
    // Require oppositely charged muons in the event
    if(!(vars.recoMuons.charge[0] != vars.recoMuons.charge[1]) && cutset.cuts[0].on){ cutset.cuts[0].passed = false; return false;}

    // Look at a certain mass range for synchronization purposes
    if(!(vars.dimuCand.recoCandMassPF > cDimuMassMin && vars.dimuCand.recoCandMassPF < cDimuMassMax) && cutset.cuts[1].on){ cutset.cuts[1].passed = false; return false;}

    // One muon in the pair must pass one of the HLT triggers. This muon have the appropriate pt and eta.
    // Should probably make this into a function so that we can look at a larger number of triggers without cluttering this too much.
    if(!cutset.cuts[2].on) ;
    else if(vars.recoMuons.isHltMatched[0][4] && vars.recoMuons.pt[0] > cTrigMuPtMin && TMath::Abs(vars.recoMuons.eta[0]) < cTrigMuEtaMax) ; // recoMuon0 passes trigger0
    else if(vars.recoMuons.isHltMatched[0][5] && vars.recoMuons.pt[0] > cTrigMuPtMin && TMath::Abs(vars.recoMuons.eta[0]) < cTrigMuEtaMax) ; // recoMuon0 passes trigger1
    else if(vars.recoMuons.isHltMatched[1][4] && vars.recoMuons.pt[1] > cTrigMuPtMin && TMath::Abs(vars.recoMuons.eta[1]) < cTrigMuEtaMax) ; // recoMuon1 passes trigger0
    else if(vars.recoMuons.isHltMatched[1][5] && vars.recoMuons.pt[1] > cTrigMuPtMin && TMath::Abs(vars.recoMuons.eta[1]) < cTrigMuEtaMax) ; // recoMuon1 passes trigger1
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
    cutset.cuts[0].name = "recoMuons.charge[0] != recoMuons.charge[1]";
    cutset.cuts[0].tstring = "recoMuons.charge[0] != recoMuons.charge[1]";
    cutset.cuts[0].bins = 2;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 2;

    cutset.cuts[1].name = "dimuCand.recoCandMassPF";
    cutset.cuts[1].tstring.Form("dimuCand.recoCandMassPF > %f && dimuCand.recoCandMassPF < %f", cDimuMassMin, cDimuMassMax);
    cutset.cuts[1].bins = 140;
    cutset.cuts[1].min = 60;
    cutset.cuts[1].max = 200;

    cutset.cuts[2].name = "passes HLT Trigger Selection";
    cutset.cuts[2].tstring.Form("recoMu.isHltMatched[4||5] && recoMu.pt > %5.2f && TMath::Abs(recoMu.eta) < %4.2f", cTrigMuPtMin, cTrigMuEtaMax);
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
    // recoMuon0
    cutset.cuts[0].value = vars.recoMuons.pt[0];
    cutset.cuts[1].value = vars.recoMuons.eta[0];
    cutset.cuts[2].value = vars.recoMuons.trackIsoSumPt[0]/vars.recoMuons.pt[0];

    // recoMuon1
    cutset.cuts[3].value = vars.recoMuons.pt[1];
    cutset.cuts[4].value = vars.recoMuons.eta[1];
    cutset.cuts[5].value = vars.recoMuons.trackIsoSumPt[1]/vars.recoMuons.pt[1];

    // if either muon fails the selection return false
    if(!(evaluate(vars.recoMuons, 0))) return false;
    if(!(evaluate(vars.recoMuons, 1))) return false;
    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool SynchMuonSelectionCuts::evaluate(_MuonInfo& recoMuons, int m)
{
    // if the event fails a single cut return false
    // only return false if the cut is activated for the specific muon m=1or2
    if(!(recoMuons.pt[m] > cMinPt) && ((m==0 && cutset.cuts[0].on) || (m==1 && cutset.cuts[3].on)) ) return false;

    // eta cuts
    if(!(TMath::Abs(recoMuons.eta[m]) < cMaxEta)  && ((m==0 && cutset.cuts[1].on) || (m==1 && cutset.cuts[4].on))) return false;

    // isolation cuts
    if(!(recoMuons.trackIsoSumPt[m]/recoMuons.pt[m] < cMaxRelIso)  && ((m==0 && cutset.cuts[2].on) || (m==1 && cutset.cuts[5].on))) return false;
    //if(!( (recoMuons.sumChargedHadronPtR03 + max(0.0,recoMuons.sumNeutralHadronEtR03+recoMuons.sumPhotonEtR03 - 0.5*recoMuons.sumPUPtR03))/recoMuons.pt <= cMaxRelIso) );

    return true;
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
    cutset.cuts[0].name = "recoMuons.pt[0]";
    cutset.cuts[0].tstring.Form("recoMuons.pt[0] > %4.1f", cMinPt);
    cutset.cuts[0].bins = 200;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 200;

    cutset.cuts[1].name = "recoMuons.eta[0]";
    cutset.cuts[1].tstring.Form("TMath::Abs(recoMuons.eta[0]) < %4.2f", cMaxEta);
    cutset.cuts[1].bins = 50;
    cutset.cuts[1].min = -3;
    cutset.cuts[1].max = 3;

    cutset.cuts[2].name = "recoMuons.iso[0]";
    cutset.cuts[2].tstring.Form("recoMuons.trackIsoSumPt[0]/recoMuons.pt[0] < %4.2f", cMaxRelIso);
    cutset.cuts[2].bins = 100;
    cutset.cuts[2].min = 0;
    cutset.cuts[2].max = 5;

    cutset.cuts[3].name = "recoMuons.pt[1]";
    cutset.cuts[3].tstring.Form("recoMuons.pt[1] > %f", cMinPt);
    cutset.cuts[3].bins = 200;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 200;

    cutset.cuts[4].name = "recoMuons.eta[1]";
    cutset.cuts[4].tstring.Form("TMath::Abs(recoMuons.eta[1]) < %4.2f", cMaxEta);
    cutset.cuts[4].bins = 50;
    cutset.cuts[4].min = -3;
    cutset.cuts[4].max = 3;

    cutset.cuts[5].name = "recoMuons.iso[1]";
    cutset.cuts[5].tstring.Form("recoMuons.trackIsoSumPt[1]/recoMuons.pt[1] < %4.2f", cMaxRelIso);
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

    cTrigMuPtMin = 24;
                                
    cDimuMassMin = 60;
    cutset.cuts = std::vector<CutInfo>(3, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Run1EventSelectionCuts::Run1EventSelectionCuts(float trigMuPtMin, float dimuMassMin)
{
// Default values for the modified run 1 event selection

    cTrigMuPtMin = trigMuPtMin;      // >
    cDimuMassMin = dimuMassMin;      // >
    cutset.cuts = std::vector<CutInfo>(3, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool Run1EventSelectionCuts::evaluate(VarSet& vars)
{
    // Make sure passed is true at first and the value is reset
    cutset.resetCuts();

    cutset.cuts[0].value = vars.recoMuons.charge[0] != vars.recoMuons.charge[1];
    cutset.cuts[1].value = -1;
    cutset.cuts[2].value = vars.dimuCand.recoCandMassPF;

    // Set cuts[1].value correctly
    if(vars.recoMuons.isHltMatched[0][4]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons.pt[0]);
    else if(vars.recoMuons.isHltMatched[0][5]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons.pt[0]);
    else if(vars.recoMuons.isHltMatched[1][4]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons.pt[1]);
    else if(vars.recoMuons.isHltMatched[1][5]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons.pt[1]);
    else
    {
         cutset.cuts[1].value = -1;
    }

    // if the event fails a single cut return false
    // Require oppositely charged muons in the event
    if(!(vars.recoMuons.charge[0] != vars.recoMuons.charge[1]) && cutset.cuts[0].on) return false;

    // One muon in the pair must pass one of the HLT triggers. This muon have the appropriate pt and eta.
    // Should probably make this into a function so that we can look at a larger number of triggers without cluttering this too much.
    if(!cutset.cuts[1].on) ;
    else if(vars.recoMuons.isHltMatched[0][4] && TMath::Abs(vars.recoMuons.pt[0]) > cTrigMuPtMin) ;
    else if(vars.recoMuons.isHltMatched[0][5] && TMath::Abs(vars.recoMuons.pt[0]) > cTrigMuPtMin) ;
    else if(vars.recoMuons.isHltMatched[1][4] && TMath::Abs(vars.recoMuons.pt[1]) > cTrigMuPtMin) ;
    else if(vars.recoMuons.isHltMatched[1][5] && TMath::Abs(vars.recoMuons.pt[1]) > cTrigMuPtMin) ;
    else
    {
         return false;
    }

    if(!(vars.dimuCand.recoCandMassPF > cDimuMassMin) && cutset.cuts[2].on) return false;

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
    cutset.cuts[0].name = "recoMuons.charge[0] != recoMuons.charge[1]";
    cutset.cuts[0].tstring = "recoMuons.charge[0] != recoMuons.charge[1]";
    cutset.cuts[0].bins = 2;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 2;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = "trigMatchedRecoMu.pt";
    cutset.cuts[1].tstring.Form("recoMuons.isHltMatched[0][4||5] && recoMuons.pt[0] > %5.2f", cTrigMuPtMin);
    cutset.cuts[1].bins = 201;
    cutset.cuts[1].min = -1;
    cutset.cuts[1].max = 200;
    cutset.cuts[1].cutvalue = &cTrigMuPtMin;
    cutset.cuts[1].ismin = true;

    cutset.cuts[2].name = "dimuCand.recoCandMassPF_Min";
    cutset.cuts[2].tstring.Form("dimuCand.recoCandMassPF > %5.2f", cDimuMassMin);
    cutset.cuts[2].bins = 150;
    cutset.cuts[2].min = 50;
    cutset.cuts[2].max = 200;
    cutset.cuts[2].cutvalue = &cDimuMassMin;
    cutset.cuts[2].ismin = true;

}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________Run1EventSelectionSig__________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

Run1EventSelectionSigCuts::Run1EventSelectionSigCuts()
{
// Default values for the modified run 1 event selection

    cTrigMuPtMin = 24;
                                
    cDimuMassMin = 122.5;
    cDimuMassMax = 127.5;
    cutset.cuts = std::vector<CutInfo>(4, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Run1EventSelectionSigCuts::Run1EventSelectionSigCuts(float trigMuPtMin, float dimuMassMin, float dimuMassMax)
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

bool Run1EventSelectionSigCuts::evaluate(VarSet& vars)
{
    // Make sure passed is true at first and the value is reset
    cutset.resetCuts();

    cutset.cuts[0].value = vars.recoMuons.charge[0] != vars.recoMuons.charge[1];
    cutset.cuts[1].value = -1;
    cutset.cuts[2].value = vars.dimuCand.recoCandMassPF;
    cutset.cuts[3].value = vars.dimuCand.recoCandMassPF;

    // Set cuts[1].value correctly
    if(vars.recoMuons.isHltMatched[0][4]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons.pt[0]);
    else if(vars.recoMuons.isHltMatched[0][5]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons.pt[0]);
    else if(vars.recoMuons.isHltMatched[1][4]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons.pt[1]);
    else if(vars.recoMuons.isHltMatched[1][5]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons.pt[1]);
    else
    {
         cutset.cuts[1].value = -1;
    }

    // if the event fails a single cut return false
    // Require oppositely charged muons in the event
    if(!(vars.recoMuons.charge[0] != vars.recoMuons.charge[1]) && cutset.cuts[0].on) return false;

    // One muon in the pair must pass one of the HLT triggers. This muon have the appropriate pt and eta.
    // Should probably make this into a function so that we can look at a larger number of triggers without cluttering this too much.
    if(!cutset.cuts[1].on) ;
    else if(vars.recoMuons.isHltMatched[0][4] && TMath::Abs(vars.recoMuons.pt[0]) > cTrigMuPtMin) ;
    else if(vars.recoMuons.isHltMatched[0][5] && TMath::Abs(vars.recoMuons.pt[0]) > cTrigMuPtMin) ;
    else if(vars.recoMuons.isHltMatched[1][4] && TMath::Abs(vars.recoMuons.pt[1]) > cTrigMuPtMin) ;
    else if(vars.recoMuons.isHltMatched[1][5] && TMath::Abs(vars.recoMuons.pt[1]) > cTrigMuPtMin) ;
    else
    {
         return false;
    }

    if(!(vars.dimuCand.recoCandMassPF > cDimuMassMin) && cutset.cuts[2].on) return false;
    if(!(vars.dimuCand.recoCandMassPF < cDimuMassMax) && cutset.cuts[3].on) return false;

    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString Run1EventSelectionSigCuts::string()
{
    return TString("Run1 Event Selection Cuts");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void Run1EventSelectionSigCuts::makeCutSet()
{
    //cTrigMuPtMin; 
    //cDimuMassMin; 

    cutset.cuts[0].name = "recoMuons.charge[0] != recoMuons.charge[1]";
    cutset.cuts[0].tstring = "recoMuons.charge[0] != recoMuons.charge[1]";
    cutset.cuts[0].bins = 2;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 2;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = "trigMatchedRecoMu.pt";
    cutset.cuts[1].tstring.Form("recoMuons.isHltMatched[0][4||5] && recoMuons.pt[0] > %5.2f", cTrigMuPtMin);
    cutset.cuts[1].bins = 201;
    cutset.cuts[1].min = -1;
    cutset.cuts[1].max = 200;
    cutset.cuts[1].cutvalue = &cTrigMuPtMin;
    cutset.cuts[1].ismin = true;

    cutset.cuts[2].name = "dimuCand.recoCandMassPF_Min";
    cutset.cuts[2].tstring.Form("dimuCand.recoCandMassPF > %5.2f", cDimuMassMin);
    cutset.cuts[2].bins = 150;
    cutset.cuts[2].min = 50;
    cutset.cuts[2].max = 200;
    cutset.cuts[2].cutvalue = &cDimuMassMin;
    cutset.cuts[2].ismin = true;

    cutset.cuts[3].name = "dimuCand.recoCandMassPF_Max";
    cutset.cuts[3].tstring.Form("dimuCand.recoCandMassPF < %5.2f", cDimuMassMax);
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
    //cutset.resetCuts();

    //// keep track of the values that were cut on
    //// recoMuon0
    //cutset.cuts[0].value = TMath::Abs(vars.recoMuons.pt[0]);
    //cutset.cuts[1].value = TMath::Abs(vars.recoMuons.eta[0]);
    //cutset.cuts[2].value = (vars.recoMuons.sumChargedHadronPtR03[0] + TMath::Max(0.0,vars.recoMuons.sumNeutralHadronEtR03[0]+
    //                        vars.recoMuons.sumPhotonEtR03[0] - 0.5*vars.recoMuons.sumPUPtR03[0]))/vars.recoMuons.pt[0];

    //// recoMuon1
    //cutset.cuts[3].value = TMath::Abs(vars.recoMuons.pt[1]);
    //cutset.cuts[4].value = TMath::Abs(vars.recoMuons.eta[1]);
    //cutset.cuts[5].value = (vars.recoMuons.sumChargedHadronPtR03[1] + TMath::Max(0.0,vars.recoMuons.sumNeutralHadronEtR03[1]+
    //                        vars.recoMuons.sumPhotonEtR03[1] - 0.5*vars.recoMuons.sumPUPtR03[1]))/vars.recoMuons.pt[1];

    // if either muon fails the selection return false
    if(!evaluate(vars.recoMuons, 0)) return false;
    if(!evaluate(vars.recoMuons, 1)) return false;
    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool Run1MuonSelectionCuts::evaluate(_MuonInfo& recoMuons, int m)
{
// Test a single muon 

    // if the event fails a single cut return false
    // only return false if the cut is activated for the specific muon m=0or1
    if(!(TMath::Abs(recoMuons.pt[m]) > cMinPt) && ((m==0 && cutset.cuts[0].on) || (m==1 && cutset.cuts[3].on)) ) return false;

    // eta cuts
    if(!(TMath::Abs(recoMuons.eta[m]) < cMaxEta)  && ((m==0 && cutset.cuts[1].on) || (m==1 && cutset.cuts[4].on))) return false;

    // isolation cuts
    if(!((recoMuons.sumChargedHadronPtR03[m] + TMath::Max(0.0,recoMuons.sumNeutralHadronEtR03[m]+recoMuons.sumPhotonEtR03[m] 
          - 0.5*recoMuons.sumPUPtR03[m]))/recoMuons.pt[m] <= cMaxRelIso) && ((m==0 && cutset.cuts[2].on) || (m==1 && cutset.cuts[5].on))) return false;
    
    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString Run1MuonSelectionCuts::string()
{
    return TString("Run1_Muon_Selection_Cuts");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void Run1MuonSelectionCuts::makeCutSet()
{
    //cMinPt 
    //cMaxEta
    //cMaxRelIso

    // recoMuon0 cuts
    cutset.cuts[0].name = Form("recoMu0.pt > %4.1f", cMinPt);
    cutset.cuts[0].tstring.Form("recoMuons.pt[0] > %4.1f", cMinPt);
    cutset.cuts[0].bins = 200;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 200;
    cutset.cuts[0].cutvalue = &cMinPt;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = Form("|recoMu0.eta| < %4.2f", cMaxEta);
    cutset.cuts[1].tstring.Form("abs(recoMuons.eta[0]) < %4.2f", cMaxEta);
    cutset.cuts[1].bins = 50;
    cutset.cuts[1].min = 0;
    cutset.cuts[1].max = 3;
    cutset.cuts[1].cutvalue = &cMaxEta;
    cutset.cuts[1].ismin = false;

    cutset.cuts[2].name = Form("recoMu0.iso <= %4.2f", cMaxRelIso);
    cutset.cuts[2].tstring.Form("((recoMuons.sumChargedHadronPtR03[0] + max(0.0,recoMuons.sumNeutralHadronEtR03[0]+recoMuons.sumPhotonEtR03[0]"+ 
                                 TString("- 0.5*recoMuons.sumPUPtR03[0]))/recoMuons.pt[0] <= %4.2f)"), cMaxRelIso);
    cutset.cuts[2].bins = 50;
    cutset.cuts[2].min = 0;
    cutset.cuts[2].max = 1;
    cutset.cuts[2].cutvalue = &cMaxRelIso;
    cutset.cuts[2].ismin = false;

     // recoMuon1 cuts
    cutset.cuts[3].name = Form("recoMu1.pt > %4.1f", cMinPt);
    cutset.cuts[3].tstring.Form("recoMuons.pt[1] > %f", cMinPt);
    cutset.cuts[3].bins = 200;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 200;
    cutset.cuts[3].cutvalue = &cMinPt;
    cutset.cuts[3].ismin = true;

    cutset.cuts[4].name = Form("|recoMu1.eta| < %4.2f", cMaxEta);
    cutset.cuts[4].tstring.Form("abs(recoMuons.eta[1]) < %4.2f", cMaxEta);
    cutset.cuts[4].bins = 50;
    cutset.cuts[4].min = 0;
    cutset.cuts[4].max = 3;
    cutset.cuts[4].cutvalue = &cMaxEta;
    cutset.cuts[4].ismin = false;

    cutset.cuts[5].name = Form("recoMu1.iso <= %4.2f", cMaxRelIso);
    cutset.cuts[5].tstring.Form("((recoMuons.sumChargedHadronPtR03[1] + max(0.0,recoMuons.sumNeutralHadronEtR03[1]+recoMuons.sumPhotonEtR03[1]"+ 
                                TString("- 0.5*recoMuons.sumPUPtR03[1]))/recoMuons.pt[1] <= %4.2f)"), cMaxRelIso);
    cutset.cuts[5].bins = 50;
    cutset.cuts[5].min = 0;
    cutset.cuts[5].max = 1;
    cutset.cuts[5].cutvalue = &cMaxRelIso;
    cutset.cuts[5].ismin = false;

    cutset.concatCuts();
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________Run1EventSelectionCuts80X______________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

Run1EventSelectionCuts80X::Run1EventSelectionCuts80X()
{
// Default values for the modified run 1 event selection

    isData = 0;
    cTrigMuPtMin = 24; 
                                
    cDimuMassMin = 60;
    cutset.cuts = std::vector<CutInfo>(4, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Run1EventSelectionCuts80X::Run1EventSelectionCuts80X(bool isData)
{
// Default values for the modified run 1 event selection

    this->isData = isData;;
    cTrigMuPtMin = 24; 
                                
    cDimuMassMin = 60;
    cutset.cuts = std::vector<CutInfo>(4, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Run1EventSelectionCuts80X::Run1EventSelectionCuts80X(float trigMuPtMin, float dimuMassMin)
{
// Default values for the modified run 1 event selection

    isData = 0;
    cTrigMuPtMin = trigMuPtMin;      // >
    cDimuMassMin = dimuMassMin;      // >
    cutset.cuts = std::vector<CutInfo>(4, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Run1EventSelectionCuts80X::Run1EventSelectionCuts80X(bool isData, float trigMuPtMin, float dimuMassMin)
{
// Default values for the modified run 1 event selection

    this->isData = isData;
    cTrigMuPtMin = trigMuPtMin;      // >
    cDimuMassMin = dimuMassMin;      // >
    cutset.cuts = std::vector<CutInfo>(4, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool Run1EventSelectionCuts80X::evaluate(VarSet& vars)
{
    // Make sure passed is true at first and the value is reset
    //cutset.resetCuts();

    //cutset.cuts[0].value = vars.recoMuons.charge[0] != vars.recoMuons.charge[1];
    //cutset.cuts[1].value = -1;
    //cutset.cuts[2].value = vars.dimuCand.recoCandMassPF;

    // Set cuts[1].value correctly
    //if(vars.recoMuons.pt[0] >= vars.recoMuons.pt[1]) cutset.cuts[1].value = vars.recoMuons.pt[0];
    //else cutset.cuts[1].value = vars.recoMuons.pt[1];

    if(!(vars.dimuCand.recoCandMassPF > cDimuMassMin) && cutset.cuts[2].on) return false;

    // Pt cuts on leading and subleading muons
    if(!cutset.cuts[1].on) ;
    else if(vars.recoMuons.pt[0] >= vars.recoMuons.pt[1] && TMath::Abs(vars.recoMuons.pt[0]) > cTrigMuPtMin) ;
    else if(vars.recoMuons.pt[0] < vars.recoMuons.pt[1] && TMath::Abs(vars.recoMuons.pt[1]) > cTrigMuPtMin) ;
    else
    {
         return false;
    }

    // trigger matching
    if(isData)
    {
        if(!(vars.recoMuons.isHltMatched[0][4] || vars.recoMuons.isHltMatched[0][5] || 
             vars.recoMuons.isHltMatched[1][4] || vars.recoMuons.isHltMatched[1][5]))
            return false;
    }

    // Require oppositely charged muons in the event
    if(!(vars.recoMuons.charge[0] != vars.recoMuons.charge[1]) && cutset.cuts[0].on) return false;

    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString Run1EventSelectionCuts80X::string()
{
    return TString("Run1_Event_Selection_80X");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void Run1EventSelectionCuts80X::makeCutSet()
{
    cutset.cuts[0].name = "recoMu0.charge != recoMu1.charge";
    cutset.cuts[0].tstring = "recoMuons.charge[0] != recoMuons.charge[1]";
    cutset.cuts[0].bins = 2;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 2;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = Form("leadingRecoMu.pt > %6.2f", cTrigMuPtMin);
    cutset.cuts[1].tstring.Form("max(recoMuons.pt[0], recoMuons.pt[1]) > %6.2f", cTrigMuPtMin);
    cutset.cuts[1].bins = 201;
    cutset.cuts[1].min = -1;
    cutset.cuts[1].max = 200;
    cutset.cuts[1].cutvalue = &cTrigMuPtMin;
    cutset.cuts[1].ismin = true;

    cutset.cuts[2].name = Form("recoCandMassPF > %6.2f", cDimuMassMin);
    cutset.cuts[2].tstring.Form("dimuCand.recoCandMassPF > %6.2f", cDimuMassMin);
    cutset.cuts[2].bins = 150;
    cutset.cuts[2].min = 50;
    cutset.cuts[2].max = 200;
    cutset.cuts[2].cutvalue = &cDimuMassMin;
    cutset.cuts[2].ismin = true;

    cutset.cuts[3].name = "recoMu0 or recoMu1 matched to (hlt4 or hlt5)";
    cutset.cuts[3].tstring = "(recoMuons.isHltMatched[0][4] || recoMuons.isHltMatched[0][5] || recoMuons.isHltMatched[1][4] || recoMuons.isHltMatched[1][5])";
    cutset.cuts[3].bins = 2;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 2;
    cutset.cuts[3].ismin = true;

    cutset.concatCuts();
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________FEWZCompareCuts________________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

FEWZCompareCuts::FEWZCompareCuts()
{
// Default values for the fewz comparison selection cuts

    useReco = false;
    cLeadPtMin = 20;   
    cSubleadPtMin = 10;   
    cMaxEta = 2.4;     
    cDimuMassMin = 110; 
    cDimuMassMax = 310; 
    cMaxRelIso = 0.12;
    cutset.cuts = std::vector<CutInfo>(9, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

FEWZCompareCuts::FEWZCompareCuts(bool useReco)
{
// Default values for the fewz comparison selection cuts
    this->useReco = useReco;
    cLeadPtMin = 20;   
    cSubleadPtMin = 10;   
    cMaxEta = 2.4;     
    cDimuMassMin = 110; 
    cDimuMassMax = 310; 
    cMaxRelIso = 0.12;      
    cutset.cuts = std::vector<CutInfo>(9, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

FEWZCompareCuts::FEWZCompareCuts(bool useReco, float leadPtMin, float subleadPtMin, float maxEta, float dimuMassMin, float dimuMassMax, float maxRelIso)
{
// Custom selection cuts for the fewz comparison
    this->useReco = useReco;
    cLeadPtMin = leadPtMin;
    cSubleadPtMin = subleadPtMin;
    cMaxEta = maxEta;
    cDimuMassMin = dimuMassMin;
    cDimuMassMax = dimuMassMax;
    cMaxRelIso = maxRelIso; 
    cutset.cuts = std::vector<CutInfo>(9, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool FEWZCompareCuts::evaluate(VarSet& vars)
{
    // Make sure passed is true at first and the value is reset
    cutset.resetCuts();

    float leadPt = -999;
    float subleadPt = -999; 
    float eta0 = -999;
    float eta1 = -999;
    float dimu_mass = -999;
    int charge0 = -999;
    int charge1 = -999;

    if(useReco)
    {
        leadPt = TMath::Max(TMath::Abs(vars.recoMuons.pt[0]), TMath::Abs(vars.recoMuons.pt[1]));
        subleadPt = TMath::Min(TMath::Abs(vars.recoMuons.pt[0]), TMath::Abs(vars.recoMuons.pt[1]));
        eta0 = vars.recoMuons.eta[0];
        eta1 = vars.recoMuons.eta[1];
        dimu_mass = vars.dimuCand.recoCandMass;
        charge0 = vars.recoMuons.charge[0];
        charge1 = vars.recoMuons.charge[1];
    }
    // use gen values
    else
    {
        // get first postFSR DY gen muon
        _TrackInfo mu0 = ParticleTools::getGenMuDY(0, 1, vars);
        // get second postFSR DY gen muon
        _TrackInfo mu1 = ParticleTools::getGenMuDY(1, 1, vars);

        if(mu0.charge == mu1.charge)
        {
            std::cout << "### ERROR: DYmu0.charge == DYmu1.charge in FEWZ SELECTION. Muons decayed from Z or gamma* should have opposite sign." << std::endl;
            std::cout << "### And this shouldn't happen. This event will be cut due to charge selection anyways." << std::endl;
        }

        eta0 = mu0.eta;
        eta1 = mu1.eta;
        charge0 = mu0.charge;
        charge1 = mu1.charge;

        leadPt = TMath::Max(TMath::Abs(mu0.pt), TMath::Abs(mu1.pt));
        subleadPt = TMath::Min(TMath::Abs(mu0.pt), TMath::Abs(mu1.pt));

        if(mu0.pt > 0 && mu1.pt > 0) dimu_mass = ParticleTools::getMotherPtEtaPhiM(mu0.pt, mu0.eta, mu0.phi, MASS_MUON, mu1.pt, mu1.eta, mu1.phi, MASS_MUON).M();
        else
        {
            std::cout << "### ERROR: gen_dimu_mass < 0 in FEWZ SELECTION" << std::endl;
            std::cout << vars.eventInfo.run << "," << vars.eventInfo.lumi << "," << vars.eventInfo.event << std::endl;
            std::cout << mu0.pt << "," << mu0.eta << "," << mu0.phi << std::endl;
            std::cout << mu1.pt << "," << mu1.eta << "," << mu1.phi << std::endl;
            std::cout << std::endl;
        }
    }

    // keep track of the values that were cut on
    cutset.cuts[0].value = TMath::Abs(leadPt);
    cutset.cuts[1].value = TMath::Abs(subleadPt);
    cutset.cuts[2].value = TMath::Abs(eta0);
    cutset.cuts[3].value = TMath::Abs(eta1);
    cutset.cuts[4].value = TMath::Abs(dimu_mass);
    cutset.cuts[5].value = TMath::Abs(dimu_mass);
    cutset.cuts[6].value = charge0 != charge1;

    // if a muon fails any of the criterea that are turned on then reutrn false
    if(cutset.cuts[0].on)
    {
        if(leadPt <= cLeadPtMin) return false;
    }
    if(cutset.cuts[1].on)
    {
        if(subleadPt <= cSubleadPtMin) return false;
    }
    if(cutset.cuts[2].on)
    {
        if(TMath::Abs(eta0) >= cMaxEta) return false;
    }
    if(cutset.cuts[3].on)
    {
        if(TMath::Abs(eta1) >= cMaxEta) return false;
    }
    if(cutset.cuts[4].on)
    {
        if(dimu_mass <= cDimuMassMin) return false;
    }
    if(cutset.cuts[5].on)
    {
        if(dimu_mass >= cDimuMassMax) return false;
    }
    if(cutset.cuts[6].on)
    {
        if(charge0 == charge1) return false;
    }
    if(cutset.cuts[7].on)
    {
        if(!((vars.recoMuons.sumChargedHadronPtR03[0] + TMath::Max(0.0,vars.recoMuons.sumNeutralHadronEtR03[0]+vars.recoMuons.sumPhotonEtR03[0] 
              - 0.5*vars.recoMuons.sumPUPtR03[0]))/vars.recoMuons.pt[0] <= cMaxRelIso) ) return false;
    }
    if(cutset.cuts[8].on)
    {
        if(!((vars.recoMuons.sumChargedHadronPtR03[1] + TMath::Max(0.0,vars.recoMuons.sumNeutralHadronEtR03[1]+vars.recoMuons.sumPhotonEtR03[1] 
              - 0.5*vars.recoMuons.sumPUPtR03[1]))/vars.recoMuons.pt[1] <= cMaxRelIso) ) return false;
    }
    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString FEWZCompareCuts::string()
{
    return TString("FEWZ Comparison Cuts");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void FEWZCompareCuts::makeCutSet()
{
    //cLeadPtMin 
    //cSubleadPtMin
    //cMaxEta1
    //cMaxEta2
    //cDimuMassMin
    //cDimuMassMax
    //charge
    //iso0
    //iso1

    // recoMuon0 cuts
    cutset.cuts[0].name = "recoLead.pt";
    cutset.cuts[0].tstring.Form("recoLead.pt > %4.1f", cLeadPtMin);
    cutset.cuts[0].bins = 200;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 200;
    cutset.cuts[0].cutvalue = &cLeadPtMin;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = "recoSublead.pt";
    cutset.cuts[1].tstring.Form("recoSublead.pt > %4.1f", cSubleadPtMin);
    cutset.cuts[1].bins = 200;
    cutset.cuts[1].min = 0;
    cutset.cuts[1].max = 200;
    cutset.cuts[1].cutvalue = &cSubleadPtMin;
    cutset.cuts[1].ismin = true;

    cutset.cuts[2].name = "recoMuons.eta[0]";
    cutset.cuts[2].tstring.Form("TMath::Abs(recoMuons.eta[0]) < %4.2f", cMaxEta);
    cutset.cuts[2].bins = 50;
    cutset.cuts[2].min = 0;
    cutset.cuts[2].max = 3;
    cutset.cuts[2].cutvalue = &cMaxEta;
    cutset.cuts[2].ismin = false;

    cutset.cuts[3].name = "recoMuons.eta[1]";
    cutset.cuts[3].tstring.Form("TMath::Abs(recoMuons.eta[1]) < %4.2f", cMaxEta);
    cutset.cuts[3].bins = 50;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 3;
    cutset.cuts[3].cutvalue = &cMaxEta;
    cutset.cuts[3].ismin = false;

    cutset.cuts[4].name = "dimuCand.recoCandMass_Min";
    cutset.cuts[4].tstring.Form("dimuCand.recoCandMass < %5.2f", cDimuMassMin);
    cutset.cuts[4].bins = 100;
    cutset.cuts[4].min = 110;
    cutset.cuts[4].max = 310;
    cutset.cuts[4].cutvalue = &cDimuMassMin;
    cutset.cuts[4].ismin = true;

    cutset.cuts[5].name = "dimuCand.recoCandMass_Max";
    cutset.cuts[5].tstring.Form("dimuCand.recoCandMass < %5.2f", cDimuMassMax);
    cutset.cuts[5].bins = 100;
    cutset.cuts[5].min = 110;
    cutset.cuts[5].max = 310;
    cutset.cuts[5].cutvalue = &cDimuMassMax;
    cutset.cuts[5].ismin = false;

    cutset.cuts[6].name = "recoMuons.charge[0] != recoMuons.charge[1]";
    cutset.cuts[6].tstring = "recoMuons.charge[0] != recoMuons.charge[1]";
    cutset.cuts[6].bins = 2;
    cutset.cuts[6].min = 0;
    cutset.cuts[6].max = 2;
    cutset.cuts[6].ismin = true;

    cutset.cuts[7].name = "recoMuons.iso[0]";
    cutset.cuts[7].tstring.Form("recoMuons.iso[0] < %4.2f", cMaxRelIso);
    cutset.cuts[7].bins = 50;
    cutset.cuts[7].min = 0;
    cutset.cuts[7].max = 1;
    cutset.cuts[7].cutvalue = &cMaxRelIso;
    cutset.cuts[7].ismin = false;

    cutset.cuts[8].name = "recoMuons.iso[1]";
    cutset.cuts[8].tstring.Form("recoMuons.iso[1] < %4.2f", cMaxRelIso);
    cutset.cuts[8].bins = 50;
    cutset.cuts[8].min = 0;
    cutset.cuts[8].max = 1;
    cutset.cuts[8].cutvalue = &cMaxRelIso;
    cutset.cuts[8].ismin = false;
}

