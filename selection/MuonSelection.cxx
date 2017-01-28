///////////////////////////////////////////////////////////////////////////
//                             MuonSelection.cxx                         //
//=======================================================================//
//                                                                       //
//        Define the muon selection cuts for the analysis.               //
//        The dimuon candidate muons must pass these criteria.           //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "MuonSelection.h"
#include "ParticleTools.h"
#include "TMath.h"
#include <iostream>

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________Run1MuonSelection______________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

// muon selection based on run 1, still used for run 2

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
    //cutset.cuts[0].value = TMath::Abs(vars.recoMuons->at(0).pt);
    //cutset.cuts[1].value = TMath::Abs(vars.recoMuons->at(0).eta);
    //cutset.cuts[2].value = (vars.recoMuons->at(0).sumChargedHadronPtR03 + TMath::Max(0.0,vars.recoMuons->at(0).sumNeutralHadronEtR03+
    //                        vars.recoMuons->at(0).sumPhotonEtR03 - 0.5*vars.recoMuons->at(0).sumPUPtR03))/vars.recoMuons->at(0).pt;

    //// recoMuon1
    //cutset.cuts[3].value = TMath::Abs(vars.recoMuons->at(1).pt);
    //cutset.cuts[4].value = TMath::Abs(vars.recoMuons->at(1).eta);
    //cutset.cuts[5].value = (vars.recoMuons->at(1).sumChargedHadronPtR03 + TMath::Max(0.0,vars.recoMuons->at(1).sumNeutralHadronEtR03+
    //                        vars.recoMuons->at(1).sumPhotonEtR03 - 0.5*vars.recoMuons->at(1).sumPUPtR03))/vars.recoMuons->at(1).pt;

    // if either muon fails the selection return false
    if(!evaluate(vars.recoMuons, 0)) return false;
    if(!evaluate(vars.recoMuons, 1)) return false;
    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool Run1MuonSelectionCuts::evaluate(std::vector<MuonInfo>* recoMuons, int m)
{
// Test a single muon 

    // if the event fails a single cut return false
    // only return false if the cut is activated for the specific muon m=0or1
    if(!(TMath::Abs(recoMuons->at(m).pt) > cMinPt) && ((m==0 && cutset.cuts[0].on) || (m==1 && cutset.cuts[3].on)) ) return false;

    // eta cuts
    if(!(TMath::Abs(recoMuons->at(m).eta) < cMaxEta)  && ((m==0 && cutset.cuts[1].on) || (m==1 && cutset.cuts[4].on))) return false;

    // isolation cuts
    if(!((recoMuons->at(m).sumChargedHadronPtR03 + TMath::Max(0.0,recoMuons->at(m).sumNeutralHadronEtR03+recoMuons->at(m).sumPhotonEtR03 
          - 0.5*recoMuons->at(m).sumPUPtR03))/recoMuons->at(m).pt <= cMaxRelIso) && 
          ((m==0 && cutset.cuts[2].on) || (m==1 && cutset.cuts[5].on))) return false;
    
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
    cutset.cuts[0].tstring.Form("recoMuons->at(0).pt > %4.1f", cMinPt);
    cutset.cuts[0].bins = 200;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 200;
    cutset.cuts[0].cutvalue = &cMinPt;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = Form("|recoMu0.eta| < %4.2f", cMaxEta);
    cutset.cuts[1].tstring.Form("abs(recoMuons->at(0).eta) < %4.2f", cMaxEta);
    cutset.cuts[1].bins = 50;
    cutset.cuts[1].min = 0;
    cutset.cuts[1].max = 3;
    cutset.cuts[1].cutvalue = &cMaxEta;
    cutset.cuts[1].ismin = false;

    cutset.cuts[2].name = Form("recoMu0.iso <= %4.2f", cMaxRelIso);
    cutset.cuts[2].tstring.Form("((recoMuons->at(0).sumChargedHadronPtR03 + max(0.0,recoMuons->at(0).sumNeutralHadronEtR03+"+
                                 TString("recoMuons->at(0).sumPhotonEtR03")+ 
                                 TString("- 0.5*recoMuons->at(0).sumPUPtR03))/recoMuons->at(0).pt <= %4.2f)"), cMaxRelIso);
    cutset.cuts[2].bins = 50;
    cutset.cuts[2].min = 0;
    cutset.cuts[2].max = 1;
    cutset.cuts[2].cutvalue = &cMaxRelIso;
    cutset.cuts[2].ismin = false;

     // recoMuon1 cuts
    cutset.cuts[3].name = Form("recoMu1.pt > %4.1f", cMinPt);
    cutset.cuts[3].tstring.Form("recoMuons->at(1).pt > %f", cMinPt);
    cutset.cuts[3].bins = 200;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 200;
    cutset.cuts[3].cutvalue = &cMinPt;
    cutset.cuts[3].ismin = true;

    cutset.cuts[4].name = Form("|recoMu1.eta| < %4.2f", cMaxEta);
    cutset.cuts[4].tstring.Form("abs(recoMuons->at(1).eta) < %4.2f", cMaxEta);
    cutset.cuts[4].bins = 50;
    cutset.cuts[4].min = 0;
    cutset.cuts[4].max = 3;
    cutset.cuts[4].cutvalue = &cMaxEta;
    cutset.cuts[4].ismin = false;

    cutset.cuts[5].name = Form("recoMu1.iso <= %4.2f", cMaxRelIso);
    cutset.cuts[5].tstring.Form("((recoMuons->at(1).sumChargedHadronPtR03 + max(0.0,recoMuons->at(1).sumNeutralHadronEtR03+"+
                                TString("recoMuons->at(1).sumPhotonEtR03")+ 
                                TString("- 0.5*recoMuons->at(1).sumPUPtR03))/recoMuons->at(1).pt <= %4.2f)"), cMaxRelIso);
    cutset.cuts[5].bins = 50;
    cutset.cuts[5].min = 0;
    cutset.cuts[5].max = 1;
    cutset.cuts[5].cutvalue = &cMaxRelIso;
    cutset.cuts[5].ismin = false;

    cutset.concatCuts();
}
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________SynchMuonSelection_____________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

// Muon selection cuts used for synchronization exercise.


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
    cutset.cuts[0].value = vars.recoMuons->at(0).pt;
    cutset.cuts[1].value = vars.recoMuons->at(0).eta;
    cutset.cuts[2].value = vars.recoMuons->at(0).trackIsoSumPt/vars.recoMuons->at(0).pt;

    // recoMuon1
    cutset.cuts[3].value = vars.recoMuons->at(1).pt;
    cutset.cuts[4].value = vars.recoMuons->at(1).eta;
    cutset.cuts[5].value = vars.recoMuons->at(1).trackIsoSumPt/vars.recoMuons->at(1).pt;

    // if either muon fails the selection return false
    if(!(evaluate(vars.recoMuons, 0))) return false;
    if(!(evaluate(vars.recoMuons, 1))) return false;
    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool SynchMuonSelectionCuts::evaluate(std::vector<MuonInfo>* recoMuons, int m)
{
    // if the event fails a single cut return false
    // only return false if the cut is activated for the specific muon m=1or2
    if(!(recoMuons->at(m).pt > cMinPt) && ((m==0 && cutset.cuts[0].on) || (m==1 && cutset.cuts[3].on)) ) return false;

    // eta cuts
    if(!(TMath::Abs(recoMuons->at(m).eta) < cMaxEta)  && ((m==0 && cutset.cuts[1].on) || (m==1 && cutset.cuts[4].on))) return false;

    // isolation cuts
    if(!(recoMuons->at(m).trackIsoSumPt/recoMuons->at(m).pt < cMaxRelIso)  && ((m==0 && cutset.cuts[2].on) || (m==1 && cutset.cuts[5].on))) return false;
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
    cutset.cuts[0].name = "recoMuons->at(0).pt";
    cutset.cuts[0].tstring.Form("recoMuons->at(0).pt > %4.1f", cMinPt);
    cutset.cuts[0].bins = 200;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 200;

    cutset.cuts[1].name = "recoMuons->at(0).eta";
    cutset.cuts[1].tstring.Form("TMath::Abs(recoMuons->at(0).eta) < %4.2f", cMaxEta);
    cutset.cuts[1].bins = 50;
    cutset.cuts[1].min = -3;
    cutset.cuts[1].max = 3;

    cutset.cuts[2].name = "recoMuons.iso[0]";
    cutset.cuts[2].tstring.Form("recoMuons->at(0).trackIsoSumPt/recoMuons->at(0).pt < %4.2f", cMaxRelIso);
    cutset.cuts[2].bins = 100;
    cutset.cuts[2].min = 0;
    cutset.cuts[2].max = 5;

    cutset.cuts[3].name = "recoMuons->at(1).pt";
    cutset.cuts[3].tstring.Form("recoMuons->at(1).pt > %f", cMinPt);
    cutset.cuts[3].bins = 200;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 200;

    cutset.cuts[4].name = "recoMuons->at(1).eta";
    cutset.cuts[4].tstring.Form("TMath::Abs(recoMuons->at(1).eta) < %4.2f", cMaxEta);
    cutset.cuts[4].bins = 50;
    cutset.cuts[4].min = -3;
    cutset.cuts[4].max = 3;

    cutset.cuts[5].name = "recoMuons.iso[1]";
    cutset.cuts[5].tstring.Form("recoMuons->at(1).trackIsoSumPt/recoMuons->at(1).pt < %4.2f", cMaxRelIso);
    cutset.cuts[5].bins = 100;
    cutset.cuts[5].min = 0;
    cutset.cuts[5].max = 5;
}
