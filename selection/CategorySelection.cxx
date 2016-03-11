///////////////////////////////////////////////////////////////////////////
//                             CategorySelection.cxx                     //
//=======================================================================//
//                                                                       //
//        Define the selection critera for the different categories.     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "CategorySelection.h"

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________CategorySelection______________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

CategorySelection::CategorySelection()
{
// Standard values for the tight muon id cuts

    // Preselection
    cLeadPtMin = 40;
    cSubleadPtMin = 30;
    cMETMax = 40;
    isPreselected = false;

    // VBF Tight
    cDijetMassMinVBFT = 650;
    cDijetDeltaEtaMinVBFT = 3.5;
    isVBFTight = false;

    // VBF Loose
    isVBFLoose = false;

    // GGF Tight
    cDijetMassMinGGFT = 250;
    cDimuPtMinGGFT = 50;
    isGGFTight = false;

    // 01Tight
    cDimuPtMin01T = 10;
    isTight01 = false;

    // 01Loose
    isLoose01 = false;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelection::CategorySelection(float leadPtMin, float subleadPtMin, float METMax, float dijetMassMinVBFT, float dijetDeltaEtaMinVBFT, float dijetMassMinGGFT,
                                     float dimuPtMinGGFT, float dimuPtMin01T)
{
// Standard values for the tight muon id cuts

    // Preselection
    cLeadPtMin = leadPtMin;
    cSubleadPtMin = subleadPtMin;
    cMETMax = METMax;
    isPreselected = false;

    // VBF Tight
    cDijetMassMinVBFT = dijetMassMinVBFT;
    cDijetDeltaEtaMinVBFT = dijetDeltaEtaMinVBFT;
    isVBFTight = false;

    // VBF Loose
    isVBFLoose = false;

    // GGF Tight
    cDijetMassMinGGFT = dijetMassMinGGFT;
    cDimuPtMinGGFT = dimuPtMinGGFT;
    isGGFTight = false;

    // 01Tight
    cDimuPtMin01T = dimuPtMin01T;
    isTight01 = false;

    // 01Loose
    isLoose01 = false;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void CategorySelection::evaluate(VarSet& vars)
{
// Determine which category the event belongs to

    // preselection
    if(vars.validJets.size() == 2)
    {
        TLorentzVector leadJet = vars.validJets[0];
        TLorentzVector subleadJet = vars.validJets[1];
        TLorentzVector dijet = leadJet + subleadJet;

        float dEta = leadJet.Eta() - subleadJet.Eta();
        float dijetMass = dijet.M();

        if(leadJet.Pt() > cLeadPtMin && subleadJet.Pt() > cSubleadPtMin && vars.met.pt < cMETMax)
        {
            isPreselected = true;
            if(dijetMass > cDijetMassMinVBFT && TMath::Abs(dEta) > cDijetDeltaEtaMinVBFT){ isVBFTight = true; return; }
            else if(dijetMass > cDijetMassMinGGFT && vars.recoCandPt > cDimuPtMinGGFT){ isGGFTight = true; return; }
            else{ isVBFLoose = true; return; }
        }
    }
    if(!isPreselected) // fails 2jet preselection enters 01 categories
    {
        if(vars.recoCandPt > cDimuPtMin01T){ isTight01 = true; return; }
        else{ isLoose01 = true; return; }
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void CategorySelection::reset()
{
// Reset the boolean values for the next iteration

    isPreselected = false;
    isVBFTight = false;
    isVBFLoose = false;
    isGGFTight = false;
    isTight01 = false;
    isLoose01 = false;
}
