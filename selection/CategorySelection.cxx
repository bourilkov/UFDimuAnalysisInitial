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

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________CategorySelectionFEWZ__________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

CategorySelectionFEWZ::CategorySelectionFEWZ()
{
// Standard values for the FEWZ categorization

    // init cut values
    cMassSplit = 160;
    cEtaCentralSplit = 0.8;
    cJetPtMin = 30;
    cJetEtaMax = 4.7;

    // Intitial Tests
    isWide = false;
    isNarrow = false;
    isCentralCentral = false;
    isCentralNotCentral = false;
    isOneJetInclusive = false;

    // Final Categories
    isCentralCentralWide = false;
    isCentralCentralNarrow = false;

    isCentralNotCentralWide = false;
    isCentralNotCentralNarrow = false;

    isOneJetInclusiveWide = false;
    isOneJetInclusiveNarrow = false;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelectionFEWZ::CategorySelectionFEWZ(float massSplit, float etaCentralSplit, float jetPtMin, float jetEtaMax)
{
// Custom values for the cuts

    // init cut values
    cMassSplit = massSplit;
    cEtaCentralSplit = etaCentralSplit;
    cJetPtMin = jetPtMin;
    cJetEtaMax = jetEtaMax;

    // Intitial Tests
    isWide = false;
    isNarrow = false;
    isCentralCentral = false;
    isCentralNotCentral = false;
    isOneJetInclusive = false;

    // Final Categories
    isCentralCentralWide = false;
    isCentralCentralNarrow = false;

    isCentralNotCentralWide = false;
    isCentralNotCentralNarrow = false;

    isOneJetInclusiveWide = false;
    isOneJetInclusiveNarrow = false;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void CategorySelectionFEWZ::evaluate(VarSet& vars)
{
// Determine which category the event belongs to

    // Should cut out all events that don't fall into the wide mass window in earlier selection stage
    // All events that pass are in window of min to max
    isWide = true;

    // Narrow goes from min to cMassSplit
    if(vars.recoCandMass < cMassSplit) isNarrow = true;

    // Both central
    if(TMath::Abs(vars.reco1.eta) < 0.8 && TMath::Abs(vars.reco2.eta) < 0.8) isCentralCentral = true;

    // Not both, but at least one is central
    else if(TMath::Abs(vars.reco1.eta) < 0.8 || TMath::Abs(vars.reco2.eta) < 0.8) isCentralNotCentral = true;

    // One category that passes basic selections and has exactly one jet
    if(vars.validJets.size() == 1) isOneJetInclusive = true; 

    // Final Categories ///////////////////////////////////////////////////////
    if(isWide && isCentralCentral) isCentralCentralWide = true;
    if(isNarrow && isCentralCentral) isCentralCentralNarrow = true;

    if(isWide && isCentralNotCentral) isCentralNotCentralWide = true;
    if(isNarrow && isCentralNotCentral) isCentralNotCentralNarrow = true;

    if(isWide && isOneJetInclusive) isOneJetInclusiveWide = true;
    if(isNarrow && isOneJetInclusive) isOneJetInclusiveNarrow = true;

    return;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void CategorySelectionFEWZ::reset()
{
// Reset the boolean values for the next iteration

    // Intitial Tests
    isWide = false;
    isNarrow = false;
    isCentralCentral = false;
    isCentralNotCentral = false;
    isOneJetInclusive = false;

    // Final Categories
    isCentralCentralWide = false;
    isCentralCentralNarrow = false;

    isCentralNotCentralWide = false;
    isCentralNotCentralNarrow = false;

    isOneJetInclusiveWide = false;
    isOneJetInclusiveNarrow = false;
}
