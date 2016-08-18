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

#include "CategorySelection_v2.h"

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________CategorySelection______________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void CategorySelection::initCategoryMap()
{
// Initialize the categories
    categoryMap["ALL"] = Category("ALL");
    categoryMap["Preselection"] = Category("Preselection");
    categoryMap["VBF_Tight"] = Category("VBF_Tight");
    categoryMap["VBF_Loose"] = Category("VBF_Loose");
    categoryMap["GGF_Tight"] = Category("GGF_Tight");
    categoryMap["01_Jet_Tight"] = Category("01_Jet_Tight");
    categoryMap["01_Jet_Loose"] = Category("01_Jet_Loose");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelection::CategorySelection()
{
// initialize the default cut values in the constructor

    initCategoryMap();

    // Preselection
    cLeadPtMin = 40;
    cSubleadPtMin = 30;
    cMETMax = 40;

    // VBF Tight
    cDijetMassMinVBFT = 650;
    cDijetDeltaEtaMinVBFT = 3.5;

    // VBF Loose

    // GGF Tight
    cDijetMassMinGGFT = 250;
    cDimuPtMinGGFT = 50;

    // 01Tight
    cDimuPtMin01T = 10;

    // 01Loose
    
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelection::CategorySelection(float leadPtMin, float subleadPtMin, float METMax, float dijetMassMinVBFT, float dijetDeltaEtaMinVBFT, float dijetMassMinGGFT,
                                     float dimuPtMinGGFT, float dimuPtMin01T)
{
// Initialize custom cut values

    initCategoryMap();

    // Preselection
    cLeadPtMin = leadPtMin;
    cSubleadPtMin = subleadPtMin;
    cMETMax = METMax;

    // VBF Tight
    cDijetMassMinVBFT = dijetMassMinVBFT;
    cDijetDeltaEtaMinVBFT = dijetDeltaEtaMinVBFT;

    // VBF Loose

    // GGF Tight
    cDijetMassMinGGFT = dijetMassMinGGFT;
    cDimuPtMinGGFT = dimuPtMinGGFT;

    // 01Tight
    cDimuPtMin01T = dimuPtMin01T;

    // 01Loose
}


///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void CategorySelection::evaluate(VarSet& vars)
{
// Determine which category the event belongs to

    // Inclusive category, all events that passed the selection cuts
    categoryMap["ALL"].inCategory = true;

    // preselection
    if(vars.validJets.size() >= 2)
    {
        TLorentzVector leadJet = vars.validJets[0];
        TLorentzVector subleadJet = vars.validJets[1];
        TLorentzVector dijet = leadJet + subleadJet;

        float dEta = leadJet.Eta() - subleadJet.Eta();
        float dijetMass = dijet.M();

        if(leadJet.Pt() > cLeadPtMin && subleadJet.Pt() > cSubleadPtMin && vars.met.pt < cMETMax)
        {
            categoryMap["Preselection"].inCategory = true;
            if(dijetMass > cDijetMassMinVBFT && TMath::Abs(dEta) > cDijetDeltaEtaMinVBFT){ categoryMap["VBF_Tight"].inCategory = true; return; }
            else if(dijetMass > cDijetMassMinGGFT && vars.recoCandPtPF > cDimuPtMinGGFT){ categoryMap["GGF_Tight"].inCategory = true; return; }
            else{ categoryMap["VBF_Loose"].inCategory = true; return; }
        }
    }
    if(!categoryMap["Preselection"].inCategory) // fails 2jet preselection enters 01 categories
    {
        if(vars.recoCandPtPF > cDimuPtMin01T){ categoryMap["01_Jet_Tight"].inCategory = true; return; }
        else{ categoryMap["01_Jet_Loose"].inCategory = true; return; }
    }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________CategorySelectionFEWZ__________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void CategorySelectionFEWZ::initCategoryMap()
{
// Initialize the categories
    categoryMap["Wide"] = Category("Wide");
    categoryMap["Narrow"] = Category("Narrow");
    categoryMap["Central_Central"] = Category("Central_Central");
    categoryMap["Central_Not_Central"] = Category("Central_Not_Central");
    categoryMap["1Jet"] = Category("1Jet");
    categoryMap["Central_Central_Wide"] = Category("Central_Central_Wide");
    categoryMap["Central_Not_Central_Wide"] = Category("Central_Not_Central_Wide");
    categoryMap["1Jet_Wide"] = Category("1Jet_Wide");
    categoryMap["Central_Central_Narrow"] = Category("Central_Central_Narrow");
    categoryMap["Central_Not_Central_Narrow"] = Category("Central_Not_Central_Narrow");
    categoryMap["1Jet_Narrow"] = Category("1Jet_Narrow");
}

CategorySelectionFEWZ::CategorySelectionFEWZ()
{
// Standard values for the FEWZ categorization

    initCategoryMap();

    // init cut values
    cMassSplit = 160;
    cEtaCentralSplit = 0.8;
    cJetPtMin = 30;
    cJetEtaMax = 4.7;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelectionFEWZ::CategorySelectionFEWZ(float massSplit, float etaCentralSplit, float jetPtMin, float jetEtaMax)
{
// Custom values for the cuts

    initCategoryMap();

    // init cut values
    cMassSplit = massSplit;
    cEtaCentralSplit = etaCentralSplit;
    cJetPtMin = jetPtMin;
    cJetEtaMax = jetEtaMax;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void CategorySelectionFEWZ::evaluate(VarSet& vars)
{
// Determine which category the event belongs to

    // Should cut out all events that don't fall into the wide mass window in earlier selection stage
    // All events that pass are in window of min to max
    categoryMap["Wide"].inCategory = true;

    // Narrow goes from min to cMassSplit
    if(vars.recoCandMass < cMassSplit) categoryMap["Narrow"].inCategory = true;

    // Both central
    if(TMath::Abs(vars.reco1.eta) < 0.8 && TMath::Abs(vars.reco2.eta) < 0.8) categoryMap["Central_Central"].inCategory = true;

    // Not both, but at least one is central
    else if(TMath::Abs(vars.reco1.eta) < 0.8 || TMath::Abs(vars.reco2.eta) < 0.8) categoryMap["Central_Not_Central"].inCategory = true;

    // One category that passes basic selections and has exactly one jet
    if(vars.validJets.size() == 1) categoryMap["1Jet"].inCategory = true; 

    // Final Categories ///////////////////////////////////////////////////////
    if(categoryMap["Wide"].inCategory && categoryMap["Central_Central"].inCategory) categoryMap["Central_Central_Wide"].inCategory = true;
    if(categoryMap["Narrow"].inCategory && categoryMap["Central_Central"].inCategory) categoryMap["Central_Central_Narrow"].inCategory = true;

    if(categoryMap["Wide"].inCategory && categoryMap["Central_Not_Central"].inCategory) categoryMap["Central_Not_Central_Wide"].inCategory = true;
    if(categoryMap["Narrow"].inCategory && categoryMap["Central_Not_Central"].inCategory) categoryMap["Central_Not_Central_Narrow"].inCategory = true;

    if(categoryMap["Wide"].inCategory && categoryMap["1Jet"].inCategory) categoryMap["1Jet_Wide"].inCategory = true;
    if(categoryMap["Narrow"].inCategory && categoryMap["1Jet"].inCategory) categoryMap["1Jet_Narrow"].inCategory = true;

    return;
}
