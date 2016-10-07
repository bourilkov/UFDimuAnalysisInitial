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
#include "ParticleTools.h"

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

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelectionFEWZ::CategorySelectionFEWZ(bool useRecoMu, bool useRecoJets)
{
// Standard values for the FEWZ categorization

    initCategoryMap();
    this->useRecoMu = useRecoMu;
    this->useRecoJets = useRecoJets;

    // init cut values
    cMassSplit = 160;
    cEtaCentralSplit = 0.8;
    cJetPtMin = 30;
    cJetEtaMax = 4.7;
}


///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelectionFEWZ::CategorySelectionFEWZ()
{
// Standard values for the FEWZ categorization

    initCategoryMap();
    useRecoMu = false;
    useRecoJets = false;

    // init cut values
    cMassSplit = 160;
    cEtaCentralSplit = 0.8;
    cJetPtMin = 30;
    cJetEtaMax = 4.7;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelectionFEWZ::CategorySelectionFEWZ(bool useRecoMu, bool useRecoJets, float massSplit, float etaCentralSplit, float jetPtMin, float jetEtaMax)
{
// Custom values for the cuts

    initCategoryMap();
    this->useRecoMu = useRecoMu;
    this->useRecoJets = useRecoJets;

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

    float dimu_mass = vars.recoCandMass;
    float eta0 = vars.recoMuons.eta[0];
    float eta1 = vars.recoMuons.eta[1];
    unsigned int njets = vars.validJets.size();

    // use gen jets for categorization
    if(!useRecoJets) njets = vars.validGenJets.size();

    // use gen muons for categorization 
    if(!useRecoMu)
    {
        // get first postFSR DY gen muon
        _TrackInfo mu0 = ParticleTools::getGenMuDY(0, 1, vars);
        // get second postFSR DY gen muon
        _TrackInfo mu1 = ParticleTools::getGenMuDY(1, 1, vars);

        if(mu0.charge == mu1.charge)
            std::cout << "### ERROR: DYmu0.charge == DYmu1.charge in FEWZ CATEGORY SELECTION. Muons decayed from Z or gamma* should have opposite sign." << std::endl;

        eta0 = mu0.eta;
        eta1 = mu1.eta;

        if(mu0.pt > 0 && mu1.pt > 0) dimu_mass = ParticleTools::getMotherPtEtaPhiM(mu0.pt, mu0.eta, mu0.phi, MASS_MUON, mu1.pt, mu1.eta, mu1.phi, MASS_MUON).M();
        else
        {
            std::cout << "### ERROR: gen_dimu_mass < 0 in FEWZ CATEGORY SELECTION" << std::endl;
            std::cout << vars.eventInfo.run << "," << vars.eventInfo.lumi << "," << vars.eventInfo.event << std::endl;
            std::cout << mu0.pt << "," << mu0.eta << "," << mu0.phi << std::endl;
            std::cout << mu1.pt << "," << mu1.eta << "," << mu1.phi << std::endl;
            std::cout << std::endl;
        }


    }

    // Should cut out all events that don't fall into the wide mass window in earlier selection stage
    // All events that pass are in window of min to max
    categoryMap["Wide"].inCategory = true;

    // Narrow goes from min to cMassSplit
    if(dimu_mass < cMassSplit) categoryMap["Narrow"].inCategory = true;

    // Both central
    if(TMath::Abs(eta0) < 0.8 && TMath::Abs(eta1) < 0.8) categoryMap["Central_Central"].inCategory = true;

    // Not both, but at least one is central
    else if(TMath::Abs(eta0) < 0.8 || TMath::Abs(eta1) < 0.8) categoryMap["Central_Not_Central"].inCategory = true;

    // One category that passes basic selections and has exactly one jet
    if(njets == 1) categoryMap["1Jet"].inCategory = true; 

    // Final Categories ///////////////////////////////////////////////////////
    if(categoryMap["Wide"].inCategory && categoryMap["Central_Central"].inCategory) categoryMap["Central_Central_Wide"].inCategory = true;
    if(categoryMap["Narrow"].inCategory && categoryMap["Central_Central"].inCategory) categoryMap["Central_Central_Narrow"].inCategory = true;

    if(categoryMap["Wide"].inCategory && categoryMap["Central_Not_Central"].inCategory) categoryMap["Central_Not_Central_Wide"].inCategory = true;
    if(categoryMap["Narrow"].inCategory && categoryMap["Central_Not_Central"].inCategory) categoryMap["Central_Not_Central_Narrow"].inCategory = true;

    if(categoryMap["Wide"].inCategory && categoryMap["1Jet"].inCategory) categoryMap["1Jet_Wide"].inCategory = true;
    if(categoryMap["Narrow"].inCategory && categoryMap["1Jet"].inCategory) categoryMap["1Jet_Narrow"].inCategory = true;

    return;
}
