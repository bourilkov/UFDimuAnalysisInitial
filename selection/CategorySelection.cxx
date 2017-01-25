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
// _______________________CategorySelectionRun1__________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void CategorySelectionRun1::initCategoryMap()
{
// Initialize the categories
    categoryMap["ALL"] = Category("ALL", true);

    // intermediate categories to make things easier
    categoryMap["2_Jet"] = Category("2_Jet", true);
    categoryMap["01_Jet"] = Category("01_Jet", true);

    categoryMap["2_Jet_VBF_Tight"] = Category("2_Jet_VBF_Tight");
    categoryMap["2_Jet_VBF_Loose"] = Category("2_Jet_VBF_Loose");
    categoryMap["2_Jet_GGF_Tight"] = Category("2_Jet_GGF_Tight");
    categoryMap["01_Jet_Tight"] = Category("01_Jet_Tight", true);
    categoryMap["01_Jet_Loose"] = Category("01_Jet_Loose", true);

    // intermediate categories to make things easier, don't plot these, hence the true
    categoryMap["BB"] = Category("BB", true);
    categoryMap["BO"] = Category("BO", true);
    categoryMap["BE"] = Category("BE", true);
    categoryMap["OO"] = Category("OO", true);
    categoryMap["OE"] = Category("OE", true);
    categoryMap["EE"] = Category("EE", true);

    categoryMap["01_Jet_Tight_BB"] = Category("01_Jet_Tight_BB");
    categoryMap["01_Jet_Tight_BO"] = Category("01_Jet_Tight_BO");
    categoryMap["01_Jet_Tight_BE"] = Category("01_Jet_Tight_BE");
    categoryMap["01_Jet_Tight_OO"] = Category("01_Jet_Tight_OO");
    categoryMap["01_Jet_Tight_OE"] = Category("01_Jet_Tight_OE");
    categoryMap["01_Jet_Tight_EE"] = Category("01_Jet_Tight_EE");

    categoryMap["01_Jet_Loose_BB"] = Category("01_Jet_Loose_BB");
    categoryMap["01_Jet_Loose_BO"] = Category("01_Jet_Loose_BO");
    categoryMap["01_Jet_Loose_BE"] = Category("01_Jet_Loose_BE");
    categoryMap["01_Jet_Loose_OO"] = Category("01_Jet_Loose_OO");
    categoryMap["01_Jet_Loose_OE"] = Category("01_Jet_Loose_OE");
    categoryMap["01_Jet_Loose_EE"] = Category("01_Jet_Loose_EE");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

CategorySelectionRun1::CategorySelectionRun1()
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

CategorySelectionRun1::CategorySelectionRun1(float leadPtMin, float subleadPtMin, float METMax, float dijetMassMinVBFT, float dijetDeltaEtaMinVBFT, float dijetMassMinGGFT,
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

void CategorySelectionRun1::evaluate(VarSet& vars)
{
// Determine which category the event belongs to

    // Inclusive category, all events that passed the selection cuts
    categoryMap["ALL"].inCategory = true;

    // Geometric Categories
    // Barrel Barrel
    if(TMath::Abs(vars.recoMuons.eta[0]) < 0.8 && TMath::Abs(vars.recoMuons.eta[1]) < 0.8) 
        categoryMap["BB"].inCategory = true;
    // Overlap Overlap
    if(TMath::Abs(vars.recoMuons.eta[0])>=0.8 && TMath::Abs(vars.recoMuons.eta[0])<1.6 && TMath::Abs(vars.recoMuons.eta[1])>=0.8 && TMath::Abs(vars.recoMuons.eta[1])<1.6) 
        categoryMap["OO"].inCategory = true;
    // Endcap Endcap
    if(TMath::Abs(vars.recoMuons.eta[0]) >= 1.6 && TMath::Abs(vars.recoMuons.eta[1]) >= 1.6) 
        categoryMap["EE"].inCategory = true;

    // Barrel Overlap
    if(TMath::Abs(vars.recoMuons.eta[0]) < 0.8 && TMath::Abs(vars.recoMuons.eta[1]) >= 0.8 && TMath::Abs(vars.recoMuons.eta[1]) < 1.6) 
        categoryMap["BO"].inCategory = true;
    if(TMath::Abs(vars.recoMuons.eta[1]) < 0.8 && TMath::Abs(vars.recoMuons.eta[0]) >= 0.8 && TMath::Abs(vars.recoMuons.eta[0]) < 1.6) 
        categoryMap["BO"].inCategory = true;

    // Barrel Endcap
    if(TMath::Abs(vars.recoMuons.eta[0]) < 0.8 && TMath::Abs(vars.recoMuons.eta[1]) >= 1.6) 
        categoryMap["BE"].inCategory = true;
    if(TMath::Abs(vars.recoMuons.eta[1]) < 0.8 && TMath::Abs(vars.recoMuons.eta[0]) >= 1.6) 
        categoryMap["BE"].inCategory = true;

    // Overlap Endcap
    if(TMath::Abs(vars.recoMuons.eta[0]) >= 0.8 && TMath::Abs(vars.recoMuons.eta[0]) < 1.6 && TMath::Abs(vars.recoMuons.eta[1]) >= 1.6) 
        categoryMap["OE"].inCategory = true;
    if(TMath::Abs(vars.recoMuons.eta[1]) >= 0.8 && TMath::Abs(vars.recoMuons.eta[1]) < 1.6 && TMath::Abs(vars.recoMuons.eta[0]) >= 1.6) 
        categoryMap["OE"].inCategory = true;

    // jet category selection
    if(vars.validJets.size() >= 2)
    {
        TLorentzVector leadJet = vars.validJets[0];
        TLorentzVector subleadJet = vars.validJets[1];
        TLorentzVector dijet = leadJet + subleadJet;

        float dEta = leadJet.Eta() - subleadJet.Eta();
        float dijetMass = dijet.M();

        if(leadJet.Pt() > cLeadPtMin && subleadJet.Pt() > cSubleadPtMin && vars.met.pt < cMETMax)
        {
            categoryMap["2_Jet"].inCategory = true;
            if(dijetMass > cDijetMassMinVBFT && TMath::Abs(dEta) > cDijetDeltaEtaMinVBFT){ categoryMap["2_Jet_VBF_Tight"].inCategory = true; return; }
            else if(dijetMass > cDijetMassMinGGFT && vars.dimuCand.recoCandPtPF > cDimuPtMinGGFT){ categoryMap["2_Jet_GGF_Tight"].inCategory = true; return; }
            else{ categoryMap["2_Jet_VBF_Loose"].inCategory = true; return; }
        }
    }
    if(!categoryMap["2_Jet"].inCategory) // fails 2jet preselection enters 01 categories
    {
        categoryMap["01_Jet"].inCategory = true;
        if(vars.dimuCand.recoCandPtPF > cDimuPtMin01T){ categoryMap["01_Jet_Tight"].inCategory = true;}
        else{ categoryMap["01_Jet_Loose"].inCategory = true; }

        // Geometric categories for 01_Jet categories
        // tight
        if(categoryMap["01_Jet_Tight"].inCategory && categoryMap["BB"].inCategory) categoryMap["01_Jet_Tight_BB"].inCategory = true;
        if(categoryMap["01_Jet_Tight"].inCategory && categoryMap["BO"].inCategory) categoryMap["01_Jet_Tight_BO"].inCategory = true;
        if(categoryMap["01_Jet_Tight"].inCategory && categoryMap["BE"].inCategory) categoryMap["01_Jet_Tight_BE"].inCategory = true;
        if(categoryMap["01_Jet_Tight"].inCategory && categoryMap["OO"].inCategory) categoryMap["01_Jet_Tight_OO"].inCategory = true;
        if(categoryMap["01_Jet_Tight"].inCategory && categoryMap["OE"].inCategory) categoryMap["01_Jet_Tight_OE"].inCategory = true;
        if(categoryMap["01_Jet_Tight"].inCategory && categoryMap["EE"].inCategory) categoryMap["01_Jet_Tight_EE"].inCategory = true;

        // loose
        if(categoryMap["01_Jet_Loose"].inCategory && categoryMap["BB"].inCategory) categoryMap["01_Jet_Loose_BB"].inCategory = true;
        if(categoryMap["01_Jet_Loose"].inCategory && categoryMap["BO"].inCategory) categoryMap["01_Jet_Loose_BO"].inCategory = true;
        if(categoryMap["01_Jet_Loose"].inCategory && categoryMap["BE"].inCategory) categoryMap["01_Jet_Loose_BE"].inCategory = true;
        if(categoryMap["01_Jet_Loose"].inCategory && categoryMap["OO"].inCategory) categoryMap["01_Jet_Loose_OO"].inCategory = true;
        if(categoryMap["01_Jet_Loose"].inCategory && categoryMap["OE"].inCategory) categoryMap["01_Jet_Loose_OE"].inCategory = true;
        if(categoryMap["01_Jet_Loose"].inCategory && categoryMap["EE"].inCategory) categoryMap["01_Jet_Loose_EE"].inCategory = true;
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

    // intermediate categories to make things easier, don't plot these, hence the true
    categoryMap["Central_Central"] = Category("Central_Central", true);
    categoryMap["Central_Not_Central"] = Category("Central_Not_Central", true);

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

    float dimu_mass = vars.dimuCand.recoCandMass;
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

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________LotsOfCategoriesRun2__________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void LotsOfCategoriesRun2::initCategoryMap()
{
// Initialize the categories

    ///////////////// INCLUSIVE //////////////////////////////
    categoryMap["ALL"] = Category("ALL");

    ///////////////// GEOMETRY //////////////////////////////
    categoryMap["BB"] = Category("BB", true);
    categoryMap["BO"] = Category("BO", true);
    categoryMap["BE"] = Category("BE", true);
    categoryMap["OO"] = Category("OO", true);
    categoryMap["OE"] = Category("OE", true);
    categoryMap["EE"] = Category("EE", true);

    ///////////////// PRESELECTION //////////////////////////////
    categoryMap["Preselection_Pass"] = Category("Preselection_Pass");

      ///////////////// AT LEAST ONE B-JET //////////////////////////////
      categoryMap["1b"] = Category("1b");

        ///////////////// TTH 2 EXTRA LEPTONS //////////////////////////////
        categoryMap["1b_TTH"]        = Category("1b_TTH");
        categoryMap["1b_TTH_2e"]     = Category("1b_TTH_2e",     true);
        categoryMap["1b_TTH_1e_1mu"] = Category("1b_TTH_1e_1mu", true);
        categoryMap["1b_TTH_2mu"]    = Category("1b_TTH_2mu",    true);
  
        ///////////////// TTH/BBH 0 EXTRA LEPTONS (HADRONS)///////////////////
        categoryMap["1b_TTH_BBH"]              = Category("1b_TTH_BBH");
        categoryMap["1b_TTH_BBH_Tight"]        = Category("1b_TTH_BBH_Tight",        true);
        categoryMap["1b_TTH_BBH_V_Hadronic_H"] = Category("1b_TTH_BBH_V_Hadronic_H", true);
  
        ///////////////// 1 B-JET EVENTS THAT DONT FIT ELSEWHERE////////////////
        categoryMap["1b_Leftovers"] = Category("1b_Leftovers");
  
      ///////////////// NO B-JETS //////////////////////////////
      categoryMap["0b"] = Category("0b");
  
        ///////////////// NOT V(lept)H (VBF, gF, V(had)H, ZvvH) /////////////////////
        categoryMap["0b_nonVlH"] = Category("0b_nonVlH");

        ///////////////// 2-jet (VBF, V(had)H, gF) //////////////////////////////
        categoryMap["0b_nonVlH_2j"]              = Category("0b_nonVlH_2j");
        categoryMap["0b_nonVlH_2j_VBF_Tight"]    = Category("0b_nonVlH_2j_VBF_Tight");
        categoryMap["0b_nonVlH_2j_VBF_Loose"]    = Category("0b_nonVlH_2j_VBF_Loose");
        categoryMap["0b_nonVlH_2j_V_Hadronic_H"] = Category("0b_nonVlH_2j_V_Hadronic_H");
        categoryMap["0b_nonVlH_2j_gF"]           = Category("0b_nonVlH_2j_gF");

        ///////////////// 01-jet (gF Tight, gF Loose, ZvvH) //////////////////////////////
        categoryMap["0b_nonVlH_01j"]          = Category("0b_nonVlH_01j");

          // ZvvH
          categoryMap["0b_nonVlH_01j_ZvvH"]     = Category("0b_nonVlH_01j_ZvvH");

          // gF Tight
          categoryMap["0b_nonVlH_01j_gF_Tight"] = Category("0b_nonVlH_01j_gF_Tight");

            // gF Tight Geometrized
            categoryMap["0b_nonVlH_01j_gF_Tight_BB"] = Category("0b_nonVlH_01j_gF_Tight_BB", true);
            categoryMap["0b_nonVlH_01j_gF_Tight_BO"] = Category("0b_nonVlH_01j_gF_Tight_BO", true);
            categoryMap["0b_nonVlH_01j_gF_Tight_BE"] = Category("0b_nonVlH_01j_gF_Tight_BE", true);
            categoryMap["0b_nonVlH_01j_gF_Tight_OO"] = Category("0b_nonVlH_01j_gF_Tight_OO", true);
            categoryMap["0b_nonVlH_01j_gF_Tight_OE"] = Category("0b_nonVlH_01j_gF_Tight_OE", true);
            categoryMap["0b_nonVlH_01j_gF_Tight_EE"] = Category("0b_nonVlH_01j_gF_Tight_EE", true);

          // gF Loose
          categoryMap["0b_nonVlH_01j_gF_Loose"] = Category("0b_nonVlH_01j_gF_Loose");

            // gF Loose Geometrized
            categoryMap["0b_nonVlH_01j_gF_Loose_BB"] = Category("0b_nonVlH_01j_gF_Loose_BB", true);
            categoryMap["0b_nonVlH_01j_gF_Loose_BO"] = Category("0b_nonVlH_01j_gF_Loose_BO", true);
            categoryMap["0b_nonVlH_01j_gF_Loose_BE"] = Category("0b_nonVlH_01j_gF_Loose_BE", true);
            categoryMap["0b_nonVlH_01j_gF_Loose_OO"] = Category("0b_nonVlH_01j_gF_Loose_OO", true);
            categoryMap["0b_nonVlH_01j_gF_Loose_OE"] = Category("0b_nonVlH_01j_gF_Loose_OE", true);
            categoryMap["0b_nonVlH_01j_gF_Loose_EE"] = Category("0b_nonVlH_01j_gF_Loose_EE", true);

        ///////////////// V(lept)H (...) /////////////////////
        categoryMap["0b_VlH"]    = Category("0b_VlH");

            // VlH according to the V decays
            categoryMap["0b_VlH_We"]         = Category("0b_VlH_We");
            categoryMap["0b_VlH_Wmu"]        = Category("0b_VlH_Wmu");
            categoryMap["0b_VlH_Ztautau"]        = Category("0b_VlH_Ztautau");
            categoryMap["0b_VlH_Zmumu"]      = Category("0b_VlH_Zmumu");
            categoryMap["0b_VlH_Zee"]        = Category("0b_VlH_Zee");
            categoryMap["0b_VlH_Leftovers"]  = Category("0b_VlH_Leftovers");

    ///////////////// FAIL PRESELECTION //////////////////////////////
    categoryMap["Preselection_Fail"] = Category("Preselection_Fail");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

LotsOfCategoriesRun2::LotsOfCategoriesRun2()
{
// initialize the default cut values in the constructor

    initCategoryMap();

        // Preselection
        c_pre_numExtraLeptonsMax = 2;
        c_pre_useTau = false;
        c_pre_numBJetsMin = 1;

        // b-jet categories
        c_1b_numExtraLeptons_tth     = 2; 
        c_1b_numExtraLeptons_tth_bbh = 0; 

            // ttH categories
            c_tth_numExtraElectrons = -999;

            // tth-bbh categories
            c_tth_bbh_mbb      = -999;
            c_tth_bbh_mt_bMET  = -999;
            c_tth_bbh_MET      = -999;

        // right side of b-jet test, fewer b-jets
        // 0 extra leptons, 1 or 2 extra leptons
        c_0b_numExtraLeptonsMin = 1; 

            // nonVlH (0 extra leptons)
            c_0b_nonVlH_njetsMin = 2;
    
                // >= c_0b_nonVlH_njets
                c_0b_nonVlH_2j_mjj_min_vbfTight    = 500;
                c_0b_nonVlH_2j_dEtajj_min_vbfTight = 2.5;

                c_0b_nonVlH_2j_mjj_min_vbfLoose    = 250;
                c_0b_nonVlH_2j_dEtajj_min_vbfLoose = 2.5;

                c_0b_nonVlH_2j_mjj_min_VhH = 60;
                c_0b_nonVlH_2j_mjj_max_VhH = 110;
                c_0b_nonVlH_2j_dEtajjMuMu_max_VhH = 1.5;

                // < c_0b_nonVlH_njets
                c_0b_nonVlH_01j_MET_min_ZvvH = 40;
                c_0b_nonVlH_01j_dimuPt_min_gfTight = 25;

                // muon geometry categorization for gf categories
                c_geo_bmax = 0.8;
                c_geo_omax = 1.6;
                c_geo_emax = 2.4;

            // VlH (1 or 2 extra leptons)
            c_0b_VlH_MET_min = 40;

                // VlH according to the V decays
                c_0b_VlH_We_num_e   = 1;         
                c_0b_VlH_We_num_mu  = 0;         

                c_0b_VlH_Wmu_num_e  = 0;        
                c_0b_VlH_Wmu_num_mu = 1;        

                c_0b_VlH_Ztautau_num_e  = 1;        
                c_0b_VlH_Ztautau_num_mu = 1;        

                c_0b_VlH_Zmumu_num_e  = 0;
                c_0b_VlH_Zmumu_num_mu = 2;

                c_0b_VlH_Zee_num_e  = 2; 
                c_0b_VlH_Zee_num_mu = 0; 

}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void LotsOfCategoriesRun2::evaluateMuonGeometry(VarSet& vars)
{
// Determine which geometric category the event belongs to

    // Geometric Categories
    // Barrel Barrel
    if(TMath::Abs(vars.recoMuons.eta[0]) < c_geo_bmax && TMath::Abs(vars.recoMuons.eta[1]) < c_geo_bmax) 
        categoryMap["BB"].inCategory = true;
    // Overlap Overlap
    if(TMath::Abs(vars.recoMuons.eta[0])>=c_geo_bmax && TMath::Abs(vars.recoMuons.eta[0])<c_geo_omax && TMath::Abs(vars.recoMuons.eta[1])>=c_geo_bmax && TMath::Abs(vars.recoMuons.eta[1])<c_geo_omax) 
        categoryMap["OO"].inCategory = true;
    // Endcap Endcap
    if(TMath::Abs(vars.recoMuons.eta[0]) >= c_geo_omax && TMath::Abs(vars.recoMuons.eta[1]) >= c_geo_omax) 
        categoryMap["EE"].inCategory = true;

    // Barrel Overlap
    if(TMath::Abs(vars.recoMuons.eta[0]) < c_geo_bmax && TMath::Abs(vars.recoMuons.eta[1]) >= c_geo_bmax && TMath::Abs(vars.recoMuons.eta[1]) < c_geo_omax) 
        categoryMap["BO"].inCategory = true;
    if(TMath::Abs(vars.recoMuons.eta[1]) < c_geo_bmax && TMath::Abs(vars.recoMuons.eta[0]) >= c_geo_bmax && TMath::Abs(vars.recoMuons.eta[0]) < c_geo_omax) 
        categoryMap["BO"].inCategory = true;

    // Barrel Endcap
    if(TMath::Abs(vars.recoMuons.eta[0]) < c_geo_bmax && TMath::Abs(vars.recoMuons.eta[1]) >= c_geo_omax) 
        categoryMap["BE"].inCategory = true;
    if(TMath::Abs(vars.recoMuons.eta[1]) < c_geo_bmax && TMath::Abs(vars.recoMuons.eta[0]) >= c_geo_omax) 
        categoryMap["BE"].inCategory = true;

    // Overlap Endcap
    if(TMath::Abs(vars.recoMuons.eta[0]) >= c_geo_bmax && TMath::Abs(vars.recoMuons.eta[0]) < c_geo_omax && TMath::Abs(vars.recoMuons.eta[1]) >= c_geo_omax) 
        categoryMap["OE"].inCategory = true;
    if(TMath::Abs(vars.recoMuons.eta[1]) >= c_geo_bmax && TMath::Abs(vars.recoMuons.eta[1]) < c_geo_omax && TMath::Abs(vars.recoMuons.eta[0]) >= c_geo_omax) 
        categoryMap["OE"].inCategory = true;
}
void LotsOfCategoriesRun2::evaluate(VarSet& vars)
{
    ///////////////// INCLUSIVE //////////////////////////////
    categoryMap["ALL"].inCategory = true;

    // figure out bb,oo,ee,bo,be,oe
    ///////////////// MUON GEOMETRY //////////////////////////////
    evaluateMuonGeometry(vars);

    ///////////////// PRESELECTION //////////////////////////////
    if(vars.validExtraMuons.size() + vars.validElectrons.size() <= c_pre_numExtraLeptonsMax) 
        categoryMap["Preselection_Pass"].inCategory = true;
    else
        categoryMap["Preselection_Fail"].inCategory = true;

   // Determine whether we are in the at least 1b-jet categories or 0b-jet categories
   if(categoryMap["Preselection_Pass"].inCategory)
   {
       //std::cout << "    pass preselection..." << std::endl;
       if(vars.validBJets.size() >= c_pre_numBJetsMin) 
           categoryMap["1b"].inCategory = true;
       else
           categoryMap["0b"].inCategory = true;
   }

       ///////////////// 1b CATEGORIES //////////////////////////////
       if(categoryMap["1b"].inCategory)
       {
           //std::cout << "    pass 1b..." << std::endl;
           if(vars.validExtraMuons.size() + vars.validElectrons.size() == c_1b_numExtraLeptons_tth)
               categoryMap["1b_TTH"].inCategory = true;
           else if(vars.validExtraMuons.size() + vars.validElectrons.size() == c_1b_numExtraLeptons_tth_bbh)
               categoryMap["1b_TTH_BBH"].inCategory = true;
           else 
               categoryMap["1b_Leftovers"].inCategory = true;
       }
           ///////////////// 1b-TTH (2 extra lept) CATEGORIES //////////////////////////////
           if(categoryMap["1b_TTH"].inCategory)
           {
               //std::cout << "    pass 1b TTH..." << std::endl;
           }
           ///////////////// 1b_TTH_BBH (0 extra lept) CATEGORIES //////////////////////////////
           if(categoryMap["1b_TTH_BBH"].inCategory)
           {
               //std::cout << "    pass 1b TTH_BBH..." << std::endl;
           }

       ///////////////// 0b CATEGORIES //////////////////////////////
       if(categoryMap["0b"].inCategory)
       {
           //std::cout << "    pass 0b..." << std::endl;
           // output event information  here to debug... we have 70% in this category for N_valid_whatevers in data
           if(vars.validExtraMuons.size() + vars.validElectrons.size() >= c_0b_numExtraLeptonsMin)
               categoryMap["0b_VlH"].inCategory = true;
           else
               categoryMap["0b_nonVlH"].inCategory = true;
       }
           ///////////////// 0b-VlH    (1,2 extra lept) CATEGORIES //////////////////////////////
           if(categoryMap["0b_VlH"].inCategory)
           {
               //std::cout << "    pass 0b_VlH..." << std::endl;
               if(vars.met.pt >= c_0b_VlH_MET_min)
               {
                    if(vars.validElectrons.size() == c_0b_VlH_We_num_e && vars.validExtraMuons.size() == c_0b_VlH_We_num_mu)    
                        categoryMap["0b_VlH_We"].inCategory = true;

                    else if(vars.validElectrons.size() == c_0b_VlH_Wmu_num_e && vars.validExtraMuons.size() == c_0b_VlH_Wmu_num_mu)    
                        categoryMap["0b_VlH_Wmu"].inCategory = true;

                    else if(vars.validElectrons.size() == c_0b_VlH_Ztautau_num_e && vars.validExtraMuons.size() == c_0b_VlH_Ztautau_num_mu)    
                        categoryMap["0b_VlH_Ztautau"].inCategory = true;

                    else    
                        categoryMap["0b_VlH_Leftovers"].inCategory = true;
               }
               else
               {
                    if(vars.validElectrons.size() == c_0b_VlH_Zmumu_num_e && vars.validExtraMuons.size() == c_0b_VlH_Zmumu_num_mu)    
                        categoryMap["0b_VlH_Zmumu"].inCategory = true;

                    else if(vars.validElectrons.size() == c_0b_VlH_Zee_num_e && vars.validExtraMuons.size() == c_0b_VlH_Zee_num_mu)    
                        categoryMap["0b_VlH_Zee"].inCategory = true;

                    else    
                        categoryMap["0b_VlH_Leftovers"].inCategory = true;
               }
           }

           ///////////////// 0b-nonVlH (0 extra lept) CATEGORIES //////////////////////////////
           if(categoryMap["0b_nonVlH"].inCategory)
           {
               //std::cout << "    pass 0b_non_VlH..." << std::endl;
               if(vars.validJets.size() >= c_0b_nonVlH_njetsMin)
                   categoryMap["0b_nonVlH_2j"].inCategory = true;
               else
                   categoryMap["0b_nonVlH_01j"].inCategory = true;
           }
               ///////////////// 0b-nonVlH_2j (0 extra lept) CATEGORIES //////////////////////////////
               if(categoryMap["0b_nonVlH_2j"].inCategory)
               {
                   //std::cout << "    pass 0b_non_VlHi_2j..." << std::endl;
                   TLorentzVector leadJet    = vars.validJets[0];
                   TLorentzVector subleadJet = vars.validJets[1];
                   TLorentzVector dijet      = leadJet + subleadJet;
           
                   float dEta = TMath::Abs(leadJet.Eta() - subleadJet.Eta());
                   float dijetMass = dijet.M();
                   float dEtajjMuMu = TMath::Abs(dijet.Eta() - vars.dimuCand.recoCandEtaPF); 

                   if(dijetMass > c_0b_nonVlH_2j_mjj_min_vbfTight && dEta > c_0b_nonVlH_2j_dEtajj_min_vbfTight)
                       categoryMap["0b_nonVlH_2j_VBF_Tight"].inCategory = true; 

                   else if(dijetMass > c_0b_nonVlH_2j_mjj_min_vbfLoose && dEta > c_0b_nonVlH_2j_dEtajj_min_vbfLoose)
                       categoryMap["0b_nonVlH_2j_VBF_Loose"].inCategory = true; 

                   else if(dijetMass > c_0b_nonVlH_2j_mjj_min_VhH && dijetMass < c_0b_nonVlH_2j_mjj_max_VhH && dEtajjMuMu < c_0b_nonVlH_2j_dEtajjMuMu_max_VhH)
                       categoryMap["0b_nonVlH_2j_V_Hadronic_H"].inCategory = true; 

                   else
                       categoryMap["0b_nonVlH_2j_gF"].inCategory = true; 
               }

               ///////////////// 0b-nonVlH_01j (0 extra lept) CATEGORIES //////////////////////////////
               if(categoryMap["0b_nonVlH_01j"].inCategory)
               {
                   //std::cout << "    pass 0b_nonVlH_01j..." << std::endl;
                   if(vars.met.pt > c_0b_nonVlH_01j_MET_min_ZvvH)
                       categoryMap["0b_nonVlH_01j_ZvvH"].inCategory = true; 

                   else if(vars.dimuCand.recoCandPtPF >= c_0b_nonVlH_01j_dimuPt_min_gfTight)
                       categoryMap["0b_nonVlH_01j_gF_Tight"].inCategory = true; 

                   else
                       categoryMap["0b_nonVlH_01j_gF_Loose"].inCategory = true; 
               }
}
