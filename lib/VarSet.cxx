///////////////////////////////////////////////////////////////////////////
// ======================================================================//
// VarSet.cxx                                                            //
// ======================================================================//
// TTree variables are loaded into these objects. Each sample has        //
// its own VarSet data structure.                                        //
// We also set up a map<TString,FUNCTION> to get feature values based    //
// upon the name of the feature.                                         //
// Has muons, dimuons, jets, electrons, gen muons, etc...                //
// ======================================================================//
///////////////////////////////////////////////////////////////////////////

// most of the methods are short and defined in the header

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "VarSet.h"

///////////////////////////////////////////////////////////////////////////
// _______________________Constructor/Destructor_________________________//
///////////////////////////////////////////////////////////////////////////

VarSet::VarSet() 
{
    varMap["dimu_pt"]     = &VarSet::dimu_pt;
    varMap["mu0_pt"]      = &VarSet::mu1_pt;
    varMap["mu1_pt"]      = &VarSet::mu2_pt;
    varMap["mu0_eta"]     = &VarSet::mu1_eta;
    varMap["mu1_eta"]     = &VarSet::mu2_eta;
    varMap["mu_res_eta"]  = &VarSet::mu_res_eta;
    varMap["phi_star"]    = &VarSet::phi_star;
    
    varMap["jet0_pt"]   = &VarSet::jet0_pt;
    varMap["jet1_pt"]   = &VarSet::jet1_pt;
    varMap["jet0_eta"]  = &VarSet::jet0_eta;
    varMap["jet1_eta"]  = &VarSet::jet1_eta;
    
    varMap["m_jj"]          = &VarSet::m_jj;
    varMap["dEta_jj"]       = &VarSet::dEta_jj;
    varMap["dEta_jj_mumu"]  = &VarSet::dEta_jj_mumu;
    varMap["dPhi_jj_mumu"]  = &VarSet::dPhi_jj_mumu;
    varMap["zep"]           = &VarSet::zep;

    varMap["vbf_jet0_pt"]   = &VarSet::vbf_jet0_pt;
    varMap["vbf_jet1_pt"]   = &VarSet::vbf_jet1_pt;
    varMap["vbf_jet0_eta"]  = &VarSet::vbf_jet0_eta;
    varMap["vbf_jet1_eta"]  = &VarSet::vbf_jet1_eta;
    
    varMap["vbf_m_jj"]          = &VarSet::vbf_m_jj;
    varMap["vbf_dEta_jj"]       = &VarSet::vbf_dEta_jj;
    varMap["vbf_dEta_jj_mumu"]  = &VarSet::vbf_dEta_jj_mumu;
    varMap["vbf_dPhi_jj_mumu"]  = &VarSet::vbf_dPhi_jj_mumu;
    varMap["vbf_zep"]           = &VarSet::vbf_zep;
    
    varMap["bjet0_pt"]   = &VarSet::bjet0_pt;
    varMap["bjet1_pt"]   = &VarSet::bjet1_pt;
    varMap["bjet0_eta"]  = &VarSet::bjet0_eta;
    varMap["bjet1_eta"]  = &VarSet::bjet1_eta;
    
    varMap["m_bb"]     = &VarSet::m_bb;
    varMap["dEta_bb"]  = &VarSet::dEta_bb;
    
    varMap["N_valid_jets"]           = &VarSet::N_valid_jets;
    varMap["N_valid_bjets"]          = &VarSet::N_valid_bjets;
    varMap["N_valid_extra_muons"]    = &VarSet::N_valid_extra_muons;
    varMap["N_valid_electrons"]      = &VarSet::N_valid_electrons;
    varMap["N_valid_extra_leptons"]  = &VarSet::N_valid_extra_leptons;
    
    varMap["MET"]  = &VarSet::MET;
    varMap["MHT"]  = &VarSet::MHT;
    
    varMap["extra_muon0_pt"]  = &VarSet::extra_muon0_pt; 
    varMap["extra_muon1_pt"]  = &VarSet::extra_muon1_pt; 
    varMap["extra_muon0_eta"] = &VarSet::extra_muon0_eta; 
    varMap["extra_muon1_eta"] = &VarSet::extra_muon1_eta; 

    varMap["electron0_pt"]  = &VarSet::electron0_pt;  
    varMap["electron1_pt"]  = &VarSet::electron1_pt; 
    varMap["electron0_eta"] = &VarSet::electron0_eta; 
    varMap["electron1_eta"] = &VarSet::electron1_eta; 

    varMap["mT_b_MET"] = &VarSet::mT_b_MET;

}
