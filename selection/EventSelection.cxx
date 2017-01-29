///////////////////////////////////////////////////////////////////////////
//                             EventSelection.cxx                        //
//=======================================================================//
//                                                                       //
//        Some selections for the event. Consider the event in the       //
//        if it passes.                                                  //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "EventSelection.h"
#include "ParticleTools.h"
#include "TMath.h"
#include <iostream>

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________Run1EventSelectionCuts80X______________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

// Event selection based on the run1 h2mu analysis, still using these
// for the run2 analysis modified for our hlt trigger matching
// criteria

// The 80X samples didn't have trigger matching, so the data requires
// hlt trigger matching but the mc doesn't, then the mc is scaled

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

    //cutset.cuts[0].value = vars.recoMuons->at(0).charge != vars.recoMuons->at(1).charge;
    //cutset.cuts[1].value = -1;
    //cutset.cuts[2].value = vars.dimuCand.mass_PF;

    // Set cuts[1].value correctly
    //if(vars.recoMuons->at(0).pt >= vars.recoMuons->at(1).pt) cutset.cuts[1].value = vars.recoMuons->at(0).pt;
    //else cutset.cuts[1].value = vars.recoMuons->at(1).pt;

    if(!(vars.dimuCand.mass_PF > cDimuMassMin) && cutset.cuts[2].on) return false;

    // Pt cuts on leading and subleading muons
    if(!cutset.cuts[1].on) ;
    else if(vars.recoMuons->at(0).pt >= vars.recoMuons->at(1).pt && TMath::Abs(vars.recoMuons->at(0).pt) > cTrigMuPtMin) ;
    else if(vars.recoMuons->at(0).pt < vars.recoMuons->at(1).pt && TMath::Abs(vars.recoMuons->at(1).pt) > cTrigMuPtMin) ;
    else
    {
         return false;
    }

    // trigger matching
    if(isData)
    {
        if(!(vars.recoMuons->at(0).isHltMatched[4] || vars.recoMuons->at(0).isHltMatched[5] || 
             vars.recoMuons->at(1).isHltMatched[4] || vars.recoMuons->at(1).isHltMatched[5]))
            return false;
    }

    // Require oppositely charged muons in the event
    if(!(vars.recoMuons->at(0).charge != vars.recoMuons->at(1).charge) && cutset.cuts[0].on) return false;

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
    cutset.cuts[0].tstring = "recoMuons->at(0).charge != recoMuons->at(1).charge";
    cutset.cuts[0].bins = 2;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 2;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = Form("leadingRecoMu.pt > %6.2f", cTrigMuPtMin);
    cutset.cuts[1].tstring.Form("max(recoMuons->at(0).pt, recoMuons->at(1).pt) > %6.2f", cTrigMuPtMin);
    cutset.cuts[1].bins = 201;
    cutset.cuts[1].min = -1;
    cutset.cuts[1].max = 200;
    cutset.cuts[1].cutvalue = &cTrigMuPtMin;
    cutset.cuts[1].ismin = true;

    cutset.cuts[2].name = Form("mass_PF > %6.2f", cDimuMassMin);
    cutset.cuts[2].tstring.Form("dimuCand.mass_PF > %6.2f", cDimuMassMin);
    cutset.cuts[2].bins = 150;
    cutset.cuts[2].min = 50;
    cutset.cuts[2].max = 200;
    cutset.cuts[2].cutvalue = &cDimuMassMin;
    cutset.cuts[2].ismin = true;

    cutset.cuts[3].name = "recoMu0 or recoMu1 matched to (hlt4 or hlt5)";
    cutset.cuts[3].tstring = "(recoMuons->at(0).isHltMatched[4] || recoMuons->at(0).isHltMatched[5] || " + 
                             TString("recoMuons->at(1).isHltMatched[4] || recoMuons->at(1).isHltMatched[5])");
    cutset.cuts[3].bins = 2;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 2;
    cutset.cuts[3].ismin = true;

    cutset.concatCuts();
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

    cutset.cuts[0].value = vars.recoMuons->at(0).charge != vars.recoMuons->at(1).charge;
    cutset.cuts[1].value = -1;
    cutset.cuts[2].value = vars.dimuCand.mass_PF;

    // Set cuts[1].value correctly
    if(vars.recoMuons->at(0).isHltMatched[4]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons->at(0).pt);
    else if(vars.recoMuons->at(0).isHltMatched[5]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons->at(0).pt);
    else if(vars.recoMuons->at(1).isHltMatched[4]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons->at(1).pt);
    else if(vars.recoMuons->at(1).isHltMatched[5]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons->at(1).pt);
    else
    {
         cutset.cuts[1].value = -1;
    }

    // if the event fails a single cut return false
    // Require oppositely charged muons in the event
    if(!(vars.recoMuons->at(0).charge != vars.recoMuons->at(1).charge) && cutset.cuts[0].on) return false;

    // One muon in the pair must pass one of the HLT triggers. This muon have the appropriate pt and eta.
    // Should probably make this into a function so that we can look at a larger number of triggers without cluttering this too much.
    if(!cutset.cuts[1].on) ;
    else if(vars.recoMuons->at(0).isHltMatched[4] && TMath::Abs(vars.recoMuons->at(0).pt) > cTrigMuPtMin) ;
    else if(vars.recoMuons->at(0).isHltMatched[5] && TMath::Abs(vars.recoMuons->at(0).pt) > cTrigMuPtMin) ;
    else if(vars.recoMuons->at(1).isHltMatched[4] && TMath::Abs(vars.recoMuons->at(1).pt) > cTrigMuPtMin) ;
    else if(vars.recoMuons->at(1).isHltMatched[5] && TMath::Abs(vars.recoMuons->at(1).pt) > cTrigMuPtMin) ;
    else
    {
         return false;
    }

    if(!(vars.dimuCand.mass_PF > cDimuMassMin) && cutset.cuts[2].on) return false;

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
    cutset.cuts[0].name = "recoMuons->at(0).charge != recoMuons->at(1).charge";
    cutset.cuts[0].tstring = "recoMuons->at(0).charge != recoMuons->at(1).charge";
    cutset.cuts[0].bins = 2;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 2;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = "trigMatchedRecoMu.pt";
    cutset.cuts[1].tstring.Form("recoMuons->at(0).isHltMatched[4||5] && recoMuons->at(0).pt > %5.2f", cTrigMuPtMin);
    cutset.cuts[1].bins = 201;
    cutset.cuts[1].min = -1;
    cutset.cuts[1].max = 200;
    cutset.cuts[1].cutvalue = &cTrigMuPtMin;
    cutset.cuts[1].ismin = true;

    cutset.cuts[2].name = "dimuCand.mass_PF_Min";
    cutset.cuts[2].tstring.Form("dimuCand.mass_PF > %5.2f", cDimuMassMin);
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

// some event selection criteria also based on the Run1Event selection
// use these to calculate the significance in a window around the higgs
// mass. So these selections have an upper and lower bound on the mass

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

    cutset.cuts[0].value = vars.recoMuons->at(0).charge != vars.recoMuons->at(1).charge;
    cutset.cuts[1].value = -1;
    cutset.cuts[2].value = vars.dimuCand.mass_PF;
    cutset.cuts[3].value = vars.dimuCand.mass_PF;

    // Set cuts[1].value correctly
    if(vars.recoMuons->at(0).isHltMatched[4]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons->at(0).pt);
    else if(vars.recoMuons->at(0).isHltMatched[5]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons->at(0).pt);
    else if(vars.recoMuons->at(1).isHltMatched[4]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons->at(1).pt);
    else if(vars.recoMuons->at(1).isHltMatched[5]) cutset.cuts[1].value = TMath::Abs(vars.recoMuons->at(1).pt);
    else
    {
         cutset.cuts[1].value = -1;
    }

    // if the event fails a single cut return false
    // Require oppositely charged muons in the event
    if(!(vars.recoMuons->at(0).charge != vars.recoMuons->at(1).charge) && cutset.cuts[0].on) return false;

    // One muon in the pair must pass one of the HLT triggers. This muon have the appropriate pt and eta.
    // Should probably make this into a function so that we can look at a larger number of triggers without cluttering this too much.
    if(!cutset.cuts[1].on) ;
    else if(vars.recoMuons->at(0).isHltMatched[4] && TMath::Abs(vars.recoMuons->at(0).pt) > cTrigMuPtMin) ;
    else if(vars.recoMuons->at(0).isHltMatched[5] && TMath::Abs(vars.recoMuons->at(0).pt) > cTrigMuPtMin) ;
    else if(vars.recoMuons->at(1).isHltMatched[4] && TMath::Abs(vars.recoMuons->at(1).pt) > cTrigMuPtMin) ;
    else if(vars.recoMuons->at(1).isHltMatched[5] && TMath::Abs(vars.recoMuons->at(1).pt) > cTrigMuPtMin) ;
    else
    {
         return false;
    }

    if(!(vars.dimuCand.mass_PF > cDimuMassMin) && cutset.cuts[2].on) return false;
    if(!(vars.dimuCand.mass_PF < cDimuMassMax) && cutset.cuts[3].on) return false;

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

    cutset.cuts[0].name = "recoMuons->at(0).charge != recoMuons->at(1).charge";
    cutset.cuts[0].tstring = "recoMuons->at(0).charge != recoMuons->at(1).charge";
    cutset.cuts[0].bins = 2;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 2;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = "trigMatchedRecoMu.pt";
    cutset.cuts[1].tstring.Form("recoMuons->at(0).isHltMatched[4||5] && recoMuons->at(0).pt > %5.2f", cTrigMuPtMin);
    cutset.cuts[1].bins = 201;
    cutset.cuts[1].min = -1;
    cutset.cuts[1].max = 200;
    cutset.cuts[1].cutvalue = &cTrigMuPtMin;
    cutset.cuts[1].ismin = true;

    cutset.cuts[2].name = "dimuCand.mass_PF_Min";
    cutset.cuts[2].tstring.Form("dimuCand.mass_PF > %5.2f", cDimuMassMin);
    cutset.cuts[2].bins = 150;
    cutset.cuts[2].min = 50;
    cutset.cuts[2].max = 200;
    cutset.cuts[2].cutvalue = &cDimuMassMin;
    cutset.cuts[2].ismin = true;

    cutset.cuts[3].name = "dimuCand.mass_PF_Max";
    cutset.cuts[3].tstring.Form("dimuCand.mass_PF < %5.2f", cDimuMassMax);
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
// _______________________FEWZCompareCuts________________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

// Yet another set of event selection criteria. These are used
// in FEWZ (matrix element based "generator" used to model Drell Yan) 
// comparison studies

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
        leadPt = TMath::Max(TMath::Abs(vars.recoMuons->at(0).pt), TMath::Abs(vars.recoMuons->at(1).pt));
        subleadPt = TMath::Min(TMath::Abs(vars.recoMuons->at(0).pt), TMath::Abs(vars.recoMuons->at(1).pt));
        eta0 = vars.recoMuons->at(0).eta;
        eta1 = vars.recoMuons->at(1).eta;
        dimu_mass = vars.dimuCand.mass_PF;
        charge0 = vars.recoMuons->at(0).charge;
        charge1 = vars.recoMuons->at(1).charge;
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
        if(!(vars.recoMuons->at(0).iso() <= cMaxRelIso) ) return false;
    }
    if(cutset.cuts[8].on)
    {
        if(!(vars.recoMuons->at(1).iso() <= cMaxRelIso) ) return false;
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

    cutset.cuts[2].name = "recoMuons->at(0).eta";
    cutset.cuts[2].tstring.Form("TMath::Abs(recoMuons->at(0).eta) < %4.2f", cMaxEta);
    cutset.cuts[2].bins = 50;
    cutset.cuts[2].min = 0;
    cutset.cuts[2].max = 3;
    cutset.cuts[2].cutvalue = &cMaxEta;
    cutset.cuts[2].ismin = false;

    cutset.cuts[3].name = "recoMuons->at(1).eta";
    cutset.cuts[3].tstring.Form("TMath::Abs(recoMuons->at(1).eta) < %4.2f", cMaxEta);
    cutset.cuts[3].bins = 50;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 3;
    cutset.cuts[3].cutvalue = &cMaxEta;
    cutset.cuts[3].ismin = false;

    cutset.cuts[4].name = "dimuCand.mass_PF_Min";
    cutset.cuts[4].tstring.Form("dimuCand.mass_PF < %5.2f", cDimuMassMin);
    cutset.cuts[4].bins = 100;
    cutset.cuts[4].min = 110;
    cutset.cuts[4].max = 310;
    cutset.cuts[4].cutvalue = &cDimuMassMin;
    cutset.cuts[4].ismin = true;

    cutset.cuts[5].name = "dimuCand.mass_PF_Max";
    cutset.cuts[5].tstring.Form("dimuCand.mass_PF < %5.2f", cDimuMassMax);
    cutset.cuts[5].bins = 100;
    cutset.cuts[5].min = 110;
    cutset.cuts[5].max = 310;
    cutset.cuts[5].cutvalue = &cDimuMassMax;
    cutset.cuts[5].ismin = false;

    cutset.cuts[6].name = "recoMuons->at(0).charge != recoMuons->at(1).charge";
    cutset.cuts[6].tstring = "recoMuons->at(0).charge != recoMuons->at(1).charge";
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

