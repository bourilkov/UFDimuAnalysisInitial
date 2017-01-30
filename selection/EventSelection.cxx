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
// _______________________Run2EventSelectionCuts80X______________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

// Event selection based on the run1 h2mu analysis, still using these
// for the run2 analysis modified for our hlt trigger matching
// criteria

// The 80X samples didn't have trigger matching, so the data requires
// hlt trigger matching but the mc doesn't, then the mc is scaled

Run2EventSelectionCuts80X::Run2EventSelectionCuts80X()
{
// Default values for the modified run 1 event selection

    isData = 0;
    cTrigMuPtMin = 26; 
    // sublead pt condition is met in muon selection
                                
    cDimuMassMin = 60;
    cutset.cuts = std::vector<CutInfo>(4, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Run2EventSelectionCuts80X::Run2EventSelectionCuts80X(bool isData)
{
// Default values for the modified run 1 event selection

    this->isData = isData;;
    cTrigMuPtMin = 26; 
    // sublead pt condition is met in muon selection
                                
    cDimuMassMin = 60;
    cutset.cuts = std::vector<CutInfo>(4, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Run2EventSelectionCuts80X::Run2EventSelectionCuts80X(float trigMuPtMin, float dimuMassMin)
{
// Default values for the modified run 1 event selection

    isData = 0;
    cTrigMuPtMin = trigMuPtMin;          // >
    // sublead pt condition is met in muon selection

    cDimuMassMin = dimuMassMin;          // >
    cutset.cuts = std::vector<CutInfo>(4, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Run2EventSelectionCuts80X::Run2EventSelectionCuts80X(bool isData, float trigMuPtMin, float dimuMassMin)
{
// Default values for the modified run 1 event selection

    this->isData = isData;
    cTrigMuPtMin = trigMuPtMin;          // >
    // sublead pt condition is met in muon selection

    cDimuMassMin = dimuMassMin;          // >
    cutset.cuts = std::vector<CutInfo>(4, CutInfo());
    makeCutSet();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool Run2EventSelectionCuts80X::evaluate(VarSet& vars)
{
    // check the mass cut first since it's the most likely to cut the event
    // saves computation time so that we don't have to check the others
    if(!(vars.dimuCand->mass_PF > cDimuMassMin) && cutset.cuts[2].on) return false;

    // Pt cuts on leading muon
    int mu1 = vars.dimuCand->iMu1;
    int mu2 = vars.dimuCand->iMu2;
    if(!cutset.cuts[1].on) ;
    else if(vars.recoMuons->at(mu1).pt >= vars.recoMuons->at(mu2).pt && vars.recoMuons->at(mu1).pt > cTrigMuPtMin) ;
    else if(vars.recoMuons->at(mu1).pt <  vars.recoMuons->at(mu2).pt && vars.recoMuons->at(mu2).pt > cTrigMuPtMin) ;
    else
    {
         return false;
    }

    // trigger matching
    if(isData)
    {
        // 2 and 3 should be isoMu24 and isoTkMu24
        if(!(vars.recoMuons->at(mu1).isHltMatched[2] || vars.recoMuons->at(mu1).isHltMatched[3] || 
             vars.recoMuons->at(mu2).isHltMatched[2] || vars.recoMuons->at(mu2).isHltMatched[3]))
            return false;
    }

    // Require oppositely charged muons in the event
    if(!(vars.recoMuons->at(mu1).charge != vars.recoMuons->at(mu2).charge) && cutset.cuts[0].on) return false;

    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString Run2EventSelectionCuts80X::string()
{
    return TString("Run2_Event_Selection_80X");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void Run2EventSelectionCuts80X::makeCutSet()
{
    cutset.cuts[0].name = "recoMu1.charge != recoMu2.charge";
    cutset.cuts[0].bins = 2;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 2;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = Form("leadingRecoMu.pt > %6.2f", cTrigMuPtMin);
    cutset.cuts[1].bins = 201;
    cutset.cuts[1].min = -1;
    cutset.cuts[1].max = 200;
    cutset.cuts[1].cutvalue = &cTrigMuPtMin;
    cutset.cuts[1].ismin = true;

    cutset.cuts[2].name = Form("mass_PF > %6.2f", cDimuMassMin);
    cutset.cuts[2].bins = 150;
    cutset.cuts[2].min = 50;
    cutset.cuts[2].max = 200;
    cutset.cuts[2].cutvalue = &cDimuMassMin;
    cutset.cuts[2].ismin = true;

    cutset.cuts[3].name = "recoMu1 or recoMu2 matched to (hlt2 or hlt3)";
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
    float leadPt = -999;
    float subleadPt = -999; 
    float eta0 = -999;
    float eta1 = -999;
    float dimu_mass = -999;
    int charge0 = -999;
    int charge1 = -999;

    // the indices for the muons of the dimu pair
    int mu1 = vars.dimuCand->iMu1;
    int mu2 = vars.dimuCand->iMu2;

    if(useReco)
    {
        leadPt = TMath::Max(vars.recoMuons->at(mu1).pt, vars.recoMuons->at(mu2).pt);
        subleadPt = TMath::Min(vars.recoMuons->at(mu1).pt, vars.recoMuons->at(mu2).pt);
        eta0 = vars.recoMuons->at(mu1).eta;
        eta1 = vars.recoMuons->at(mu2).eta;
        dimu_mass = vars.dimuCand->mass_PF;
        charge0 = vars.recoMuons->at(mu1).charge;
        charge1 = vars.recoMuons->at(mu2).charge;
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
        if(!(vars.recoMuons->at(mu1).iso() <= cMaxRelIso) ) return false;
    }
    if(cutset.cuts[8].on)
    {
        if(!(vars.recoMuons->at(mu2).iso() <= cMaxRelIso) ) return false;
    }
    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString FEWZCompareCuts::string()
{
    return TString("FEWZ_Comparison_Cuts");
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

    cutset.cuts[0].name = Form("recoLead.pt > %6.2f", cLeadPtMin);
    cutset.cuts[0].bins = 200;
    cutset.cuts[0].min = 0;
    cutset.cuts[0].max = 200;
    cutset.cuts[0].cutvalue = &cLeadPtMin;
    cutset.cuts[0].ismin = true;

    cutset.cuts[1].name = Form("recoSublead.pt > %6.2f", cSubleadPtMin);
    cutset.cuts[1].bins = 200;
    cutset.cuts[1].min = 0;
    cutset.cuts[1].max = 200;
    cutset.cuts[1].cutvalue = &cSubleadPtMin;
    cutset.cuts[1].ismin = true;

    cutset.cuts[2].name = Form("recoMu1.eta < %6.3f", cMaxEta);
    cutset.cuts[2].bins = 50;
    cutset.cuts[2].min = 0;
    cutset.cuts[2].max = 3;
    cutset.cuts[2].cutvalue = &cMaxEta;
    cutset.cuts[2].ismin = false;

    cutset.cuts[3].name = Form("recoMu2.eta < %6.3f", cMaxEta);
    cutset.cuts[3].bins = 50;
    cutset.cuts[3].min = 0;
    cutset.cuts[3].max = 3;
    cutset.cuts[3].cutvalue = &cMaxEta;
    cutset.cuts[3].ismin = false;

    cutset.cuts[4].name = Form("dimuCand.mass_PF > %6.2f", cDimuMassMin);
    cutset.cuts[4].bins = 100;
    cutset.cuts[4].min = 110;
    cutset.cuts[4].max = 310;
    cutset.cuts[4].cutvalue = &cDimuMassMin;
    cutset.cuts[4].ismin = true;

    cutset.cuts[5].name = Form("dimuCand.mass_PF < %6.2f", cDimuMassMax);
    cutset.cuts[5].bins = 100;
    cutset.cuts[5].min = 110;
    cutset.cuts[5].max = 310;
    cutset.cuts[5].cutvalue = &cDimuMassMax;
    cutset.cuts[5].ismin = false;

    cutset.cuts[6].name = "recoMu1.charge != recoMu2.charge";
    cutset.cuts[6].bins = 2;
    cutset.cuts[6].min = 0;
    cutset.cuts[6].max = 2;
    cutset.cuts[6].ismin = true;

    cutset.cuts[7].name = Form("recoMu1.iso < %5.3f");
    cutset.cuts[7].bins = 50;
    cutset.cuts[7].min = 0;
    cutset.cuts[7].max = 1;
    cutset.cuts[7].cutvalue = &cMaxRelIso;
    cutset.cuts[7].ismin = false;

    cutset.cuts[8].name = Form("recoMu2.iso < %5.3f");
    cutset.cuts[8].bins = 50;
    cutset.cuts[8].min = 0;
    cutset.cuts[8].max = 1;
    cutset.cuts[8].cutvalue = &cMaxRelIso;
    cutset.cuts[8].ismin = false;
}

