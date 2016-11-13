///////////////////////////////////////////////////////////////////////////
//                         MuSelectionTools.cxx                          //
//=======================================================================//
//                                                                       //
//        Select valid muons beyond the candidate pair.                  //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "MuSelectionTools.h"
#include "TMath.h"
#include "TLorentzVector.h"

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// ___________________MuSelectionTools______________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

MuSelectionTools::MuSelectionTools()
{
    cMuSelectionPtMin = 10;
    cMuSelectionEtaMax = 2.4;
    cMuSelectionIsoMax = 0.12;
    cMuSelectionID = 0;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

MuSelectionTools::MuSelectionTools(float muSelectionPtMin, float muSelectionEtaMax, float muSelectionIsoMax, int muSelectionID)
{
    cMuSelectionPtMin = muSelectionPtMin;
    cMuSelectionEtaMax = muSelectionEtaMax;
    cMuSelectionIsoMax = muSelectionIsoMax;
    cMuSelectionID = muSelectionID;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void MuSelectionTools::getValidMus(VarSet& vars, std::vector<TLorentzVector>& muvec)
{
    for(unsigned int j=0; j < vars.recoMuons.nMuons && j < vars.recoMuons.arraySize; ++j)
    {
        bool id = false;
        if(cMuSelectionID == 0) id = vars.recoMuons.isTightMuon[j]; 
        if(cMuSelectionID == 1) id = vars.recoMuons.isMediumMuon[j]; 
        if(cMuSelectionID == 2) id = vars.recoMuons.isLooseMuon[j]; 

        // Pt, Eta and ID
        if(!(vars.recoMuons.pt[j] > cMuSelectionPtMin && TMath::Abs(vars.recoMuons.eta[j]) < cMuSelectionEtaMax && id)) 
            continue;

        // isolation
        if(!((vars.recoMuons.sumChargedHadronPtR03[j] + TMath::Max(0.0,vars.recoMuons.sumNeutralHadronEtR03[j]+vars.recoMuons.sumPhotonEtR03[j]
          - 0.5*vars.recoMuons.sumPUPtR03[j]))/vars.recoMuons.pt[j] <= cMuSelectionIsoMax))
            continue;

        // passes all selections, add to valid extra muons
        TLorentzVector mu4vec; 
        mu4vec.SetPtEtaPhiM(vars.recoMuons.pt[j],vars.recoMuons.eta[j],vars.recoMuons.phi[j],MASS_MUON);
        muvec.push_back(mu4vec);
    }
}
