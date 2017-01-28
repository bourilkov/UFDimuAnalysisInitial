///////////////////////////////////////////////////////////////////////////
//                         MuonCollectionCleaner.cxx                     //
//=======================================================================//
//                                                                       //
//        Select valid muons beyond the candidate pair.                  //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "MuonCollectionCleaner.h"
#include "TMath.h"
#include "TLorentzVector.h"

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// ___________________MuonCollectionCleaner______________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

MuonCollectionCleaner::MuonCollectionCleaner()
{
    cMuonSelectionPtMin = 10;
    cMuonSelectionEtaMax = 2.4;
    cMuonSelectionIsoMax = 0.12;
    cMuonSelectionID = 0;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

MuonCollectionCleaner::MuonCollectionCleaner(float muSelectionPtMin, float muSelectionEtaMax, float muSelectionIsoMax, int muSelectionID)
{
    cMuonSelectionPtMin = muSelectionPtMin;
    cMuonSelectionEtaMax = muSelectionEtaMax;
    cMuonSelectionIsoMax = muSelectionIsoMax;
    cMuonSelectionID = muSelectionID;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void MuonCollectionCleaner::getValidMuons(VarSet& vars, std::vector<TLorentzVector>& muvec)
{
    for(unsigned int j=0; j < vars.recoMuons.nMuons && j < vars.recoMuons.arraySize; ++j)
    {
        bool id = false;
        if(cMuonSelectionID == 0) id = vars.recoMuons.isTightMuon[j]; 
        if(cMuonSelectionID == 1) id = vars.recoMuons.isMediumMuon[j]; 
        if(cMuonSelectionID == 2) id = vars.recoMuons.isLooseMuon[j]; 

        // Pt, Eta and ID
        if(!(vars.recoMuons.pt[j] > cMuonSelectionPtMin && TMath::Abs(vars.recoMuons.eta[j]) < cMuonSelectionEtaMax && id)) 
            continue;

        // isolation
        if(!((vars.recoMuons.sumChargedHadronPtR03[j] + TMath::Max(0.0,vars.recoMuons.sumNeutralHadronEtR03[j]+vars.recoMuons.sumPhotonEtR03[j]
          - 0.5*vars.recoMuons.sumPUPtR03[j]))/vars.recoMuons.pt[j] <= cMuonSelectionIsoMax))
            continue;

        // passes all selections, add to valid extra muons
        TLorentzVector mu4vec; 
        mu4vec.SetPtEtaPhiM(vars.recoMuons.pt[j],vars.recoMuons.eta[j],vars.recoMuons.phi[j],MASS_MUON);
        muvec.push_back(mu4vec);
    }
}
