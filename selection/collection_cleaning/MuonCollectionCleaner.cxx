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
    for(unsigned int j=0; j < vars.recoMuons->size(); ++j)
    {
        bool id = false;
        if(cMuonSelectionID == 0) id = vars.recoMuons->at(j).isTightID; 
        if(cMuonSelectionID == 1) id = vars.recoMuons->at(j).isMediumID; 
        if(cMuonSelectionID == 2) id = vars.recoMuons->at(j).isLooseID; 

        // Pt, Eta and ID
        if(!(vars.recoMuons->at(j).pt > cMuonSelectionPtMin && TMath::Abs(vars.recoMuons->at(j).eta) < cMuonSelectionEtaMax && id)) 
            continue;

        // isolation
        if(!(vars.recoMuons->at(j).iso() <= cMuonSelectionIsoMax))
            continue;

        // passes all selections, add to valid extra muons
        muvec.push_back(vars.recoMuons->at(j).get4vec());
    }
}
