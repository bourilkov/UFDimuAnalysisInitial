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
    cMuonSelectionIsoMax = 0.25;
    cMuonSelectionID = 1;
    cUseMedium2016 = false;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

MuonCollectionCleaner::MuonCollectionCleaner(float muSelectionPtMin, float muSelectionEtaMax, float muSelectionIsoMax, int muSelectionID, bool useMedium2016)
{
    cMuonSelectionPtMin = muSelectionPtMin;
    cMuonSelectionEtaMax = muSelectionEtaMax;
    cMuonSelectionIsoMax = muSelectionIsoMax;
    cMuonSelectionID = muSelectionID;
    cUseMedium2016 = useMedium2016;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void MuonCollectionCleaner::getValidMuons(VarSet& vars, std::vector<TLorentzVector>& muvec, bool exclude_pair)
{
    for(unsigned int j=0; j < vars.muons->size(); ++j)
    {
        if(exclude_pair && (j==vars.dimuCand->iMu1 || j==vars.dimuCand->iMu2)) continue;

        bool id = false;

        if(cMuonSelectionID == 0) id = vars.muons->at(j).isTightID; 
        if(cMuonSelectionID == 1 && cUseMedium2016)  id = vars.muons->at(j).isMediumID2016; 
        if(cMuonSelectionID == 1 && !cUseMedium2016) id = vars.muons->at(j).isMediumID; 
        if(cMuonSelectionID == 2) id = vars.muons->at(j).isLooseID; 

        // Pt, Eta and ID
        if(!(vars.muons->at(j).pt > cMuonSelectionPtMin && TMath::Abs(vars.muons->at(j).eta) < cMuonSelectionEtaMax && id)) 
            continue;

        // isolation
        if(!(vars.muons->at(j).iso() <= cMuonSelectionIsoMax))
            continue;

        // passes all selections, add to valid extra muons
        muvec.push_back(vars.muons->at(j).get4vec());
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void MuonCollectionCleaner::getValidMuons(VarSet& vars, std::vector<TLorentzVector>& muvec, std::vector<TLorentzVector>& xmuvec)
{
    for(unsigned int j=0; j < vars.muons->size(); ++j)
    {

        bool id = false;

        if(cMuonSelectionID == 0) id = vars.muons->at(j).isTightID; 
        if(cMuonSelectionID == 1 && cUseMedium2016)  id = vars.muons->at(j).isMediumID2016; 
        if(cMuonSelectionID == 1 && !cUseMedium2016) id = vars.muons->at(j).isMediumID; 
        if(cMuonSelectionID == 2) id = vars.muons->at(j).isLooseID; 

        // Pt, Eta and ID
        if(!(vars.muons->at(j).pt > cMuonSelectionPtMin && TMath::Abs(vars.muons->at(j).eta) < cMuonSelectionEtaMax && id)) 
            continue;

        // isolation
        if(!(vars.muons->at(j).iso() <= cMuonSelectionIsoMax))
            continue;

        // passes all selections, add to valid muons, valid extra muons
        muvec.push_back(vars.muons->at(j).get4vec());
        if(j!=vars.dimuCand->iMu1 && j!=vars.dimuCand->iMu2) xmuvec.push_back(vars.muons->at(j).get4vec());
    }
}
