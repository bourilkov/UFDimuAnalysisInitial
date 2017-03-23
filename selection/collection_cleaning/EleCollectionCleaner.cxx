///////////////////////////////////////////////////////////////////////////
//                         EleCollectionCleaner.cxx                      //
//=======================================================================//
//                                                                       //
//        Select valid electrons.                                        //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "EleCollectionCleaner.h"
#include "TMath.h"
#include "TLorentzVector.h"

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// ___________________EleCollectionCleaner_______________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

EleCollectionCleaner::EleCollectionCleaner()
{
    cElectronSelectionPtMin = 10;
    cElectronSelectionEtaMax = 2.5;
    cElectronSelectionIsoMax = 0.15;
    cElectronSelectionID = 1;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

EleCollectionCleaner::EleCollectionCleaner(float electronSelectionPtMin, float electronSelectionEtaMax, float electronSelectionIsoMax, int electronSelectionID)
{
    cElectronSelectionPtMin = electronSelectionPtMin;
    cElectronSelectionEtaMax = electronSelectionEtaMax;
    cElectronSelectionIsoMax = electronSelectionIsoMax;
    cElectronSelectionID = electronSelectionID;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void EleCollectionCleaner::getValidElectrons(VarSet& vars, std::vector<TLorentzVector>& electronvec)
{
    for(unsigned int j=0; j < vars.electrons->size(); ++j)
    {
        bool id = false;
        if(cElectronSelectionID == 0) id = vars.electrons->at(j).isTightID; 
        if(cElectronSelectionID == 1) id = vars.electrons->at(j).isMediumID; 
        if(cElectronSelectionID == 2) id = vars.electrons->at(j).isLooseID; 
        if(cElectronSelectionID == 3) id = vars.electrons->at(j).isVetoID; 

        // crack in the hcal
        double eta = TMath::Abs(vars.electrons->at(j).eta);
        if(!(eta < 1.4442 || (eta >  1.566 && eta < cElectronSelectionEtaMax)) )
            continue;

        // Pt, Eta and ID
        if(!(vars.electrons->at(j).pt > cElectronSelectionPtMin && TMath::Abs(vars.electrons->at(j).eta) < cElectronSelectionEtaMax && id)) 
            continue;

        // Conversion Veto
        if(!vars.electrons->at(j).passConversionVeto)
            continue;

        // missing inner hits
        if(!(TMath::Abs(vars.electrons->at(j).missingInnerHits) <= 1))
            continue;

        // isolation
        if(!(vars.electrons->at(j).iso() <= cElectronSelectionIsoMax))
            continue;

        // passes all selections, add to valid extra electrons
        electronvec.push_back(vars.electrons->at(j).get4vec());
    }
}
