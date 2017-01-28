///////////////////////////////////////////////////////////////////////////
//                         EleCollectionCleaner.cxx                    //
//=======================================================================//
//                                                                       //
//        Select valid electrons beyond the candidate pair.            //
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
// ___________________EleCollectionCleaner_____________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

EleCollectionCleaner::EleCollectionCleaner()
{
    cElectronSelectionPtMin = 10;
    cElectronSelectionEtaMax = 2.5;
    cElectronSelectionIsoMax = 0.15;
    cElectronSelectionID = 0;
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
    for(unsigned int j=0; j < vars.recoElectrons->size(); ++j)
    {
        bool id = false;
        if(cElectronSelectionID == 0) id = vars.recoElectrons->at(j).isTightID; 
        if(cElectronSelectionID == 1) id = vars.recoElectrons->at(j).isMediumID; 
        if(cElectronSelectionID == 2) id = vars.recoElectrons->at(j).isLooseID; 
        if(cElectronSelectionID == 3) id = vars.recoElectrons->at(j).isVetoID; 

        // Pt, Eta and ID
        if(!(vars.recoElectrons->at(j).pt > cElectronSelectionPtMin && TMath::Abs(vars.recoElectrons->at(j).eta) < cElectronSelectionEtaMax && id)) 
            continue;

        // Conversion Veto
        if(!vars.recoElectrons->at(j).passConversionVeto)
            continue;

        // missing inner hits
        if(!(TMath::Abs(vars.recoElectrons->at(j).missingInnerHits) <= 1))
            continue;

        // isolation
        if(!((vars.recoElectrons->at(j).sumChargedHadronPtR03+ TMath::Max(0.0,vars.recoElectrons->at(j).sumNeutralHadronEtR03+
              vars.recoElectrons->at(j).sumPhotonEtR03 - 0.5*vars.recoElectrons->at(j).sumPUPtR03))/vars.recoElectrons->at(j).pt 
               <= cElectronSelectionIsoMax))
            continue;

        // passes all selections, add to valid extra electrons
        TLorentzVector electron4vec; 
        electron4vec.SetPtEtaPhiM(vars.recoElectrons->at(j).pt,vars.recoElectrons->at(j).eta,vars.recoElectrons->at(j).phi,MASS_ELECTRON);
        electronvec.push_back(electron4vec);
    }
}
