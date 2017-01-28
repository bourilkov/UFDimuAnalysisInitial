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
    for(unsigned int j=0; j < vars.recoElectrons.nElectrons && j < vars.recoElectrons.arraySize; ++j)
    {
        bool id = false;
        if(cElectronSelectionID == 0) id = vars.recoElectrons.isTightElectron[j]; 
        if(cElectronSelectionID == 1) id = vars.recoElectrons.isMediumElectron[j]; 
        if(cElectronSelectionID == 2) id = vars.recoElectrons.isLooseElectron[j]; 
        if(cElectronSelectionID == 3) id = vars.recoElectrons.isVetoElectron[j]; 

        // Pt, Eta and ID
        if(!(vars.recoElectrons.pt[j] > cElectronSelectionPtMin && TMath::Abs(vars.recoElectrons.eta[j]) < cElectronSelectionEtaMax && id)) 
            continue;

        // Conversion Veto
        if(!vars.recoElectrons.passConversionVeto[j])
            continue;

        // missing inner hits
        if(!(TMath::Abs(vars.recoElectrons.missingInnerHits[j]) <= 1))
            continue;

        // isolation
        if(!((vars.recoElectrons.sumChargedHadronPtR03[j] + TMath::Max(0.0,vars.recoElectrons.sumNeutralHadronEtR03[j]+vars.recoElectrons.sumPhotonEtR03[j]
          - 0.5*vars.recoElectrons.sumPUPtR03[j]))/vars.recoElectrons.pt[j] <= cElectronSelectionIsoMax))
            continue;

        // passes all selections, add to valid extra electrons
        TLorentzVector electron4vec; 
        electron4vec.SetPtEtaPhiM(vars.recoElectrons.pt[j],vars.recoElectrons.eta[j],vars.recoElectrons.phi[j],MASS_ELECTRON);
        electronvec.push_back(electron4vec);
    }
}
