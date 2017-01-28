///////////////////////////////////////////////////////////////////////////
//                         TauCollectionCleaner.cxx                    //
//=======================================================================//
//                                                                       //
//        Select valid taus beyond the candidate pair.            //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "TauCollectionCleaner.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// ___________________TauCollectionCleaner_____________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

TauCollectionCleaner::TauCollectionCleaner()
{
    cTauSelectionPtMin = 10;
    cTauSelectionEtaMax = 2.5;
    cTauSelectionIDs = std::vector<unsigned int>();
    cTauSelectionIDs = {0, 8, 10, 11}; // ids that need to be checked
                                       // decayModeFinding, byVTightIsolationMVArun2v1DBoldDMwLT, againstMuonTight3, againstElectronVLooseMVA6
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TauCollectionCleaner::TauCollectionCleaner(float tauSelectionPtMin, float tauSelectionEtaMax, std::vector<unsigned int>& tauSelectionIDs)
{
    cTauSelectionPtMin = tauSelectionPtMin;
    cTauSelectionEtaMax = tauSelectionEtaMax;
    cTauSelectionIDs = tauSelectionIDs;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void TauCollectionCleaner::getValidTaus(VarSet& vars, std::vector<TLorentzVector>& tauvec)
{
    for(unsigned int j=0; j < vars.recoTaus.nTaus && j < vars.recoTaus.arraySize; ++j)
    {
        // Pt, Eta
        if(!(vars.recoTaus.pt[j] > cTauSelectionPtMin && TMath::Abs(vars.recoTaus.eta[j]) < cTauSelectionEtaMax)) 
            continue;

        // IDs, j is the tau, i is the id
        for(unsigned int i=0; i<cTauSelectionIDs.size(); i++)
        {
            if(!(cTauSelectionIDs[i]<vars.recoTaus.idArraySize)) std::cout << "ERROR: TAU SELECTION ID OUT OF BOUNDS in TauCollectionCleaner::getValidTaus." << std::endl;
            if(!(vars.recoTaus.tauID[j][cTauSelectionIDs[i]] > 0.5)) continue;
        }

        // passes all selections, add to valid extra taus
        TLorentzVector tau4vec; 
        tau4vec.SetPtEtaPhiM(vars.recoTaus.pt[j],vars.recoTaus.eta[j],vars.recoTaus.phi[j],MASS_TAU);
        tauvec.push_back(tau4vec);
    }
}
