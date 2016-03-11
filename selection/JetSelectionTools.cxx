///////////////////////////////////////////////////////////////////////////
//                             JetSelectionTools.cxx                     //
//=======================================================================//
//                                                                       //
//        Work with the _PFJetsInfo data structure.                      //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "JetSelectionTools.h"
#include "TMath.h"
#include "TLorentzVector.h"

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________JetSelectionTools______________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

JetSelectionTools::JetSelectionTools()
{
    cJetSelectionPtMin = 30;
    cJetSelectionEtaMax = 4.7;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

JetSelectionTools::JetSelectionTools(float jetSelectionPtMin, float jetSelectionEtaMax)
{
    cJetSelectionPtMin = jetSelectionPtMin;
    cJetSelectionEtaMax = jetSelectionEtaMax;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

int JetSelectionTools::getNValidJets(_PFJetInfo& jets) 
{
// Determine the number of valid jets using the given cuts

    int nValidJets = 0;
    for(unsigned int j=0; j < jets.nJets && j < 10; ++j)
    {   
        if(!(jets.pt[j] > cJetSelectionPtMin && TMath::Abs(jets.eta[j]) < cJetSelectionEtaMax))
            nValidJets++;
    }   
    return nValidJets;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void JetSelectionTools::getValidJets(_PFJetInfo& jets, std::vector<TLorentzVector>& jetvec) 
{
// Determine the number of valid jets using the given cuts

    for(unsigned int j=0; j < jets.nJets && j < 10; ++j)
    {   
        if(!(jets.pt[j] > cJetSelectionPtMin && TMath::Abs(jets.eta[j]) < cJetSelectionEtaMax))
        {
            TLorentzVector jet4vec; 
            jet4vec.SetPtEtaPhiM(jets.pt[j],jets.eta[j],jets.phi[j],jets.mass[j]);
            jetvec.push_back(jet4vec);
        }
    }   
}
