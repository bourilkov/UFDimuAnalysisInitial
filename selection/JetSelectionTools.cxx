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
    cJetSelectiondRMax = 0.3;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

JetSelectionTools::JetSelectionTools(float jetSelectionPtMin, float jetSelectionEtaMax, float jetSelectiondRMax)
{
    cJetSelectionPtMin = jetSelectionPtMin;
    cJetSelectionEtaMax = jetSelectionEtaMax;
    cJetSelectiondRMax = jetSelectiondRMax;
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
        if(jets.pt[j] > cJetSelectionPtMin && TMath::Abs(jets.eta[j]) < cJetSelectionEtaMax)
            nValidJets++;
    }   
    return nValidJets;
}
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

float JetSelectionTools::dR(float eta1, float phi1, float eta2, float phi2) 
{
// Determine the number of valid jets using the given cuts

    // it's easiest to use the TLorentzVector class to account for the dPhi > pi problem
    TLorentzVector v1;
    TLorentzVector v2;
    v1.SetPtEtaPhiM(10, eta1, phi1, 0); // doesn't matter what pt and mass are we only care about the dR value
    v2.SetPtEtaPhiM(10, eta2, phi2, 0); // same thing
    return v1.DeltaR(v2);
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void JetSelectionTools::getValidJets(_PFJetInfo& jets, std::vector<TLorentzVector>& jetvec) 
{
// Determine the number of valid jets using the given cuts
// Use this when the jets have already been cut by dR in CMSSW

    for(unsigned int j=0; j < jets.nJets && j < 10; ++j)
    {   
        if(jets.pt[j] > cJetSelectionPtMin && TMath::Abs(jets.eta[j]) < cJetSelectionEtaMax)
        {
            TLorentzVector jet4vec; 
            jet4vec.SetPtEtaPhiM(jets.pt[j],jets.eta[j],jets.phi[j],jets.mass[j]);
            jetvec.push_back(jet4vec);
        }
    }   
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void JetSelectionTools::getValidJetsdR(VarSet& vars, std::vector<TLorentzVector>& jetvec)
{
// Determine the number of valid jets using the given cuts
// Cut jets by dR here instead of in CMSSW

    for(unsigned int j=0; j < vars.jets.nJets && j < 10; ++j)
    {
        // Pt and Eta selections
        if(vars.jets.pt[j] > cJetSelectionPtMin && TMath::Abs(vars.jets.eta[j]) < cJetSelectionEtaMax)
        {
            // dR vs muons selection
            if(!(dR(vars.jets.eta[j], vars.jets.phi[j], vars.reco1.eta, vars.reco1.phi) < cJetSelectiondRMax) && !(dR(vars.jets.eta[j], vars.jets.phi[j], vars.reco2.eta, vars.reco2.phi) < cJetSelectiondRMax))
            {
                TLorentzVector jet4vec; 
                jet4vec.SetPtEtaPhiM(vars.jets.pt[j],vars.jets.eta[j],vars.jets.phi[j],vars.jets.mass[j]);
                // passes all selections, add to valid jets
                jetvec.push_back(jet4vec);
            }
        }
    }
}

