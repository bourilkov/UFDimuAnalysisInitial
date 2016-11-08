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
    cJetSelectiondRMin = 0.3;
    cJetSelectionBTagMin = 0.8;
    cJetSelectionBJetEtaMax = 2.4;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

JetSelectionTools::JetSelectionTools(float jetSelectionPtMin, float jetSelectionEtaMax, float jetSelectiondRMin, float jetSelectionBTagMin, float jetSelectionBJetEtaMax)
{
    cJetSelectionPtMin = jetSelectionPtMin;
    cJetSelectionEtaMax = jetSelectionEtaMax;
    cJetSelectiondRMin = jetSelectiondRMin;
    cJetSelectionBTagMin = jetSelectionBTagMin;
    cJetSelectionBJetEtaMax = jetSelectionBJetEtaMax;
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

void JetSelectionTools::getValidJetsdR(VarSet& vars, std::vector<TLorentzVector>& jetvec)
{
// Determine the number of valid jets using the given cuts
// Cut jets by dR here instead of in CMSSW
    for(unsigned int j=0; j < vars.jets.nJets && j < vars.jets.arraySize; ++j)
    {
        // Pt and Eta selections
        if(vars.jets.pt[j] > cJetSelectionPtMin && TMath::Abs(vars.jets.eta[j]) < cJetSelectionEtaMax)
        {
            // dR vs muons selection
            if(!(dR(vars.jets.eta[j], vars.jets.phi[j], vars.recoMuons.eta[0], vars.recoMuons.phi[0]) < cJetSelectiondRMin) 
                 && !(dR(vars.jets.eta[j], vars.jets.phi[j], vars.recoMuons.eta[1], vars.recoMuons.phi[1]) < cJetSelectiondRMin))
            {
                TLorentzVector jet4vec; 
                jet4vec.SetPtEtaPhiM(vars.jets.pt[j],vars.jets.eta[j],vars.jets.phi[j],vars.jets.mass[j]);
                // passes all selections, add to valid jets
                jetvec.push_back(jet4vec);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void JetSelectionTools::getValidBJetsdR(VarSet& vars, std::vector<TLorentzVector>& jetvec)
{
// Determine the number of valid jets using the given cuts
// Cut jets by dR here instead of in CMSSW
    for(unsigned int j=0; j < vars.jets.nJets && j < vars.jets.arraySize; ++j)
    {
        // Pt selection, btag-id, eta selections
        if(vars.jets.pt[j] > cJetSelectionPtMin && vars.jets.isB[j] > cJetSelectionBTagMin && TMath::Abs(vars.jets.eta[j]) < cJetSelectionBJetEtaMax)
        {
            // dR vs muons selection
            if(!(dR(vars.jets.eta[j], vars.jets.phi[j], vars.recoMuons.eta[0], vars.recoMuons.phi[0]) < cJetSelectiondRMin) 
                 && !(dR(vars.jets.eta[j], vars.jets.phi[j], vars.recoMuons.eta[1], vars.recoMuons.phi[1]) < cJetSelectiondRMin))
            {
                TLorentzVector jet4vec; 
                jet4vec.SetPtEtaPhiM(vars.jets.pt[j],vars.jets.eta[j],vars.jets.phi[j],vars.jets.mass[j]);
                // add to valid b-jets
                jetvec.push_back(jet4vec);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void JetSelectionTools::getValidJets(VarSet& vars, std::vector<TLorentzVector>& jetvec)
{
// Determine the number of valid jets using the given cuts
    for(unsigned int j=0; j < vars.jets.nJets && j < vars.jets.arraySize; ++j)
    {
        // Pt and Eta selections
        if(vars.jets.pt[j] > cJetSelectionPtMin && TMath::Abs(vars.jets.eta[j]) < cJetSelectionEtaMax)
        {
           TLorentzVector jet4vec; 
           jet4vec.SetPtEtaPhiM(vars.jets.pt[j],vars.jets.eta[j],vars.jets.phi[j],vars.jets.mass[j]);
           // passes all selections, add to valid jets
           jetvec.push_back(jet4vec);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void JetSelectionTools::getValidBJets(VarSet& vars, std::vector<TLorentzVector>& jetvec)
{
// Determine the number of valid jets using the given cuts
    for(unsigned int j=0; j < vars.jets.nJets && j < vars.jets.arraySize; ++j)
    {
        // Pt selection
        if(vars.jets.pt[j] > cJetSelectionPtMin && vars.jets.isB[j] > cJetSelectionBTagMin && TMath::Abs(vars.jets.eta[j]) < cJetSelectionBJetEtaMax)
        {
           TLorentzVector jet4vec; 
           jet4vec.SetPtEtaPhiM(vars.jets.pt[j],vars.jets.eta[j],vars.jets.phi[j],vars.jets.mass[j]);
           // add to valid b-jets
           jetvec.push_back(jet4vec);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void JetSelectionTools::getValidGenJets(VarSet& vars, std::vector<TLorentzVector>& jetvec)
{
// Determine the number of valid jets using the given cuts
// Cut jets by dR here instead of in CMSSW
    for(unsigned int j=0; j < vars.genJets.nJets && j < vars.genJets.arraySize; ++j)
    {
        // Pt and Eta selections
        if(vars.genJets.pt[j] > cJetSelectionPtMin && TMath::Abs(vars.genJets.eta[j]) < cJetSelectionEtaMax)
        {
           TLorentzVector jet4vec; 
           jet4vec.SetPtEtaPhiM(vars.genJets.pt[j],vars.genJets.eta[j],vars.genJets.phi[j],vars.genJets.mass[j]);
           // passes all selections, add to valid jets
           jetvec.push_back(jet4vec);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool JetSelectionTools::jetID(VarSet& vars, unsigned int jet, int id)
{
    bool looseJetID, tightJetID, tightLepVetoJetID;
    tightLepVetoJetID = false;

    float eta = TMath::Abs(vars.jets.eta[jet]);

    float NHF = vars.jets.nhf[jet];
    float NEMF = vars.jets.nef[jet];
    float NumConst = 2;
    // NumConst = vars.jets.cm[jet]+vars.jets.nm[jet];
    float MUF = vars.jets.muf[jet];
    float CHF = vars.jets.chf[jet];
    float CHM = vars.jets.cm[jet];
    float CEMF = vars.jets.cef[jet];
    float NumNeutralParticle = 11;
    // NumNeutralParticle = vars.jets.nm[jet];

    if(eta <= 2.7)
    {
        looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((eta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || eta>2.4) && eta<=2.7;
        tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((eta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || eta>2.4) && eta<=2.7;
        tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((eta<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || eta>2.4) && eta<=2.7;
    }
    if(eta > 2.7 && eta <= 3.0)
    {
        looseJetID = (NEMF<0.90 && NumNeutralParticle>2 && eta>2.7 && eta<=3.0 );
        tightJetID = (NEMF<0.90 && NumNeutralParticle>2 && eta>2.7 && eta<=3.0 );
    }
    if(eta > 3.0)
    {
        looseJetID = (NEMF<0.90 && NumNeutralParticle>10 && eta>3.0 );
        tightJetID = (NEMF<0.90 && NumNeutralParticle>10 && eta>3.0 );
    }

    if(id == 0) return tightLepVetoJetID;
    else if(id == 1) return tightJetID;
    else return looseJetID;
}
