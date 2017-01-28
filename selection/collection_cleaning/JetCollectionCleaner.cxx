///////////////////////////////////////////////////////////////////////////
//                             JetCollectionCleaner.cxx                     //
//=======================================================================//
//                                                                       //
//        Work with the _PFJetsInfo data structure.                      //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "JetCollectionCleaner.h"
#include "TMath.h"
#include "TLorentzVector.h"

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________JetCollectionCleaner______________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

JetCollectionCleaner::JetCollectionCleaner()
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

JetCollectionCleaner::JetCollectionCleaner(float jetSelectionPtMin, float jetSelectionEtaMax, float jetSelectiondRMin, float jetSelectionBTagMin, float jetSelectionBJetEtaMax)
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

float JetCollectionCleaner::dR(float eta1, float phi1, float eta2, float phi2) 
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

void JetCollectionCleaner::getValidJetsdR(VarSet& vars, std::vector<TLorentzVector>& jetvec)
{
// Determine the number of valid jets using the given cuts
// Cut jets by dR here instead of in CMSSW
    for(unsigned int j=0; j < vars.jets->size(); ++j)
    {
        // Pt and Eta selections
        if(vars.jets->at(j).pt > cJetSelectionPtMin && TMath::Abs(vars.jets->at(j).eta) < cJetSelectionEtaMax)
        {
            // dR vs muons selection
            if(!(dR(vars.jets->at(j).eta, vars.jets->at(j).phi, vars.recoMuons->at(0).eta, vars.recoMuons->at(0).phi) < cJetSelectiondRMin) 
                 && !(dR(vars.jets->at(j).eta, vars.jets->at(j).phi, vars.recoMuons->at(1).eta, vars.recoMuons->at(1).phi) < cJetSelectiondRMin))
            {
                TLorentzVector jet4vec; 
                jet4vec.SetPtEtaPhiM(vars.jets->at(j).pt,vars.jets->at(j).eta,vars.jets->at(j).phi,vars.jets->at(j).mass);
                // passes all selections, add to valid jets
                jetvec.push_back(jet4vec);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void JetCollectionCleaner::getValidBJetsdR(VarSet& vars, std::vector<TLorentzVector>& jetvec)
{
// Determine the number of valid jets using the given cuts
// Cut jets by dR here instead of in CMSSW
    for(unsigned int j=0; j < vars.jets->size(); ++j)
    {
        // Pt selection, btag-id, eta selections
        if(vars.jets->at(j).pt > cJetSelectionPtMin && vars.jets->at(j).CSV > cJetSelectionBTagMin && TMath::Abs(vars.jets->at(j).eta) < cJetSelectionBJetEtaMax)
        {
            // dR vs muons selection
            if(!(dR(vars.jets->at(j).eta, vars.jets->at(j).phi, vars.recoMuons->at(0).eta, vars.recoMuons->at(0).phi) < cJetSelectiondRMin) 
                 && !(dR(vars.jets->at(j).eta, vars.jets->at(j).phi, vars.recoMuons->at(1).eta, vars.recoMuons->at(1).phi) < cJetSelectiondRMin))
            {
                TLorentzVector jet4vec; 
                jet4vec.SetPtEtaPhiM(vars.jets->at(j).pt,vars.jets->at(j).eta,vars.jets->at(j).phi,vars.jets->at(j).mass);
                // add to valid b-jets
                jetvec.push_back(jet4vec);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void JetCollectionCleaner::getValidJets(VarSet& vars, std::vector<TLorentzVector>& jetvec)
{
// Determine the number of valid jets using the given cuts
    for(unsigned int j=0; j < vars.jets->size(); ++j)
    {
        // Pt and Eta selections
        if(vars.jets->at(j).pt > cJetSelectionPtMin && TMath::Abs(vars.jets->at(j).eta) < cJetSelectionEtaMax)
        {
           TLorentzVector jet4vec; 
           jet4vec.SetPtEtaPhiM(vars.jets->at(j).pt,vars.jets->at(j).eta,vars.jets->at(j).phi,vars.jets->at(j).mass);
           // passes all selections, add to valid jets
           jetvec.push_back(jet4vec);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void JetCollectionCleaner::getValidBJets(VarSet& vars, std::vector<TLorentzVector>& jetvec)
{
// Determine the number of valid jets using the given cuts
    for(unsigned int j=0; j < vars.jets->size(); ++j)
    {
        // Pt selection
        if(vars.jets->at(j).pt > cJetSelectionPtMin && vars.jets->at(j).CSV > cJetSelectionBTagMin 
           && TMath::Abs(vars.jets->at(j).eta) < cJetSelectionBJetEtaMax)
        {
           TLorentzVector jet4vec; 
           jet4vec.SetPtEtaPhiM(vars.jets->at(j).pt,vars.jets->at(j).eta,vars.jets->at(j).phi,vars.jets->at(j).mass);
           // add to valid b-jets
           jetvec.push_back(jet4vec);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void JetCollectionCleaner::getValidGenJets(VarSet& vars, std::vector<TLorentzVector>& jetvec)
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
/*
bool JetCollectionCleaner::jetID(VarSet& vars, unsigned int jet, int id)
{
    bool looseJetID, tightJetID, tightLepVetoJetID;
    tightLepVetoJetID = false;

    float eta = TMath::Abs(vars.jets->at(jet).eta);

    int NHF = vars.jets->at(jet).nhf;
    int NEMF = vars.jets->at(jet).nef;
    int NumConst = vars.jets->at(jet).cm+vars.jets->at(jet).nm; // not available in awb ttree
    int MUF = vars.jets->at(jet).muf;
    int CHF = vars.jets->at(jet).chf;
    int CHM = vars.jets->at(jet).cm;
    int CEMF = vars.jets->at(jet).cef;
    int NumNeutralParticle = vars.jets->at(jet).nm;            // same nm not available in awb ttree

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
*/
