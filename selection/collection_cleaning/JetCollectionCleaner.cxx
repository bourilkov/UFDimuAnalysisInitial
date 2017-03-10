///////////////////////////////////////////////////////////////////////////
//                             JetCollectionCleaner.cxx                  //
//=======================================================================//
//                                                                       //
//        Select good jets from the main collection.                     //
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
// _______________________JetCollectionCleaner___________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

JetCollectionCleaner::JetCollectionCleaner()
{
    cJetSelectionPtMin = 30;
    cJetSelectionEtaMax = 4.7;
    cJetSelectiondRMin = 0.3;
    cJetSelectionBTagMin = 0.8484;
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

void JetCollectionCleaner::getValidJets(VarSet& vars, std::vector<TLorentzVector>& jetvec, std::vector<TLorentzVector>& bjetvec, bool print)
{
// Determine the number of valid jets using the given cuts
    for(unsigned int j=0; j < vars.jets->size(); ++j)
    {
        if(print) std::cout << Form("Checking > %s\n", vars.jets->at(j).outputInfo().Data());
        // Pt and Eta selections for a regular jet
        if(vars.jets->at(j).pt > cJetSelectionPtMin && TMath::Abs(vars.jets->at(j).eta) < cJetSelectionEtaMax)
        {
           if(print) std::cout << Form("Adding to jets > %s\n", vars.jets->at(j).outputInfo().Data());
           TLorentzVector jet4vec = vars.jets->at(j).get4vec();
           jetvec.push_back(jet4vec);

           // further selections for a bjet, eta should be tighter since we need the tracker
           if(vars.jets->at(j).CSV > cJetSelectionBTagMin && TMath::Abs(vars.jets->at(j).eta) < cJetSelectionBJetEtaMax)
           {
               if(print) std::cout << Form("Adding to bjets > %s\n", vars.jets->at(j).outputInfo().Data());
               bjetvec.push_back(jet4vec);
           }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void JetCollectionCleaner::getValidJets(VarSet& vars, std::vector<TLorentzVector>& jetvec, bool require_b)
{
// Determine the number of valid jets using the given cuts
    for(unsigned int j=0; j < vars.jets->size(); ++j)
    {
        // bjet selection
        if(require_b && vars.jets->at(j).pt > cJetSelectionPtMin && vars.jets->at(j).CSV > cJetSelectionBTagMin 
           && TMath::Abs(vars.jets->at(j).eta) < cJetSelectionBJetEtaMax)
        {
           TLorentzVector jet4vec = vars.jets->at(j).get4vec();
           jetvec.push_back(jet4vec);
        }
        // regular jet selection
        if(!require_b && vars.jets->at(j).pt > cJetSelectionPtMin && TMath::Abs(vars.jets->at(j).eta) < cJetSelectionEtaMax)
        {
           TLorentzVector jet4vec = vars.jets->at(j).get4vec();
           jetvec.push_back(jet4vec);
        }
    }
}
