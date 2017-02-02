//JetCollectionCleaner.h

#ifndef ADD_JETSELECTIONTOOLS
#define ADD_JETSELECTIONTOOLS

#include "VarSet.h"
#include <vector>
#include "CollectionCleaner.hxx"

class JetCollectionCleaner : public CollectionCleaner
{
    public:
        JetCollectionCleaner();
        JetCollectionCleaner(float cJetSelectionPtMin, float cJetSelectionEtaMax, float cJetSelectiondRMin, float cJetSelectionBTagMin, float cJetSelectionBJetEtaMax);

        float cJetSelectionPtMin; 
        float cJetSelectionEtaMax;       
        float cJetSelectiondRMin;
        float cJetSelectionBTagMin;
        float cJetSelectionBJetEtaMax;
 
        static float dR(float eta1, float phi1, float eta2, float phi2);

        void getValidJetsdR(VarSet& vars, std::vector<TLorentzVector>& jetvec);
        void getValidJets(VarSet& vars, std::vector<TLorentzVector>& jetvec);

        void getValidBJets(VarSet& vars, std::vector<TLorentzVector>& jetvec);
        void getValidBJetsdR(VarSet& vars, std::vector<TLorentzVector>& jetvec);

        void getValidGenJets(VarSet& vars, std::vector<TLorentzVector>& jetvec);

        bool jetID(VarSet& vars, unsigned int jet, int id);
};

#endif
