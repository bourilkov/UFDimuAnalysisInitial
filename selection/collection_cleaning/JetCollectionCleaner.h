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

        void getValidJets(VarSet& vars, std::vector<TLorentzVector>& jetvec, bool require_b = false);
        void getValidJets(VarSet& vars, std::vector<TLorentzVector>& jetvec, std::vector<TLorentzVector>& bjetvec, bool print=false);

        bool jetID(VarSet& vars, unsigned int jet, int id);
};

#endif
