//TauCollectionCleaner.h

#ifndef ADD_TAUSELECTIONTOOLS
#define ADD_TAUSELECTIONTOOLS

#include "VarSet.h"
#include <vector>
#include "CollectionCleaner.hxx"

class TauCollectionCleaner : public CollectionCleaner
{
    public:
        TauCollectionCleaner();
        TauCollectionCleaner(float cTauSelectionPtMin, float cTauSelectionEtaMax, std::vector<unsigned int>& cTauSelectionIDs);

        float cTauSelectionPtMin; 
        float cTauSelectionEtaMax;       

        std::vector<unsigned int> cTauSelectionIDs;

        void getValidTaus(VarSet& vars, std::vector<TLorentzVector>& tauvec);
};

#endif
