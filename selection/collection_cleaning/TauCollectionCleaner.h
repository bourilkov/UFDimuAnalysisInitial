//TauCollectionCleaner.h

#ifndef ADD_TAUSELECTIONTOOLS
#define ADD_TAUSELECTIONTOOLS

#include "VarSet.h"
#include <vector>

class TauCollectionCleaner
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
