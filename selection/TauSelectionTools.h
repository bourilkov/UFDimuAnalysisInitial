//TauSelectionTools.h

#ifndef ADD_TAUSELECTIONTOOLS
#define ADD_TAUSELECTIONTOOLS

#include "VarSet.h"
#include <vector>

class TauSelectionTools
{
    public:
        TauSelectionTools();
        TauSelectionTools(float cTauSelectionPtMin, float cTauSelectionEtaMax, std::vector<unsigned int>& cTauSelectionIDs);

        float cTauSelectionPtMin; 
        float cTauSelectionEtaMax;       

        std::vector<unsigned int> cTauSelectionIDs;

        void getValidTaus(VarSet& vars, std::vector<TLorentzVector>& tauvec);
};

#endif
