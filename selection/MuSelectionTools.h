//MuSelectionTools.h

#ifndef ADD_MUSELECTIONTOOLS
#define ADD_MUSELECTIONTOOLS

#include "VarSet.h"
#include <vector>

class MuSelectionTools
{
    public:
        MuSelectionTools();
        MuSelectionTools(float cMuSelectionPtMin, float cMuSelectionEtaMax, float cMuSelectionIsoMax, int cMuSelectionID);

        float cMuSelectionPtMin; 
        float cMuSelectionEtaMax;       
        float cMuSelectionIsoMax;
        int   cMuSelectionID;

        void getValidMus(VarSet& vars, std::vector<TLorentzVector>& muvec);
};

#endif
