//MuonSelectionTools.h

#ifndef ADD_MUSELECTIONTOOLS
#define ADD_MUSELECTIONTOOLS

#include "VarSet.h"
#include <vector>

class MuonSelectionTools
{
    public:
        MuonSelectionTools();
        MuonSelectionTools(float cMuonSelectionPtMin, float cMuonSelectionEtaMax, float cMuonSelectionIsoMax, int cMuonSelectionID);

        float cMuonSelectionPtMin; 
        float cMuonSelectionEtaMax;       
        float cMuonSelectionIsoMax;
        int   cMuonSelectionID;

        void getValidMuons(VarSet& vars, std::vector<TLorentzVector>& muvec);
};

#endif
