//MuonCollectionCleaner.h

#ifndef ADD_MUSELECTIONTOOLS
#define ADD_MUSELECTIONTOOLS

#include "VarSet.h"
#include <vector>

class MuonCollectionCleaner
{
    public:
        MuonCollectionCleaner();
        MuonCollectionCleaner(float cMuonSelectionPtMin, float cMuonSelectionEtaMax, float cMuonSelectionIsoMax, int cMuonSelectionID);

        float cMuonSelectionPtMin; 
        float cMuonSelectionEtaMax;       
        float cMuonSelectionIsoMax;
        int   cMuonSelectionID;

        void getValidMuons(VarSet& vars, std::vector<TLorentzVector>& muvec);
};

#endif
