//MuonCollectionCleaner.h

#ifndef ADD_MUSELECTIONTOOLS
#define ADD_MUSELECTIONTOOLS

#include "VarSet.h"
#include <vector>
#include "CollectionCleaner.hxx"

class MuonCollectionCleaner : public CollectionCleaner
{
    public:
        MuonCollectionCleaner();
        MuonCollectionCleaner(float cMuonSelectionPtMin, float cMuonSelectionEtaMax, float cMuonSelectionIsoMax, int cMuonSelectionID, bool useMedium2016);

        float cMuonSelectionPtMin; 
        float cMuonSelectionEtaMax;       
        float cMuonSelectionIsoMax;
        int   cMuonSelectionID;
        bool  cUseMedium2016;

        void getValidMuons(VarSet& vars, std::vector<TLorentzVector>& muvec, bool exclude_pair=false);
        void getValidMuons(VarSet& vars, std::vector<TLorentzVector>& muvec, std::vector<TLorentzVector>& xmuvec);
};

#endif
