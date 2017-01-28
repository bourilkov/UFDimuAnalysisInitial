//MuonSelection.h

#ifndef ADD_MUONSELECTIONCUTS
#define ADD_MUONSELECTIONCUTS

#include "Cut.h"

class SynchMuonSelectionCuts : public Cut
{
    public:
        SynchMuonSelectionCuts();
        SynchMuonSelectionCuts(float minPt, float maxEta, float maxRelIso);

        float cMinPt;               // >
        float cMaxEta;              // <
        float cMaxRelIso;           // <

        void makeCutSet();
        bool evaluate(VarSet& vars);
        bool evaluate(_MuonInfo& recoMu, int m);
        TString string();
};

class Run1MuonSelectionCuts : public Cut
{
    public:
        Run1MuonSelectionCuts();
        Run1MuonSelectionCuts(float minPt, float maxEta, float maxRelIso);

        float cMinPt;               // >
        float cMaxEta;              // <
        float cMaxRelIso;           // <

        void makeCutSet();
        bool evaluate(VarSet& vars);
        bool evaluate(_MuonInfo& recoMu, int m);
        TString string();
};

#endif
