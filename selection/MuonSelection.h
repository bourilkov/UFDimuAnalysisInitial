///////////////////////////////////////////////////////////////////////////
//                             MuonSelection.h                           //
//=======================================================================//
//                                                                       //
//        Eliminate events by cutting on individual muon info.           //
//        The dimuon candidate muons must pass these criteria.           //
//        Usually cuts based on muon eta, pt, and iso.                   //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ADD_MUONSELECTIONCUTS
#define ADD_MUONSELECTIONCUTS

#include "Cut.h"

class Run2MuonSelectionCuts : public Cut
{
    public:
        Run2MuonSelectionCuts();
        Run2MuonSelectionCuts(float minPt, float maxEta, float maxRelIso);

        float cMinPt;               // >
        float cMaxEta;              // <
        float cMaxRelIso;           // <

        void makeCutSet();
        bool evaluate(VarSet& vars);
        TString string();
};

#endif
