//SelectionCuts.h

#ifndef ADD_CUTS
#define ADD_CUTS

#include "Cut.h"

// Define the different cuts
class TightMuonIdCuts : public Cut
{
    public:

        TightMuonIdCuts();
        TightMuonIdCuts(int numTrackerLayers, int numValidMuonHits, int numValidPixelHits, int numOfMatchedStations, 
                        int normChiSquare, float d0_pV, float dz_PV);

        int cNumTrackerLayers;     // >
        int cNumValidMuonHits;     // >
        int cNumValidPixelHits;    // >
        int cNumOfMatchedStations; // >
        int cNormChiSquare;        // <
        float cd0_PV;              // <
        float cdz_PV;              // <
        
        bool evaluate(VarSet& vars);
        TString string();
};

class SynchEventSelectionCuts : public Cut
{
    public:
        SynchEventSelectionCuts();
        SynchEventSelectionCuts(float cDimuMassMin, float cDimuMassMax, float cTrigMuPtMin, float cTrigMuEtaMax, 
                                float cPVzMax, int cNDFpv, int cNPV);

        float cDimuMassMin;        // >
        float cDimuMassMax;        // <
        float cTrigMuPtMin;        // >
        float cTrigMuEtaMax;       // <
        float cPVzMax;             // < 
        int cNDFpv;                // >
        int cNPV;                  // > 

        bool evaluate(VarSet& vars);
        TString string();
};

class SynchMuonSelectionCuts : public Cut
{
    public:
        SynchMuonSelectionCuts();
        SynchMuonSelectionCuts(float minPt, float maxEta, float maxRelIso);

        float cMinPt;               // >
        float cMaxEta;              // <
        float cMaxRelIso;           // <

        bool evaluate(VarSet& vars);
        TString string();
};

class Run1EventSelectionCuts : public Cut
{

};

class Run1MuonSelectionCuts : public Cut
{

};
#endif
