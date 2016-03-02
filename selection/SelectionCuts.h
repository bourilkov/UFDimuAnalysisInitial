//SelectionCuts.h

#ifndef ADD_SELECTIONCUTS
#define ADD_SELECTIONCUTS

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
                                float cPVzMax, int cNDFpv, int cNPV, int nJets);

        float cDimuMassMin;        // >
        float cDimuMassMax;        // <
        float cTrigMuPtMin;        // >
        float cTrigMuEtaMax;       // <
        float cPVzMax;             // < 
        int cNDFpv;                // >
        int cNPV;                  // > 
        int cNJets;                // <=

        bool evaluate(VarSet& vars);
        bool passesVertexSelection(_VertexInfo& vertices);
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
        bool evaluate(_MuonInfo& recoMu, float rho);
        TString string();
};

class Run1EventSelectionCuts : public Cut
{
    public:
        Run1EventSelectionCuts();
        Run1EventSelectionCuts(float trigMuPtMin);

        float cTrigMuPtMin;        // >

        bool evaluate(VarSet& vars);
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

        bool evaluate(VarSet& vars);
        bool evaluate(_MuonInfo& recoMu);
        TString string();
};
#endif
