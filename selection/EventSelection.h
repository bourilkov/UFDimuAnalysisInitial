//EventSelection.h

#ifndef ADD_EVENTSELECTIONCUTS
#define ADD_EVENTSELECTIONCUTS

#include "Cut.h"

class Run1EventSelectionCuts80X : public Cut
{
// CMSSW 8_0_X is missing HLT Info so we can't use HLT in event selection
// We have to apply other cuts and scale for trigger efficiency

    public:
        Run1EventSelectionCuts80X();
        Run1EventSelectionCuts80X(bool isData);
        Run1EventSelectionCuts80X(float trigMuPtMin, float dimuMassMin);
        Run1EventSelectionCuts80X(bool isData, float trigMuPtMin, float dimuMassMin);

        bool isData = 0;
        float cTrigMuPtMin;        // >
        float cDimuMassMin;        // >

        void makeCutSet();
        bool evaluate(VarSet& vars);
        TString string();
};

class Run1EventSelectionCuts : public Cut
{
    public:
        Run1EventSelectionCuts();
        Run1EventSelectionCuts(float trigMuPtMin, float dimuMassMin);

        float cTrigMuPtMin;        // >
        float cDimuMassMin;        // >

        void makeCutSet();
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

        void makeCutSet();
        bool evaluate(VarSet& vars);
        bool passesVertexSelection(_VertexInfo& vertices);
        TString string();
};

class Run1EventSelectionSigCuts : public Cut
{
    public:
        Run1EventSelectionSigCuts();
        Run1EventSelectionSigCuts(float trigMuPtMin, float dimuMassMin, float dimuMassMax);

        float cTrigMuPtMin;        // >
        float cDimuMassMin;        // >
        float cDimuMassMax;        // <

        void makeCutSet();
        bool evaluate(VarSet& vars);
        TString string();
};

class FEWZCompareCuts : public Cut
{
    public:
        FEWZCompareCuts();
        FEWZCompareCuts(bool useReco);
        FEWZCompareCuts(bool useReco, float leadPtMin, float subleadPtMin, float maxEta, float dimuMassMin, float dimuMassMax, float maxRelIso);

        bool useReco;

        float cLeadPtMin;           
        float cSubleadPtMin;        
        float cMaxEta;              
        float cDimuMassMin;         
        float cDimuMassMax;         
        float cMaxRelIso;         

        void makeCutSet();
        bool evaluate(VarSet& vars);
        bool evaluate(_MuonInfo& recoMu, int m);
        TString string();
};
#endif
