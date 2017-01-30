//EventSelection.h

#ifndef ADD_EVENTSELECTIONCUTS
#define ADD_EVENTSELECTIONCUTS

#include "Cut.h"

class Run2EventSelectionCuts80X : public Cut
{
// if data apply hlt trigger matching
// if MC do not apply hlt trigger matching
    public:
        Run2EventSelectionCuts80X();
        Run2EventSelectionCuts80X(bool isData);
        Run2EventSelectionCuts80X(float trigMuPtMin, float dimuMassMin);
        Run2EventSelectionCuts80X(bool isData, float trigMuPtMin, float dimuMassMin);

        bool isData = 0;
        float cTrigMuPtMin;        // >
        float cDimuMassMin;        // >

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
