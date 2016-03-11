//CategorySelection.h

#ifndef ADD_CATEGORYSELECTION
#define ADD_CATEGORYSELECTION

#include "VarSet.h"
#include "JetSelectionTools.h"

// Define the different cuts
class CategorySelection
{
    public:
        CategorySelection(); 
        CategorySelection(float cLeadPtMin, float cSubleadPtMin, float cMETMax, float cDijetMassMinVBFT, float cDijetDeltaEtaMinVBFT, float cDijetMassMinGGFT,
                          float cDimuPtMinGGFT, float cDimuPtMin01T); 

        // Preselection
        float cLeadPtMin;
        float cSubleadPtMin;
        float cMETMax;
        bool isPreselected;

        // VBF Tight
        float cDijetMassMinVBFT;
        float cDijetDeltaEtaMinVBFT;
        bool isVBFTight;

        // VBF Loose
        bool isVBFLoose;
     
        // GGF Tight
        float cDijetMassMinGGFT;
        float cDimuPtMinGGFT;
        bool isGGFTight;

        // 01Tight
        float cDimuPtMin01T; 
        bool isTight01;

        // 01Loose
        bool isLoose01;

        // Determine which category the event belongs to
        // result stored in isVBFTight, isGGFTight, etc 
        void evaluate(VarSet& vars);
        void reset();
};
#endif
