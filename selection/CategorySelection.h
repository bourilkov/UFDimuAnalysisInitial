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
        CategorySelection(float cLeadPtMin, float cSubleadPtMin, float cMETMax, float cDijetMassMinVBFT, float cDijetDeltaEtaMin, float cDijetMassMinGGFT,
                          float cDimuPtMinGGFT, float cDimuPtMin01T); 

        // Preselection
        float cLeadPtMin;
        float cSubleadPtMin;
        float cMETMax;

        // Jet selection
        float cJetSelectionPtMin;  
        float cJetSelectionEtaMax;      

        bool isPreselected;

        // VBF Tight
        float cDijetMassMinVBFT;
        float cDijetDeltaEtaMin;
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
        void evaluate(VarSet& vars);
        
        // boolean tests for the different categories
        bool preselection(VarSet& vars);
        bool VBFTight(VarSet& vars);
        bool GGFTight(VarSet& vars);
        bool VBFLoose(VarSet& vars);
        bool Tight01(VarSet& vars);
        bool Loose01(VarSet& vars);
};
#endif
