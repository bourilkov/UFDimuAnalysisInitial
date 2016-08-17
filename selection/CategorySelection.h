//CategorySelection.h

#ifndef ADD_CATEGORYSELECTION
#define ADD_CATEGORYSELECTION

#include "VarSet.h"
#include "JetSelectionTools.h"


class Categorizer
{
    public:
        // Determine which category the event belongs to
        virtual void evaluate(VarSet& vars) = 0;
        // reset the boolean values for the categories
        virtual void reset() = 0;
};

class CategorySelection : public Categorizer
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

class CategorySelectionFEWZ : public Categorizer
{
    public:
        CategorySelectionFEWZ(); 
        CategorySelectionFEWZ(float massSplit, float etaCentralSplit, float jetPtMin, float jetEtaMax); 

        // Selections
        float cMassSplit;
        float cEtaCentralSplit;
        float cJetPtMin;
        float cJetEtaMax;

        // Initial Tests
        bool isWide;
        bool isNarrow;
        bool isCentralCentral;
        bool isCentralNotCentral;
        bool isOneJetInclusive;

        // Final Categories
        bool isCentralCentralWide;
        bool isCentralCentralNarrow;

        bool isCentralNotCentralWide;
        bool isCentralNotCentralNarrow;

        bool isOneJetInclusiveWide;
        bool isOneJetInclusiveNarrow;

        // Determine which category the event belongs to
        void evaluate(VarSet& vars);
        void reset();
};
#endif
