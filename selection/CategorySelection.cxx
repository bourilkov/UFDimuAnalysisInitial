///////////////////////////////////////////////////////////////////////////
//                             CategorySelection.cxx                     //
//=======================================================================//
//                                                                       //
//        Define the selection critera for the different categories.     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "CategorySelection.h"

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________CategorySelection______________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

CategorySelection::CategorySelection()
{
// Standard values for the tight muon id cuts

    // Preselection
    cLeadPtMin = 40;
    cSubleadPtMin = 30;
    cMETMax = 40;
    isPreselected = false;

    // VBF Tight
    cDijetMassMinVBFT = 650;
    cDijetDeltaEtaMin = 3.5;
    isVBFTight = false;

    // VBF Loose
    isVBFLoose = false;

    // GGF Tight
    cDijetMassMinGGFT = 250;
    cDimuPtMinGGFT = 50;
    isGGFTight = false;

    // 01Tight
    cDimuPtMin01T = 10;
    isTight01 = false;

    // 01Loose
    isLoose01 = false;
}
