///////////////////////////////////////////////////////////////////////////
//                             SelectionCuts.cxx                         //
//=======================================================================//
//                                                                       //
//        Define the Event and Muon Selection cuts for the analysis.     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "SelectionCuts.h"

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________TightMuonID ___________________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

TightMuonIdCuts::TightMuonIdCuts()
{
// Standard values for the tight muon id cuts

    cNumTrackerLayers = 5;     // >
    cNumValidMuonHits = 0;     // >
    cNumValidPixelHits = 0;    // >
    cNumOfMatchedStations = 1; // >
    cNormChiSquare = 10;       // <
    cd0_PV = 0.2;              // <
    cdz_PV = 0.5;              // <
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TightMuonIdCuts::TightMuonIdCuts(int numTrackerLayers, int numValidMuonHits, int numValidPixelHits, int numOfMatchedStations, 
int normChiSquare, float d0_pV, float dz_PV)
{
// Custom values for the tight muon id cuts

    cNumTrackerLayers = numTrackerLayers;
    cNumValidMuonHits = numValidMuonHits;
    cNuMValidPixelHits = numValidPixelHits;
    cNumOfMatchedStations = numOfMatchedStations;
    cNormChiSquare = normChiSquare;
    cd0_PV = d0_PV;
    cdz_PV = dz_PV;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool TightMuonIdCuts::evaluate(VarSet& vars)
{
    // Cuts on reco1
    // if the event fails a single cut return false
    if(!vars.reco1.isGlobal) return false;
    if(!vars.reco1.isPFMuon) return false;
    if(!(vars.reco1.numTrackerLayers > cNumTrackerLayers)) return false;
    if(!(vars.reco1.numValidMuonHits > cNumValidMuonHits)) return false;
    if(!(vars.reco1.numValidPixelHits > cNumValidPixelHits)) return false;
    if(!(vars.reco1.numOfMatchedStations > cNumOfMatchedStations)) return false;
    if(!(vars.reco1.normChiSquare < cNormChiSquare)) return false;
    if(!(vars.reco1.d0_PV < cd0_PV)) return false;
    if(!(vars.reco1.dz_PV < cdz_PV)) return false;

    // Cuts on reco2
    if(!vars.reco2.isGlobal) return false;
    if(!vars.reco2.isPFMuon) return false;
    if(!(vars.reco2.numTrackerLayers > cNumTrackerLayers)) return false;
    if(!(vars.reco2.numValidMuonHits > cNumValidMuonHits)) return false;
    if(!(vars.reco2.numValidPixelHits > cNumValidPixelHits)) return false;
    if(!(vars.reco2.numOfMatchedStations > cNumOfMatchedStations)) return false;
    if(!(vars.reco2.normChiSquare < cNormChiSquare)) return false;
    if(!(vars.reco2.d0_PV < cd0_PV)) return false;
    if(!(vars.reco2.dz_PV < cdz_PV)) return false;

    // The event passed all the cuts return true
    return true;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TString TightMuonIdCuts::string
{
    return TString("Tight Muon ID");
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________SynchEventSelection____________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

SynchEventSelectionCuts::SynchEventSelectionCuts()
{

}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

#endif
