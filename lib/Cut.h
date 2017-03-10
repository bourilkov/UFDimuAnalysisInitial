///////////////////////////////////////////////////////////////////////////
// ======================================================================//
// Cut.h                                                                 //
// ======================================================================//
// Abstract class to be inherited by specific cuts.                      //
// Check out ../selection/ for implementations. We inherit this class    //
// to make EventSelection and MuonSelection objects.                     //
// The object takes in the dimuons, muons, electrons, etc via VarSet     //
// then evaluates this information to see if they pass the set of cuts   //
// for the Cut object.                                                   //
// The cutset object has book keeping information for the different cuts.//
// ======================================================================//
///////////////////////////////////////////////////////////////////////////

#ifndef ADD_CUT
#define ADD_CUT

#include "CutSet.hxx"
#include "VarSet.h"
#include "TString.h"

// Define the Interface
class Cut
{
    public:
        CutSet cutset;                              // bookkeeping for all of the cuts in the cut object

        virtual bool evaluate(VarSet& vars) = 0;    // see whether the set of variables passes the cut
        virtual TString string() = 0;               // represent the cut as a tstring, may be used in ttree->Draw("...", cut)
                                                    // if possible 
        virtual void makeCutSet() = 0;              // book keeping for all of the cuts, useful for the N-1 plots, significance optimization
                                                    // and debugging
};

#endif
