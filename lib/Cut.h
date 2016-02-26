//Cut.h

#ifndef ADD_CUT
#define ADD_CUT


#include "VarSet.h"
#include "TString.h"

// Define the Interface
class Cut
{
    public:

        virtual bool evaluate(VarSet& vars) = 0;  // see whether the set of variables passes the cut
        TString tstring;                          // represent the cut as a tstring, may be used in ttree->Draw("...", cut);
};

#endif
