//CutSet.h

#ifndef ADD_CUTSET
#define ADD_CUTSET

#include "Cut.h"

// Define the Interface
class CutSet
{
    public:

        TString tstring;                         // represent the set of cuts as a tstring
        std::vector<Cut> cuts;                   // the collection of cuts

        bool evaluate(VarSet& vars)              // see whether the event passes the set of cuts
        {
            for(unsigned int i=0; i<cuts.size(); i++)
            {
                // if the event fails any of the cuts then return false
                if(!cuts[i].evaluate(vars)) return false;
            } 
            // if the event passes all of the cuts return true
            return true;
        };
};

#endif
