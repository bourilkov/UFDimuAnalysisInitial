//JetSelectionTools.h

#ifndef ADD_JETSELECTIONTOOLS
#define ADD_JETSELECTIONTOOLS

#include "VarSet.h"

class JetSelectionTools
{
    public:
        JetSelectionTools();
        JetSelectionTools(float cJetSelectionPtMin, float cJetSelectionEtaMax);

        float cJetSelectionPtMin; 
        float cJetSelectionEtaMax;       
 
        int getNValidJets(_PFJetInfo& jets);
};

#endif
