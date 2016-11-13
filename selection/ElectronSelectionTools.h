//ElectronSelectionTools.h

#ifndef ADD_ELECTRONSELECTIONTOOLS
#define ADD_ELECTRONSELECTIONTOOLS

#include "VarSet.h"
#include <vector>

class ElectronSelectionTools
{
    public:
        ElectronSelectionTools();
        ElectronSelectionTools(float cElectronSelectionPtMin, float cElectronSelectionEtaMax, float cElectronSelectionIsoMax, int cElectronSelectionID);

        float cElectronSelectionPtMin; 
        float cElectronSelectionEtaMax;       
        float cElectronSelectionIsoMax;
        int   cElectronSelectionID;

        void getValidElectrons(VarSet& vars, std::vector<TLorentzVector>& evec);
};

#endif
