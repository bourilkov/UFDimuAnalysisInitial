///////////////////////////////////////////////////////////////////////////
//                         EleCollectionCleaner.h                        //
//=======================================================================//
//                                                                       //
//        Select valid electrons.                                        //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ADD_ELECTRONSELECTIONTOOLS
#define ADD_ELECTRONSELECTIONTOOLS

#include "VarSet.h"
#include <vector>
#include "CollectionCleaner.hxx"

class EleCollectionCleaner : public CollectionCleaner
{
    public:
        EleCollectionCleaner();
        EleCollectionCleaner(float cElectronSelectionPtMin, float cElectronSelectionEtaMax, float cElectronSelectionIsoMax, int cElectronSelectionID);

        float cElectronSelectionPtMin; 
        float cElectronSelectionEtaMax;       
        float cElectronSelectionIsoMax;
        int   cElectronSelectionID;

        void getValidElectrons(VarSet& vars, std::vector<TLorentzVector>& evec);
};

#endif
