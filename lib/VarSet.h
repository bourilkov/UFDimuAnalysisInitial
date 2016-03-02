// VarSet.h
// Load the variables from the ttree into the structs in this class
// The structs are defined in DataFormats.h

#ifndef ADD_VARSET
#define ADD_VARSET

#include "DataFormats.h"
#include "TLorentzVector.h"

class VarSet
{
    public:
        VarSet(){};
        ~VarSet(){};

        _EventInfo eventInfo;
        _VertexInfo vertices;
        _MuonInfo reco1;
        _MuonInfo reco2;
        _MetInfo met;
        _PFJetInfo jets;
        std::vector<TLorentzVector> validJets;
 
        float recoCandMass, recoCandMassPF;
        float recoCandPt, recoCandPtPF;
        float rho;
        int genWeight;
        int nPU;
};

#endif
