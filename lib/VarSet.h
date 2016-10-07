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

        // reco info
        _EventInfo eventInfo;
        _VertexInfo vertices;
        _MuonInfo recoMuons;
        _MetInfo met;
        _PFJetInfo jets;
        std::vector<TLorentzVector> validJets;

        float recoCandMass, recoCandMassPF;
        float recoCandPt, recoCandPtPF;
        float rho;
        int nPU;

        // gen info
        _GenJetInfo  genJets;
        std::vector<TLorentzVector> validGenJets;
        _genPartInfo genZpreFSR, genZpostFSR, genGpreFSR, genGpostFSR;
        _TrackInfo genM1GpreFSR, genM1GpostFSR, genM2GpreFSR, genM2GpostFSR; // muons from virtual photon
        _TrackInfo genM1ZpreFSR, genM1ZpostFSR, genM2ZpreFSR, genM2ZpostFSR; // muons from Z

        int genWeight;
 
};

#endif
