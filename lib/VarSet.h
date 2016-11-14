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
        _DimuCandInfo dimuCand;
        _ElectronInfo recoElectrons;
        _TauInfo recoTaus;
        _MetInfo met;
        _PFJetInfo jets;

        std::vector<TLorentzVector> validMuonsDecoy;
        std::vector<TLorentzVector> validMuons;
        std::vector<TLorentzVector> validExtraMuons;
        std::vector<TLorentzVector> validElectrons;
        std::vector<TLorentzVector> validTaus;
        std::vector<TLorentzVector> validJets;
        std::vector<TLorentzVector> validBJets;

        // gen info
        _GenJetInfo  genJets;
        std::vector<TLorentzVector> validGenJets;
        _genPartInfo genZpreFSR, genZpostFSR, genGpreFSR, genGpostFSR;
        _TrackInfo genM1GpreFSR, genM1GpostFSR, genM2GpreFSR, genM2GpostFSR; // muons from virtual photon
        _TrackInfo genM1ZpreFSR, genM1ZpostFSR, genM2ZpreFSR, genM2ZpostFSR; // muons from Z

        int genWeight;
        int nPU;
 
};

#endif
