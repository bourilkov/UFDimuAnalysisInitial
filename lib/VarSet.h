// VarSet.h
// Load the variables from the ttree into the structs in this class
// The structs are defined in DataFormats.h

#ifndef ADD_VARSET
#define ADD_VARSET

#include "DataFormats.h"
#include "EventInfo.h"
#include "VertexInfo.h"
#include "MuonInfo.h"
#include "PairInfo.h"
#include "EleInfo.h"
#include "MetInfo.h"
#include "JetInfo.h"
#include "TLorentzVector.h"

class VarSet
{
    public:
        VarSet(){};
        ~VarSet(){};

        // reco info
        EventInfo eventInfo;
        PairInfo dimuCand;
        MetInfo met;

        std::vector<VertexInfo>* vertices = 0;
        std::vector<MuonInfo>* recoMuons = 0;
        std::vector<EleInfo>* recoElectrons = 0;
        std::vector<JetInfo>* jets = 0;

        std::vector<TLorentzVector> validMuonsDecoy;
        std::vector<TLorentzVector> validMuons;
        std::vector<TLorentzVector> validExtraMuons;
        std::vector<TLorentzVector> validElectrons;
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
