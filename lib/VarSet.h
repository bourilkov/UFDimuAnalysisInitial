// VarSet.h
// Load the variables from the ttree into the structs in this class
// The structs are defined in DataFormats.h

#ifndef ADD_VARSET
#define ADD_VARSET

#include "DataFormats.h"
#include "EventInfo.h"
#include "MuonInfo.h"
#include "PairInfo.h"
#include "EleInfo.h"
#include "MhtInfo.h"
#include "SlimJetInfo.h"
#include "TLorentzVector.h"

class VarSet
{
    public:
        VarSet(){};
        ~VarSet(){};

        // reco weights
        float eff_wgt;
        float pu_wgt;

        // reco info
        Int_t nVertices;
        EventInfo eventInfo;
        PairInfo* dimuCand = 0; // this is a pointer to one of the dimu candidates in the vector
                                // we don't want to copy the object for ~40 million events
                                // the other objects/primitives are loaded via TBranch->Get();
        MhtInfo* mht = 0;

        std::vector<PairInfo>* recoDimuCands = 0;
        std::vector<MuonInfo>* recoMuons = 0;
        std::vector<EleInfo>* recoElectrons = 0;
        std::vector<SlimJetInfo>* jets = 0;

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

        int nPU;

        // gen weights
        int gen_wgt;
 
};

#endif
