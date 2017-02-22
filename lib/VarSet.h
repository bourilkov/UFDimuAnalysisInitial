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
#include <unordered_map>
#include <map>

class VarSet
{
    public:
        VarSet();
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

        // map variable string to variable value via function pointers
        // uses functions below
        std::unordered_map<std::string, double(VarSet::*)()> varMap;

        // get the the value for some variable in one of the structs above
        // by name (string). Use the varMap to get the appropriate function
        // of those below.
        double getValue(const std::string& name) 
        {
          // The specific object doesn't have its own function until
          // some method in the object is called.
          if(varMap[name])
            return (this->*varMap[name])();
          else return -999;
        }

        // muon variables
        double dimu_pt(){  return dimuCand->pt_PF;                     };
        double mu1_pt() {  return recoMuons->at(dimuCand->iMu1).pt_PF; };
        double mu2_pt() {  return recoMuons->at(dimuCand->iMu2).pt_PF; };
        double mu1_eta(){  return recoMuons->at(dimuCand->iMu1).eta;   };
        double mu2_eta(){  return recoMuons->at(dimuCand->iMu2).eta;   };

        // jet variables
        double jet0_pt() { return (validJets.size()>=1)?validJets[0].Pt():-999;  };
        double jet1_pt() { return (validJets.size()>=2)?validJets[1].Pt():-999;  };
        double jet0_eta(){ return (validJets.size()>=1)?validJets[0].Eta():-999; };
        double jet1_eta(){ return (validJets.size()>=2)?validJets[1].Eta():-999; };

        double m_jj()        { return (validJets.size()>=2)?(validJets[0]+validJets[1]).M():-999; };
        double dEta_jj()     { return (validJets.size()>=2)?TMath::Abs(validJets[0].Eta()-validJets[1].Eta()):-999; };
        double dEta_jj_mumu(){ return (validJets.size()>=2)?TMath::Abs((validJets[0]+validJets[1]).Eta()-dimuCand->eta):-999; };

        // bjet variables
        double bjet0_pt() { return (validBJets.size()>=1)?validBJets[0].Pt():-999;  };
        double bjet1_pt() { return (validBJets.size()>=2)?validBJets[1].Pt():-999;  };
        double bjet0_eta(){ return (validBJets.size()>=1)?validBJets[0].Eta():-999; };
        double bjet1_eta(){ return (validBJets.size()>=2)?validBJets[1].Eta():-999; };

        double m_bb()   { return (validBJets.size()>=2)?(validBJets[0]+validBJets[1]).M():-999; };
        double dEta_bb(){ return (validBJets.size()>=2)?TMath::Abs(validBJets[0].Eta()-validBJets[1].Eta()):-999; };

        // # variables
        double N_valid_jets()         { return validJets.size();          };
        double N_valid_bjets()        { return validBJets.size();         };
        double N_valid_extra_muons()  { return validExtraMuons.size();    };
        double N_valid_electrons()    { return validElectrons.size();     };
        double N_valid_extra_leptons(){ return validExtraMuons.size() + validElectrons.size(); };

        // MET
        double MET(){ return mht->pt; };

        // special variables, to be implemented
        double zep()     { return -999; };
        double phi_star(){ return -999; };
        double dPhi()    { return -999; };
};

#endif
