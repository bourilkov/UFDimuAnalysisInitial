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
#include "GenParentInfo.h"
#include "GenMuonInfo.h"
#include "GenPairInfo.h"
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
        EventInfo* eventInfo = 0;
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
        std::vector<GenParentInfo>* genParents = 0;
        std::vector<GenMuonInfo>* genMuons = 0;
        std::vector<GenPairInfo>* genDimuons = 0;

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
        double dimu_pt()   {  return dimuCand->pt_PF;                     };
        double mu1_pt()    {  return recoMuons->at(dimuCand->iMu1).pt_PF; };
        double mu2_pt()    {  return recoMuons->at(dimuCand->iMu2).pt_PF; };
        double mu1_eta()   {  return recoMuons->at(dimuCand->iMu1).eta;   };
        double mu2_eta()   {  return recoMuons->at(dimuCand->iMu2).eta;   };
        double mu_res_eta(){  return (TMath::Abs(mu1_eta()) + TMath::Abs(mu2_eta()))/2;   };

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

        double extra_muon0_pt() { return (validExtraMuons.size()>=1)?validExtraMuons[0].Pt():-999; };
        double extra_muon1_pt() { return (validExtraMuons.size()>=2)?validExtraMuons[1].Pt():-999; };
        double extra_muon0_eta(){ return (validExtraMuons.size()>=1)?validExtraMuons[0].Eta():-999; };
        double extra_muon1_eta(){ return (validExtraMuons.size()>=2)?validExtraMuons[1].Eta():-999; };

        double electron0_pt() { return (validElectrons.size()>=1)?validElectrons[0].Pt():-999; };
        double electron1_pt() { return (validElectrons.size()>=2)?validElectrons[1].Pt():-999; };
        double electron0_eta(){ return (validElectrons.size()>=1)?validElectrons[0].Eta():-999; };
        double electron1_eta(){ return (validElectrons.size()>=2)?validElectrons[1].Eta():-999; };

        double mT_b_MET()
        { 
            if(validBJets.size() < 1) return -999;
            TLorentzVector met(mht->pt*TMath::Cos(mht->phi), mht->pt*TMath::Cos(mht->phi), 0, mht->pt);                  
            TLorentzVector bjet = validBJets[0];
            TLorentzVector bjet_t(bjet.Px(), bjet.Py(), 0, bjet.Et());
            TLorentzVector bmet_t = met + bjet_t;
            return bmet_t.M();
        };

        // special variables, to be implemented
        double zep0()
        { 
            if(validJets.size() < 2) return -999; 
            double meanEta = (validJets[0].Eta() +validJets[1].Eta())/2;
            return validJets[0].Eta() - meanEta;
        };
        double zep1()
        { 
            if(validJets.size() < 2) return -999; 
            double meanEta = (validJets[0].Eta() +validJets[1].Eta())/2;
            return validJets[1].Eta() - meanEta;
        };
        double phi_star()
        {
            double phi_star = 0;
            double mu_dPhi = TMath::Abs(recoMuons->at(dimuCand->iMu1).phi - recoMuons->at(dimuCand->iMu2).phi);
            if(mu_dPhi > TMath::Pi()) mu_dPhi = 2*TMath::Pi() - mu_dPhi;
            double phiACOP = TMath::Pi() - mu_dPhi;
            double thetaStarEta = TMath::ACos(TMath::TanH((recoMuons->at(dimuCand->iMu1).eta - recoMuons->at(dimuCand->iMu2).eta)/2));
            phi_star = TMath::Tan(phiACOP/2)*TMath::Sin(thetaStarEta);
            return phi_star;
        };
        double dPhi_jj_mumu()    
        { 
            if(validJets.size() < 2) return -999;
            double dphi = TMath::Abs((validJets[0] + validJets[1]).Phi() - dimuCand->phi);
            if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
        };
};

#endif
