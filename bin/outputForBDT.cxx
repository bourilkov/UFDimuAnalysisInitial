// output the events in a small window around the higgs mass = 125 GeV and their variables
// so that we can train a bdt to make our categories for us

// Missing HLT trigger info in CMSSW_8_0_X MC so we have to compare Data and MC in a different manner.
// We apply triggers to data but not to MC. Then scale MC for trigger efficiency.

#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "CutSet.h"
#include "Cut.h"
#include "SelectionCuts.h"
#include "CategorySelection.h"
#include "JetSelectionTools.h"
#include "MuonSelectionTools.h"
#include "ElectronSelectionTools.h"
#include "TauSelectionTools.h"

#include "EventTools.h"
#include "PUTools.h"
#include "SignificanceMetrics.cxx"

#include "TLorentzVector.h"

#include <sstream>
#include <map>
#include <vector>
#include <utility>

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

UInt_t getNumCPUs()
{
  SysInfo_t s;
  gSystem->GetSysInfo(&s);
  UInt_t ncpu  = s.fCpus;
  return ncpu;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    // save the errors for the histogram correctly so they depend upon 
    // the number used to fill originally rather than the scaling
    TH1::SetDefaultSumw2();

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
    }   

    // Not sure that we need a map if we have a vector
    // Should use this as the main database and choose from it to make the vector
    std::map<std::string, Sample*> samples;

    // Second container so that we can have a copy sorted by cross section.
    std::vector<Sample*> samplevec;

    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    float reductionFactor = 1;
    float luminosity = 36814;      // pb-1
    GetSamples(samples, "UF");

    ///////////////////////////////////////////////////////////////////
    // PREPROCESSING---------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // Loop through all of the samples to do some pre-processing
    std::cout << std::endl;
    std::cout << "======== Preprocess the samples... " << std::endl;
    std::cout << std::endl;

    //makePUHistos(samples);
    
    for(auto &i : samples)
    {
        if(i.second->sampleType == "data") continue;

        // Output some info about the current file
        std::cout << "  /// Preprocessing " << i.second->name << std::endl;
        std::cout << std::endl;
        std::cout << "    sample name:       " << i.second->name << std::endl;
        std::cout << "    sample file:       " << i.second->filename << std::endl;
        std::cout << "    pileup file:       " << i.second->pileupfile << std::endl;
        std::cout << "    nOriginal:         " << i.second->nOriginal << std::endl;
        std::cout << "    N:                 " << i.second->N << std::endl;
        std::cout << "    nOriginalWeighted: " << i.second->nOriginalWeighted << std::endl;
        std::cout << std::endl;

        i.second->setBranchAddresses(2);
        samplevec.push_back(i.second);
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 

    std::cout << "@@@ nCPUs Available: " << getNumCPUs() << std::endl;
    std::cout << "@@@ nCPUs used     : " << nthreads << std::endl;
    std::cout << "@@@ nSamples used  : " << samplevec.size() << std::endl;

    auto outputSampleInfo = [luminosity, reductionFactor](Sample* s)
    {
      // Output some info about the current file
      std::cout << Form("  /// Processing %s \n", s->name.Data());

      std::map<TString, double> vars;
      vars["is_signal"] = -999;
      vars["weight"] = -999;
      vars["dimu_pt"] = -999;
      vars["mu0_pt"] = -999;
      vars["mu1_pt"] = -999;
      vars["mu0_eta"] = -999;
      vars["mu1_eta"] = -999;
      vars["jet0_pt"] = -999;
      vars["jet1_pt"] = -999;
      vars["jet0_eta"] = -999;
      vars["jet1_eta"] = -999;
      vars["N_valid_jets"] = -999;
      vars["m_jj"] = -999;
      vars["dEta_jj"] = -999;
      vars["dEta_jj_mumu"] = -999;
      vars["N_valid_extra_muons"] = -999;
      vars["extra_muon0_pt"] = -999;
      vars["extra_muon1_pt"] = -999;
      vars["extra_muon0_eta"] = -999;
      vars["extra_muon1_eta"] = -999;
      vars["N_valid_electrons"] = -999;
      vars["electron0_pt"] = -999;
      vars["electron1_pt"] = -999;
      vars["electron0_eta"] = -999;
      vars["electron1_eta"] = -999;
      vars["N_valid_extra_leptons"] = -999;
      vars["N_valid_bjets"] = -999;
      vars["bjet0_pt"] = -999;
      vars["bjet1_pt"] = -999;
      vars["bjet0_eta"] = -999;
      vars["bjet1_eta"] = -999;
      vars["dEta_bb"] = -999;
      vars["m_bb"] = -999;
      vars["mT_b_MET"] = -999;
      vars["MET"] = -999;
      vars["zep"] = -999;
      vars["phi_star"] = -999;
      vars["dPhi"] = -999;

      // !!!! output first line of csv to file
      std::ofstream file("csv/"+s->name+"_categorization_trianing.csv", std::ofstream::out);
      file << EventTools::outputMapKeysCSV(vars).Data() << std::endl;

      // Objects to help with the cuts and selections
      JetCollectionCleaner      jetCollectionCleaner;
      MuonCollectionCleaner     muonCollectionCleaner;
      EleCollectionCleaner      eleCollectionCleaner;

      Run2MuonSelectionCuts     run2MuonSelection;
      Run2EventSelectionCuts80X run2EventSelectionMC;

      for(unsigned int i=0; i<s->N/reductionFactor; i++)
      {
        // only load essential information for the first set of cuts 
        s->branches.recoDimuCands->GetEntry(i);
        s->branches.recoMuons->GetEntry(i);

        // loop and find a good dimuon candidate
        if(s->vars.recoDimuCands->size() < 1) continue;
        bool found_good_dimuon = false;

        // find the first good dimuon candidate and fill info
        for(auto& dimu: (*s->vars.recoDimuCands))
        {
          // the dimuon candidate and the muons that make up the pair
          s->vars.dimuCand = &dimu;
          MuonInfo& mu1 = s->vars.recoMuons->at(s->vars.dimuCand->iMu1);
          MuonInfo& mu2 = s->vars.recoMuons->at(s->vars.dimuCand->iMu2);

          // only want to train on events that are in the higgs mass window
          if(!(dimu.mass_PF > 120 && dimu.mass_PF < 130))
          {
              continue;
          }
          // usual cuts
          if(!mu1.isTightID || !mu2.isTightID)
          { 
              continue; 
          }
          if(!run2EventSelectionMC.evaluate(s->vars))
          { 
              continue; 
          }
          if(!run2MuonSelection.evaluate(s->vars)) 
          {
              continue; 
          }

          // dimuon event passes selections, set flag to true so that we only fill info for
          // the first good dimu candidate
          found_good_dimuon = true;

          // Load the rest of the information needed for run2 categories
          s->branches.jets->GetEntry(i);
          s->branches.mht->GetEntry(i);
          s->branches.nVertices->GetEntry(i);
          s->branches.recoElectrons->GetEntry(i);

          s->branches.gen_wgt->GetEntry(i);
          s->branches.nPU->GetEntry(i);
          s->branches.pu_wgt->GetEntry(i);
          s->branches.eff_wgt->GetEntry(i);

          // clear vectors for the valid collections
          s->vars.validMuons.clear();
          s->vars.validExtraMuons.clear();
          s->vars.validElectrons.clear();
          s->vars.validJets.clear();
          s->vars.validBJets.clear();

          // load valid collections from s->vars raw collections
          jetCollectionCleaner.getValidJets(s->vars, s->vars.validJets, s->vars.validBJets);
          muonCollectionCleaner.getValidMuons(s->vars, s->vars.validMuons, s->vars.validExtraMuons);
          eleCollectionCleaner.getValidElectrons(s->vars, s->vars.validElectrons);

          // Clean jets and electrons from muons, then clean remaining jets from remaining electrons
          CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validMuons, 0.3);
          CollectionCleaner::cleanByDR(s->vars.validElectrons, s->vars.validMuons, 0.3);
          CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validElectrons, 0.3);

          if(s->vars.validExtraMuons.size() + s->vars.validElectrons.size() > 0) continue;

          vars["dimu_pt"] = dimu.pt_PF;
          vars["mu0_pt"] = mu1.pt_PF;
          vars["mu1_pt"] = mu2.pt_PF;
          vars["mu0_eta"] = mu1.eta;
          vars["mu1_eta"] = mu2.eta;

          vars["N_valid_jets"] = s->vars.validJets.size();
          if(s->vars.validJets.size() >= 1) vars["jet0_pt"] = s->vars.validJets[0].Pt();
          if(s->vars.validJets.size() >= 2) vars["jet1_pt"] = s->vars.validJets[1].Pt();

          if(s->vars.validJets.size() >= 1) vars["jet0_eta"] = s->vars.validJets[0].Eta();
          if(s->vars.validJets.size() >= 2) vars["jet1_eta"] = s->vars.validJets[1].Eta();
          if(s->vars.validJets.size() >= 2)
          {
              TLorentzVector dijet = s->vars.validJets[0] + s->vars.validJets[1];
              float dEta = s->vars.validJets[0].Eta() - s->vars.validJets[1].Eta();
              vars["m_jj"] = dijet.M();
              vars["dEta_jj"] = TMath::Abs(dEta);
              float dEtajjmumu = dijet.Eta() - s->vars.dimuCand.recoCandEtaPF;
              vars["dEta_jj_mumu"] = TMath::Abs(dEtajjmumu);
          }

          vars["N_valid_extra_muons"] = s->vars.validExtraMuons.size();
          if(s->vars.validExtraMuons.size() >= 1) vars["extra_muon0_pt"] = s->vars.validExtraMuons[0].Pt();
          if(s->vars.validExtraMuons.size() >= 2) vars["extra_muon1_pt"] = s->vars.validExtraMuons[1].Pt();
          if(s->vars.validExtraMuons.size() >= 1) vars["extra_muon0_eta"] = s->vars.validExtraMuons[0].Eta();
          if(s->vars.validExtraMuons.size() >= 2) vars["extra_muon1_eta"] = s->vars.validExtraMuons[1].Eta();
          
          vars["N_valid_electrons"] = s->vars.validExtraMuons.size();
          if(s->vars.validElectrons.size() >= 1) vars["extra_electron0_pt"] = s->vars.validElectrons[0].Pt();
          if(s->vars.validElectrons.size() >= 2) vars["extra_electron1_pt"] = s->vars.validElectrons[1].Pt();
          if(s->vars.validElectrons.size() >= 1) vars["extra_electron0_eta"] = s->vars.validElectrons[0].Eta();
          if(s->vars.validElectrons.size() >= 2) vars["extra_electron1_eta"] = s->vars.validElectrons[1].Eta();

          vars["N_valid_extra_leptons"] = s->vars.validElectrons.size() + s->vars.validExtraMuons.size();

          vars["N_valid_bjets"] = s->vars.validBJets.size();
          if(s->vars.validBJets.size() >= 1) vars["bjet0_pt"] = s->vars.validBJets[0].Pt();
          if(s->vars.validBJets.size() >= 2) vars["bjet1_pt"] = s->vars.validBJets[1].Pt();

          if(s->vars.validBJets.size() >= 1) vars["bjet0_eta"] = s->vars.validBJets[0].Eta();
          if(s->vars.validBJets.size() >= 2) vars["bjet1_eta"] = s->vars.validBJets[1].Eta();
          if(s->vars.validBJets.size() >= 2)
          {
              TLorentzVector dibjet = s->vars.validBJets[0] + s->vars.validBJets[1];
              float dEta = s->vars.validBJets[0].Eta() - s->vars.validBJets[1].Eta();
              vars["m_bb"] = dibjet.M();
              vars["dEta_bb"] = TMath::Abs(dEta);
          }

          //if(s->vars.validBJets.size() > 0)
          //{
          //    TLorentzVector met(s->vars.met.px, s->vars.met.py, 0, s->vars.met.sumEt);
          //    TLorentzVector bjet = s->vars.validBJets[0];
          //    TLorentzVector bjet_t(bjet.Px(), bjet.Py(), 0, bjet.Et());
          //    TLorentzVector bmet_t = met + bjet_t;
          //    
          //    vars["mT_b_MET"] = bmet_t.M();
          //}

          vars["MET"] = s->vars.mht.pt;

          if(false)
            EventTools::outputEvent(s->vars, categorySelection);

          // need to put in is_signal and weight and the extra fancy variables
          vars["is_signal"] = s->sampleType.Contains("signal")?1:0;
          vars["weight"] = s->getScaleFactor(luminosity)*triggerSF*s->getWeight();

          // !!!! output event to file
          file << EventTools::outputMapValuesCSV(vars).Data() << std::endl;

          if(found_good_dimuon) break; // only fill one dimuon, break from dimu cand loop

        } // end dimucand loop
      } // end event loop

      file.close();
      std::cout << Form("  /// Done processing %s \n", s->name.Data());
      return 0;
    }; // end sample lambda function
}
