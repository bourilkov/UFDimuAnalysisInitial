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

int main(int argc, char* argv[])
{
    // save the errors for the histogram correctly so they depend upon 
    // the number used to fill originally rather than the scaling
    TH1::SetDefaultSumw2();

    unsigned int nPartitions = 1; // break the samples up into nPartitions to run in parallel
    unsigned int partition = 0;   // which partition to make the plots for

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) ss >> nPartitions;
        if(i==2) ss >> partition;
    }   

    // Not sure that we need a map if we have a vector
    // Should use this as the main database and choose from it to make the vector
    std::map<std::string, Sample*> samples;

    // Second container so that we can have a copy sorted by cross section.
    std::vector<Sample*> samplevec;

    // Use this to plot some things if we wish
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

  
    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    float luminosity = 33598;      // pb-1
    float triggerSF = 0.913;       // no HLT trigger info available for the samples so we scale for the trigger efficiency instead
    float signalSF = 100;          // not using this at the moment, but scale the signal samples to see them better in the plots if you want

    // ================================================================
    // Data -----------------------------------------------------------
    // ================================================================
 

    TString datafilename = 
    TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data/25ns/golden/CMSSW_8_0_X/stage_1_singleMuon_Run2016BCDEFGH_ALL_etm.root");

    Sample* datasample = new Sample(datafilename, "Data", "data");
    datasample->lumi = luminosity;
    datasample->xsec = 9999;
    datasample->pileupfile = "pu_reweight_trees/8_0_X/PU_2016BCDEFGH_xsec69p2mb_CMSSW_8_0_X.root";
    //datasample->pileupfile = "pu_reweight_trees/8_0_X/PU_2016B_xsec69mb_CMSSW_8_0_X.root";
    //datasample->pileupfile = "pu_reweight_trees/8_0_X/PU_2016BCD_ICHEP_xsec69p2mb_CMSSW_8_0_X.root";
    //datasample->pileupfile = "pu_reweight_trees/8_0_X/PU_2016BCDE_xsec69p2mb_CMSSW_8_0_X.root";
    //samples["Data"] = datasample;

    // ================================================================
    // DYJetsToLL -----------------------------------------------------
    // ================================================================

    TString dyfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/bg/dy/CMSSW_8_0_X/stage_1_dy_jetsToLL_ALL_etm.root");
    samples["DYJetsToLL"] = new Sample(dyfilename, "DYJetsToLL", "background");
    samples["DYJetsToLL"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_DYJetsToLL.root"; //nPU
    samples["DYJetsToLL"]->xsec = 6025.2; // pb

    // ================================================================
    // TTJets ---------------------------------------------------------
    // ================================================================

    TString ttbarfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/bg/ttbar/CMSSW_8_0_X/stage_1_ttJets_ALL_etm.root");
    samples["TTJets"] = new Sample(ttbarfilename, "TTJets", "background");
    samples["TTJets"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_TTJets.root"; //nPU
    samples["TTJets"]->xsec = 831.76; // pb

    // ================================================================
    // VBF ---------------------------------------------------------
    // ================================================================

    TString vbffilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/signal/CMSSW_8_0_X/stage_1_vbf_HToMuMu_ALL_etm.root");
    samples["VBF"] = new Sample(vbffilename, "VBF", "signal");
    samples["VBF"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_VBF.root"; //nPU
    samples["VBF"]->xsec = 3.727*0.00022; // pb

    // ================================================================
    // GGF ---------------------------------------------------------
    // ================================================================

    TString ggfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/signal/CMSSW_8_0_X/stage_1_gg_HToMuMu_ALL_etm.root");
    samples["GGF"] = new Sample(ggfilename, "GGF", "signal");
    samples["GGF"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_GGF.root"; //nPU
    samples["GGF"]->xsec = 43.62*0.00022; // pb

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

        if(!i.second->sampleType.Contains("data"))
        {
            // Pileup reweighting
            std::cout << "    +++ PU Reweighting " << i.second->name << "..."  << std::endl;
            std::cout << std::endl;

            i.second->lumiWeights = new reweight::LumiReWeighting(i.second->pileupfile.Data(), datasample->pileupfile.Data(), "pileup", "pileup");
            std::cout << "        " << i.first << "->lumiWeights: " << i.second->lumiWeights << std::endl;
            std::cout << std::endl;
        }
        samplevec.push_back(i.second);
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 
    

    ///////////////////////////////////////////////////////////////////
    // Cut and Categorize ---------------------------------------------
    ///////////////////////////////////////////////////////////////////
    
    // Objects to help with the cuts and selections
    JetSelectionTools      jetSelectionTools;
    MuonSelectionTools     muonSelectionTools;
    ElectronSelectionTools electronSelectionTools;
    TauSelectionTools      tauSelectionTools;

    LotsOfCategoriesRun2      categorySelection;
    Run1MuonSelectionCuts     run1MuonSelection;
    Run1EventSelectionCuts80X run1EventSelectionData(true);
    Run1EventSelectionCuts80X run1EventSelectionMC;

    TString varname;
    int bins;
    float min;
    float max;

    std::map<TString, double> vars;
    vars["is_signal"] = -999;
    vars["weight"] = -999;
    vars["dimu_pt"] = -999;
    vars["mu0_pt"] = -999;
    vars["mu1_pt"] = -999;
    vars["mu0_eta"] = -999;
    vars["mu1_eta"] = -999;
    vars["jj_jet0_pt"] = -999;
    vars["jj_jet1_pt"] = -999;
    vars["jj_jet0_eta"] = -999;
    vars["jj_jet1_eta"] = -999;
    vars["nValJets"] = -999;
    vars["m_jj"] = -999;
    vars["dEta_jj"] = -999;
    vars["dEta_jj_mumu"] = -999;
    vars["nExtraMu"] = -999;
    vars["extra_muon0_pt"] = -999;
    vars["extra_muon1_pt"] = -999;
    vars["extra_muon0_eta"] = -999;
    vars["extra_muon1_eta"] = -999;
    vars["nEle"] = -999;
    vars["electron0_pt"] = -999;
    vars["electron1_pt"] = -999;
    vars["electron0_eta"] = -999;
    vars["electron1_eta"] = -999;
    vars["nExtraLep"] = -999;
    vars["nValBTags"] = -999;
    vars["bjet0_pt"] = -999;
    vars["bjet1_pt"] = -999;
    vars["bjet0_eta"] = -999;
    vars["bjet1_eta"] = -999;
    vars["dEta_bb"] = -999;
    vars["m_bb"] = -999;
    vars["mT_b_MET"] = -999;
    vars["MET"] = -999;
    vars["zep"] = -999;
    vars["dimu_dPhiStar"] = -999;
    vars["dPhi"] = -999;

    // !!!! output first line of csv to file
    std::ofstream file("csv/categorization_trianing.csv", std::ofstream::out);
    file << EventTools::outputMapKeysCSV(vars).Data() << std::endl;
 
    std::cout << std::endl;
    std::cout << "======== Plot Configs ========" << std::endl;
    std::cout << "nPartitions : " << nPartitions << std::endl;
    std::cout << "partition   : " << partition << std::endl;
    std::cout << std::endl;

    for(auto &s : samplevec)
    {
      // Output some info about the current file
      unsigned int start = ((float)partition/(float)nPartitions)*s->N;
      unsigned int end   = ((float)(partition+1))/((float)nPartitions)*s->N;
      std::cout << std::endl;
      std::cout << "  /// Looping over " << s->name << std::endl;
      std::cout << "  ///    start       : " << start << std::endl;
      std::cout << "  ///    end         : " << end << std::endl;

      for(unsigned int i=start; i<end; i++)
      {

        ///////////////////////////////////////////////////////////////////
        // GET INFORMATION ------------------------------------------------
        ///////////////////////////////////////////////////////////////////

        //std::cout << "s->getEntry..." << std::endl;
        s->getEntry(i); 

        // initialize vectors for the good jets, bjets, muons, electrons, and taus
        s->vars.validMuons.clear();
        s->vars.validMuonsDecoy.clear();
        s->vars.validExtraMuons.clear();
        s->vars.validElectrons.clear();
        s->vars.validTaus.clear();
        s->vars.validJets.clear();
        s->vars.validBJets.clear();

        // try to catch errors in initializing the valid object vectors
        if(s->vars.validMuons.size() != 0 || s->vars.validMuonsDecoy.size() !=0 || s->vars.validExtraMuons.size() !=0 
          || s->vars.validElectrons.size() !=0 ||  s->vars.validTaus.size() !=0 || s->vars.validBJets.size() !=0 
          || s->vars.validJets.size() !=0)
        {
            std::cout << "!!! ERROR: clear valid object vectors should be zero..." << std::endl;
            std::cout << "mu size      : " << s->vars.validMuons.size() << std::endl;
            std::cout << "mu decoy size: " << s->vars.validMuonsDecoy.size() << std::endl;
            std::cout << "xmu size     : " << s->vars.validExtraMuons.size() << std::endl;
            std::cout << "e size       : " << s->vars.validElectrons.size() << std::endl;
            std::cout << "t size       : " << s->vars.validTaus.size() << std::endl;
            std::cout << "bjet size    : " << s->vars.validBJets.size() << std::endl;
            std::cout << "jet size     : " << s->vars.validJets.size() << std::endl;
            return 0;
        }
        // filter the good objects from the larger set
        // Need to clean collections by dR later

        jetSelectionTools.getValidJets(s->vars, s->vars.validJets);
        jetSelectionTools.getValidBJets(s->vars, s->vars.validBJets);
        muonSelectionTools.getValidMuons(s->vars, s->vars.validMuons);
        electronSelectionTools.getValidElectrons(s->vars, s->vars.validElectrons);
        tauSelectionTools.getValidTaus(s->vars, s->vars.validTaus);

        // segfaults at the moment, fix this
        // clean collections by dR
        // clean 4vecs in v1 by dR based upon 4vecs in v2
        // EventTools::cleanByDR(v1, v2, dRmin);
        //EventTools::cleanByDR(s->vars.validJets, s->vars.validMuons,      0.3);
        //EventTools::cleanByDR(s->vars.validJets, s->vars.validElectrons,  0.3);
        //EventTools::cleanByDR(s->vars.validBJets, s->vars.validMuons,     0.3);
        //EventTools::cleanByDR(s->vars.validBJets, s->vars.validElectrons, 0.3);

        // Store the valid muons besides the main dimuon candidate
        // AKA initialize validExtraMuons
        //std::cout << "Get extra muons..." << std::endl;
        for(unsigned int m=2; m<s->vars.validMuons.size(); m++)
            s->vars.validExtraMuons.push_back(s->vars.validMuons[m]);

        // catch errors in loading the valid object vectors
        if(s->vars.validMuons.size() <0 || s->vars.validMuonsDecoy.size() <0 || s->vars.validExtraMuons.size() <0
          || s->vars.validElectrons.size() <0 ||  s->vars.validTaus.size() <0 || s->vars.validBJets.size() <0
          || s->vars.validJets.size() < 0)
        {
            std::cout << "!!! ERROR: after loading vectors size < 0..." << std::endl;
            std::cout << "mu size      : " << s->vars.validMuons.size() << std::endl;
            std::cout << "mu decoy size: " << s->vars.validMuonsDecoy.size() << std::endl;
            std::cout << "xmu size     : " << s->vars.validExtraMuons.size() << std::endl;
            std::cout << "e size       : " << s->vars.validElectrons.size() << std::endl;
            std::cout << "t size       : " << s->vars.validTaus.size() << std::endl;
            std::cout << "bjet size    : " << s->vars.validBJets.size() << std::endl;
            std::cout << "jet size     : " << s->vars.validJets.size() << std::endl;
            return 0;
        }
        if(s->vars.validMuons.size() >=11 || s->vars.validMuonsDecoy.size() >=11 || s->vars.validExtraMuons.size() >=11
          || s->vars.validElectrons.size() >=11 ||  s->vars.validTaus.size() >=11 || s->vars.validBJets.size() >=11
          || s->vars.validJets.size() >=11)
        {
            std::cout << "!!! ERROR: after loading vectors size >= 11..." << std::endl;
            std::cout << "mu size      : " << s->vars.validMuons.size() << std::endl;
            std::cout << "mu decoy size: " << s->vars.validMuonsDecoy.size() << std::endl;
            std::cout << "xmu size     : " << s->vars.validExtraMuons.size() << std::endl;
            std::cout << "e size       : " << s->vars.validElectrons.size() << std::endl;
            std::cout << "t size       : " << s->vars.validTaus.size() << std::endl;
            std::cout << "bjet size    : " << s->vars.validBJets.size() << std::endl;
            std::cout << "jet size     : " << s->vars.validJets.size() << std::endl;
            return 0;
        }
        std::pair<int,int> e(s->vars.eventInfo.run, s->vars.eventInfo.event); // create a pair that identifies the event uniquely

        ///////////////////////////////////////////////////////////////////
        // CUTS  ----------------------------------------------------------
        ///////////////////////////////////////////////////////////////////

        // only want to train on events that are in the higgs mass window
        if(!(s->vars.dimuCand.recoCandMassPF > 120 && s->vars.dimuCand.recoCandMass < 130))
        {
            continue;
        }
        // usual cuts
        if(!s->vars.muons.isTightMuon[0] || !s->vars.muons.isTightMuon[1])
        { 
            continue; 
        }
        if(!run1EventSelectionData.evaluate(s->vars) && s->sampleType.Contains("data"))
        { 
            continue; 
        }
        if(!run1EventSelectionMC.evaluate(s->vars) && !s->sampleType.Contains("data"))
        { 
            continue; 
        }
        if(!run1MuonSelection.evaluate(s->vars)) 
        {
            continue; 
        }

        vars["dimu_pt"] = s->vars.dimuCand.recoCandPtPF;
        vars["mu0_pt"] = s->vars.muons.pt[0];
        vars["mu1_pt"] = s->vars.muons.pt[1];
        vars["mu0_eta"] = s->vars.muons.eta[0];
        vars["mu1_eta"] = s->vars.muons.eta[1];

        vars["nValJets"] = s->vars.validJets.size();
        if(s->vars.validJets.size() >= 1) vars["jj_jet0_pt"] = s->vars.validJets[0].Pt();
        if(s->vars.validJets.size() >= 2) vars["jj_jet1_pt"] = s->vars.validJets[1].Pt();

        if(s->vars.validJets.size() >= 1) vars["jj_jet0_eta"] = s->vars.validJets[0].Eta();
        if(s->vars.validJets.size() >= 2) vars["jj_jet1_eta"] = s->vars.validJets[1].Eta();
        if(s->vars.validJets.size() >= 2)
        {
            TLorentzVector dijet = s->vars.validJets[0] + s->vars.validJets[1];
            float dEta = s->vars.validJets[0].Eta() - s->vars.validJets[1].Eta();
            vars["m_jj"] = dijet.M();
            vars["dEta_jj"] = TMath::Abs(dEta);
            float dEtajjmumu = dijet.Eta() - s->vars.dimuCand.recoCandEtaPF;
            vars["dEta_jj_mumu"] = TMath::Abs(dEtajjmumu);
        }

        vars["nExtraMu"] = s->vars.validExtraMuons.size();
        if(s->vars.validExtraMuons.size() >= 1) vars["extra_muon0_pt"] = s->vars.validExtraMuons[0].Pt();
        if(s->vars.validExtraMuons.size() >= 2) vars["extra_muon1_pt"] = s->vars.validExtraMuons[1].Pt();
        if(s->vars.validExtraMuons.size() >= 1) vars["extra_muon0_eta"] = s->vars.validExtraMuons[0].Eta();
        if(s->vars.validExtraMuons.size() >= 2) vars["extra_muon1_eta"] = s->vars.validExtraMuons[1].Eta();
        
        vars["nEle"] = s->vars.validExtraMuons.size();
        if(s->vars.validElectrons.size() >= 1) vars["extra_electron0_pt"] = s->vars.validElectrons[0].Pt();
        if(s->vars.validElectrons.size() >= 2) vars["extra_electron1_pt"] = s->vars.validElectrons[1].Pt();
        if(s->vars.validElectrons.size() >= 1) vars["extra_electron0_eta"] = s->vars.validElectrons[0].Eta();
        if(s->vars.validElectrons.size() >= 2) vars["extra_electron1_eta"] = s->vars.validElectrons[1].Eta();

        vars["nExtraLep"] = s->vars.validElectrons.size() + s->vars.validExtraMuons.size();

        vars["nValBTags"] = s->vars.validBJets.size();
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

        if(s->vars.validBJets.size() > 0)
        {
            TLorentzVector met(s->vars.met.px, s->vars.met.py, 0, s->vars.met.sumEt);
            TLorentzVector bjet = s->vars.validBJets[0];
            TLorentzVector bjet_t(bjet.Px(), bjet.Py(), 0, bjet.Et());
            TLorentzVector bmet_t = met + bjet_t;
            
            vars["mT_b_MET"] = bmet_t.M();
        }

        vars["MET"] = s->vars.met.pt;

        if(false)
          EventTools::outputEvent(s->vars, categorySelection);

        // need to put in is_signal and weight and the extra fancy variables
        vars["is_signal"] = s->sampleType.Contains("signal")?1:0;
        vars["weight"] = s->getScaleFactor(luminosity)*triggerSF*s->getWeight();

        // !!!! output event to file
        file << EventTools::outputMapValuesCSV(vars).Data() << std::endl;

      } // end event loop

    } // end sample loop
    file.close();

    return 0;
}
