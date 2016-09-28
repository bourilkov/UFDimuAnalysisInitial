// Missing HLT trigger info in CMSSW_8_0_X MC so we have to compare Data and MC in a different manner.
// We apply triggers to data but not to MC. Then scale MC for trigger efficiency.


#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "CutSet.h"
#include "Cut.h"
#include "SelectionCuts.h"
#include "CategorySelection_v2.h"
#include "JetSelectionTools.h"

#include "EventTools_v2.cxx"
#include "PUTools.cxx"
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
    int input = 0;
    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) ss >> input;
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

    float luminosity = 3990;       // pb-1
    float triggerSF = 0.913;       // no HLT trigger info available for the samples so we scale for the trigger efficiency instead
    float signalSF = 100;

    // ================================================================
    // Data -----------------------------------------------------------
    // ================================================================
 

    TString datafilename = 
    TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data/25ns/golden/CMSSW_8_0_X/stage_1_singleMuon_Run2016B_ALL.root");

    Sample* datasample = new Sample(datafilename, "Data", "data");
    datasample->lumi = luminosity;
    datasample->xsec = 9999;
    //datasample->pileupfile = "/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data_from_json/25ns/golden/pileup/old/PUCalib_Golden_71mb.root";
    datasample->pileupfile = "pu_reweight_trees/8_0_X/PU_2016B_xsec69mb_CMSSW_8_0_X.root";
    samples["Data"] = datasample;

/*
    // ================================================================
    // DYJetsToLL -----------------------------------------------------
    // ================================================================

    TString dyfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/bg/dy/CMSSW_8_0_X/stage_1_dy_jetsToLL_ALL.root");
    samples["DYJetsToLL"] = new Sample(dyfilename, "DYJetsToLL", "background");
    samples["DYJetsToLL"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_DYJetsToLL.root"; //nPU
    samples["DYJetsToLL"]->xsec = 6025.2; // pb

    // ================================================================
    // TTJets ---------------------------------------------------------
    // ================================================================

    TString ttbarfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/bg/ttbar/CMSSW_8_0_X/stage_1_ttJets2_ALL.root");
    samples["TTJets"] = new Sample(ttbarfilename, "TTJets", "background");
    samples["TTJets"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_TTJets.root"; //nPU
    samples["TTJets"]->xsec = 831.76; // pb

    // ================================================================
    // VBF ---------------------------------------------------------
    // ================================================================

    TString vbffilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/signal/CMSSW_8_0_X/stage_1_vbf_HToMuMu_ALL.root");
    samples["VBF"] = new Sample(vbffilename, "VBF", "signal");
    samples["VBF"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_VBF.root"; //nPU
    samples["VBF"]->xsec = 3.727*0.00022; // pb

    // ================================================================
    // GGF ---------------------------------------------------------
    // ================================================================

    TString ggfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/signal/CMSSW_8_0_X/stage_1_gg_HToMuMu2_ALL.root");
    samples["GGF"] = new Sample(ggfilename, "GGF", "signal");
    samples["GGF"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_GGF.root"; //nPU
    samples["GGF"]->xsec = 43.62*0.00022; // pb
*/

    ///////////////////////////////////////////////////////////////////
    // PREPROCESSING---------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // Loop through all of the samples to do some pre-processing
    std::cout << std::endl;
    std::cout << "======== Preprocess the samples... " << std::endl;
    std::cout << std::endl;

    //makePUHistos(samples);
    
    for(auto const &i : samples)
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

            i.second->lumiWeights = new reweight::LumiReWeighting(i.second->pileupfile.Data(), samples["Data"]->pileupfile.Data(), "pileup", "pileup");
            std::cout << "        " << i.first << "->lumiWeights: " << i.second->lumiWeights << std::endl;
            std::cout << std::endl;
        }
        samplevec.push_back(i.second);
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 
    
    TList* histolist = new TList();   // the list of objects to save
    TList* graphlist = new TList();   // the list of objects to save

    ///////////////////////////////////////////////////////////////////
    // Cut and Categorize ---------------------------------------------
    ///////////////////////////////////////////////////////////////////
    
    // Objects to help with the cuts and selections
    JetSelectionTools jetSelectionTools;
    CategorySelection categorySelection;
    Run1MuonSelectionCuts run1MuonSelection;
    Run1EventSelectionCuts80X run1EventSelectionData(true); //true sets the slightly different data cuts
    Run1EventSelectionCuts80X run1EventSelectionMC;

    ///////////////////////////////////////////////////////////////////
    // Plot Settings Depending on Input--------------------------------
    ///////////////////////////////////////////////////////////////////
    
    TString xname;
    Float_t fitsig, massmin, massmax, xmin, xmax;
    Int_t massbins, xbins;

    massmin = 86.2;
    massmax = 96.2;
    massbins = 50;
    xbins = 25;
    fitsig = 1;

    if(input == 0)
    {
        xname = "phi_plus";
        xmin = -3.14;
        xmax = 3.14;
    }

    if(input == 1)
    {
        xname = "phi_minus";
        xmin = -3.14;
        xmax = 3.14;
    }

    if(input == 2)
    {
        xname = "eta_plus";
        xmin = -2.4;
        xmax = 2.4;
    }

    if(input == 3)
    {
        xname = "eta_minus";
        xmin = -2.4;
        xmax = 2.4;
    }


    // Not sure how to deal with the scaling correctly when using a subset of events
    float reductionFactor = 1;

    for(auto const &s : samplevec)
    {
      // Output some info about the current file
      std::cout << std::endl;
      std::cout << "  /// Looping over " << s->name << std::endl;

      ZCalibration* zcal = new ZCalibration(xname, fitsig, massmin, massmax, massbins, xmin, xmax, xbins);

      ///////////////////////////////////////////////////////////////////
      // HISTOGRAMS TO FILL ---------------------------------------------
      ///////////////////////////////////////////////////////////////////

      for(unsigned int i=0; i<s->N/reductionFactor; i++)
      {

        ///////////////////////////////////////////////////////////////////
        // GET INFORMATION ------------------------------------------------
        ///////////////////////////////////////////////////////////////////

        s->getEntry(i); 
        s->vars.validJets = std::vector<TLorentzVector>();
        jetSelectionTools.getValidJetsdR(s->vars, s->vars.validJets);
        std::pair<int,int> e(s->vars.eventInfo.run, s->vars.eventInfo.event); // create a pair that identifies the event uniquely

        ///////////////////////////////////////////////////////////////////
        // CUTS  ----------------------------------------------------------
        ///////////////////////////////////////////////////////////////////

        if(!s->vars.reco1.isTightMuon || !s->vars.reco2.isTightMuon)
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

        //phi plus
        if(s->vars.reco1.charge == 1 && input == 0) zcal->fill(s->vars.reco1.phi, s->vars.recoCandMassPF);
        if(s->vars.reco2.charge == 1 && input == 0) zcal->fill(s->vars.reco2.phi, s->vars.recoCandMassPF);

        //phi minus
        if(s->vars.reco1.charge == -1 && input == 1) zcal->fill(s->vars.reco1.phi, s->vars.recoCandMassPF);
        if(s->vars.reco2.charge == -1 && input == 1) zcal->fill(s->vars.reco2.phi, s->vars.recoCandMassPF);

        //eta plus
        if(s->vars.reco1.charge == 1 && input == 2) zcal->fill(s->vars.reco1.eta, s->vars.recoCandMassPF);
        if(s->vars.reco2.charge == 1 && input == 2) zcal->fill(s->vars.reco2.eta, s->vars.recoCandMassPF);

        //eta minus
        if(s->vars.reco1.charge == -1 && input == 3) zcal->fill(s->vars.reco1.eta, s->vars.recoCandMassPF);
        if(s->vars.reco2.charge == -1 && input == 3) zcal->fill(s->vars.reco2.eta, s->vars.recoCandMassPF);

        if(false)
          // ouput pt, mass info etc
          outputEvent(s->vars, categorySelection);

        // Reset the flags in preparation for the next event
        categorySelection.reset();
      }

      for(unsigned int i=0; i<zcal->histos.size(); i++)
          histolist->Add(zcal->histos[i]);

      graphlist->Add(zcal->plot());
    }

    // Create the stack and ratio plot    
    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile("rootfiles/"+xname+"_data_8_0_X.root", "RECREATE");
    TDirectory* graphs = savefile->mkdir("graphs");
    TDirectory* hists = savefile->mkdir("hists");

    // save the different histos in the appropriate directories in the tfile
    hists->cd();
    histolist->Write();

    graphs->cd();
    graphlist->Write();

    savefile->Close();

    return 0;
}
