// output amount of signal and background in each category
// also output the significance in each category
// store info into a CSV file, save histos into a root file as well

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

    // Use this to plot some things if we wish
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

  
    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    float luminosity = 36000;      // pb-1
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

        // Don't worry about pileup reweighting
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
    JetSelectionTools jetSelectionTools;
    CategorySelectionRun1 categorySelection;
    Run1MuonSelectionCuts run1MuonSelection;
    Run1EventSelectionCuts80X run1EventSelectionData(true);
    Run1EventSelectionCuts80X run1EventSelectionMC;

    TString varname;
    int bins;
    float min;
    float max;

    bins = 150;
    min = 50;
    max = 200;
    varname = "dimu_mass";

    std::cout << std::endl;
    std::cout << "======== Plot Configs ========" << std::endl;
    std::cout << "var         : " << varname << std::endl;
    std::cout << "min         : " << min << std::endl;
    std::cout << "max         : " << max << std::endl;
    std::cout << "bins        : " << bins << std::endl;
    std::cout << std::endl;

    // reduce the number of events you run over in case you want to debug or some such thing
    float reductionFactor = 1;

    for(auto &s : samplevec)
    {
      // Output some info about the current file
      std::cout << std::endl;
      std::cout << "  /// Looping over " << s->name << std::endl;

      ///////////////////////////////////////////////////////////////////
      // HISTOGRAMS TO FILL ---------------------------------------------
      ///////////////////////////////////////////////////////////////////

      // Fewer bins for lowstats categories if necessary
      int lowstatsbins = bins/5;

      // If we are dealing with NPV or nValJets then don't change the binning
      if(varname.Contains("N")) lowstatsbins = bins;

      // Keep track of which histogram (sample and variable) to fill in the category
      TString hkey;

      // Different categories for the analysis
      for(auto &c : categorySelection.categoryMap)
      {
          //number of bins for the histogram
          int hbins;

          // c.second is the category object, c.first is the category name
          TString hname = c.first+"_"+s->name;
          hkey = s->name;

          // The VBF categories have low stats so we use fewer bins
          // Same goes for 01jet categories with dijet variables
          if(c.first.Contains("VBF") || c.first.Contains("GGF") || (c.first.Contains("01_Jet") && varname.Contains("jj"))) hbins = lowstatsbins; 
          else hbins = bins;

          // Set up the histogram for the category and variable to plot
          c.second.histoMap[hkey] = new TH1D(hname, hname, hbins, min, max);
          c.second.histoMap[hkey]->GetXaxis()->SetTitle(varname);
          c.second.histoList->Add(c.second.histoMap[hkey]);                                        // need them ordered by xsec for the stack and ratio plot
          if(s->sampleType.Contains("data")) c.second.dataList->Add(c.second.histoMap[hkey]);      // data histo
          if(s->sampleType.Contains("signal")) c.second.signalList->Add(c.second.histoMap[hkey]);  // signal histos
          if(s->sampleType.Contains("background")) c.second.bkgList->Add(c.second.histoMap[hkey]); // bkg histos
      }

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

        if(!(s->vars.dimuCand.recoCandMassPF > 120 && s->vars.dimuCand.recoCandMassPF < 130))
        { 
            continue; 
        }
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

        // Figure out which category the event belongs to
        categorySelection.evaluate(s->vars);

        // Look at each category
        for(auto &c : categorySelection.categoryMap)
        {
           // fill the dimu mass histograms
           float varvalue = s->vars.dimuCand.recoCandMassPF;
           // blind the signal region for data but not for MC
           if(!(s->sampleType.Contains("data") && varvalue >= 110 && varvalue < 140))
               // if the event is in the current category then fill the category's histogram for the given sample and variable
               if(c.second.inCategory) c.second.histoMap[hkey]->Fill(varvalue, s->getWeight());
               //std::cout << "    " << c.first << ": " << varvalue;

        } // end category loop

        if(false)
          // ouput pt, mass info etc for the event
          EventTools::outputEvent(s->vars, categorySelection);

        // Reset the flags in preparation for the next event
        categorySelection.reset();

      } // end event loop

      // Scale according to luminosity and sample xsec now that the histograms are done being filled for that sample
      for(auto &c : categorySelection.categoryMap)
      {
          c.second.histoMap[hkey]->Scale(s->getScaleFactor(luminosity));
          if(!s->sampleType.Contains("data")) c.second.histoMap[hkey]->Scale(triggerSF);
      }
      std::cout << std::endl;

    } // end sample loop

    TList* signallist = new TList();         // list to save all of the signal histos
    TList* bglist = new TList();             // list to save all of the background histos
    TList* netlist = new TList();            // list to save all of the net histos

    for(auto &c : categorySelection.categoryMap)
    {
        // some categories are intermediate and we don't want to save the plots for those
        if(c.second.hide) continue;

        // lists will contain signal, bg, and data histos for every category
        signallist->Add(c.second.signalList);
        bglist->Add(c.second.bkgList);
       
        // we need the data histo, the net signal, and the net bkg dimu mass histos for the datacards
        // so we make these histos. Might as well make them for every variable, not just dimu_mass.
        TH1D* hNetSignal = dps->addHists(c.second.signalList, c.first+"_Net_Signal", c.first+"_Net_Signal");
        TH1D* hNetBkg    = dps->addHists(c.second.bkgList,    c.first+"_Net_Bkg",    c.first+"_Net_Bkg");

        netlist->Add(hNetSignal);
        netlist->Add(hNetBkg);

        // ////////////////////////////////////////////////////////////////////////////
        // ========= Total Counts =====================================================
        // ////////////////////////////////////////////////////////////////////////////

        double nsignal = hNetSignal->Integral(0, hNetSignal->GetSize());
        double nbackground = hNetBkg->Integral(0, hNetBkg->GetSize());

        std::cout << std::endl;
        std::cout << "=========== " << c.first << " Total Counts ============" << std::endl;
        std::cout << "Signal:     " << nsignal << std::endl;
        std::cout << "Background: " << nbackground << std::endl;
        std::cout << "Total:      " << nsignal + nbackground << std::endl;
        std::cout << std::endl;

        // ////////////////////////////////////////////////////////////////////////////
        // ========= Significance =====================================================
        // ////////////////////////////////////////////////////////////////////////////

        AsimovSignificance asimov0(0);
        AsimovSignificance asimov1(1);
        PoissonSignificance poisson0(0);
        PoissonSignificance poisson1(1);

        // unc = 0.06 for 12596 background events, scales like 1/sqrt(N)
        double asimovZ0 = asimov0.significance(nsignal, nbackground);
        double asimovZ1 = asimov1.significance(nsignal, nbackground);
        double rootbZ0 = poisson0.significance(nsignal, nbackground);
        double rootbZ1 = poisson1.significance(nsignal, nbackground);

        std::cout << std::endl;
        std::cout << "=========== " << c.first << " Significance ============" << std::endl;
        std::cout << "Asimov w/o bgfit:        " << asimovZ0 << std::endl;
        std::cout << "Asimov w/ bgfit:         " << asimovZ1 << std::endl;
        std::cout << "s/sqrt(s + b):           " << rootbZ0 << std::endl;
        std::cout << "s/sqrt(s + b + varb):    " << rootbZ1 << std::endl;
        std::cout << std::endl;
    }

    std::cout << std::endl;

    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile("rootfiles/counts_"+varname+Form("_%d_%d", (int)min, (int)max)+
                                "_8_0_X_MC_run1categories_"+Form("%d",(int)luminosity)+".root", "RECREATE");

    TDirectory* signal_histos = savefile->mkdir("signal_histos");
    TDirectory* bg_histos     = savefile->mkdir("bg_histos");
    TDirectory* net_histos    = savefile->mkdir("net_histos");

    // save the different histos and stacks in the appropriate directories in the tfile
    signal_histos->cd();
    signallist->Write();

    bg_histos->cd();
    bglist->Write();

    net_histos->cd();
    netlist->Write();

    EventTools::outputCategoryCountsToFile(categorySelection, "csv/counts_run1.csv");

    savefile->Close();

    return 0;
}
