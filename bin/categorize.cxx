// Missing HLT trigger info in CMSSW_8_0_X MC so we have to compare Data and MC in a different manner.
// We apply triggers to data but not to MC. Then scale MC for trigger efficiency.

#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "CutSet.h"
#include "Cut.h"
#include "SelectionCuts.h"
#include "CategorySelection.h"
#include "JetSelectionTools.h"
#include "MuSelectionTools.h"
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

    int input = 0;       // the variable to plot, 0 is dimu_mass for instance
    bool rebin = false;  // rebin the histograms so that the ratio plots have small errors
    int binning = 0;     // binning = 1 -> plot dimu_mass from 110 to 160 for limit setting

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) ss >> input;
        if(i==2) ss >> rebin;
        if(i==2) ss >> binning;
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

    //float luminosity = 3990;     // pb-1
    float luminosity = 33598;      // pb-1
    //float luminosity = 27217;      // pb-1
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
    samples["Data"] = datasample;

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

            i.second->lumiWeights = new reweight::LumiReWeighting(i.second->pileupfile.Data(), samples["Data"]->pileupfile.Data(), "pileup", "pileup");
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

    // dimu_mass
    if(input == 0)
    {
        if(binning == 0)
        {
            bins = 150;
            min = 50;
            max = 200;
        }

        else if(binning == 1)
        {
            bins = 50;
            min = 110;
            max = 160;
        }
        else if(binning == 2)
        {
            bins = 100;
            min = 110;
            max = 310;
        }
        else
        {
            bins = 150;
            min = 50;
            max = 200;
        }


        varname = "dimu_mass";
    }

    // dimu_pt 
    if(input == 1)
    {
        bins = 200;
        min = 0;
        max = 100;
        varname = "dimu_pt";
    }

    // mu_pt
    if(input == 2)
    {
        bins = 200;
        min = 0;
        max = 150;
        varname = "mu_pt";
    }
 
    // mu_eta
    if(input == 3)
    {
        bins = 100;
        min = -2.5;
        max = 2.5;
        varname = "mu_eta";
    }

    // NPV
    if(input == 4)
    {
        bins = 50;
        min = 0;
        max = 50;
        varname = "NPV";
    }

    // jet_pt
    if(input == 5)
    {   
        bins = 200;
        min = 0;
        max = 200;
        varname = "jet_pt";
    }   

    // jet_eta 
    if(input == 6)
    {   
        bins = 100;
        min = -5; 
        max = 5;
        varname = "jet_eta";
    }   

    // N_valid_jets
    if(input == 7)
    {   
        bins = 11;
        min = 0; 
        max = 11;
        varname = "N_valid_jets";
    }   

    // m_jj
    if(input == 8)
    {   
        bins = 200;
        min = 0; 
        max = 2000;
        varname = "m_jj";
    }   

    // dEta_jj
    if(input == 9)
    {   
        bins = 100;
        min = -10; 
        max = 10;
        varname = "dEta_jj";
    }   

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

      // Fewer bins for lowstats categories
      //int lowstatsbins = bins/5;
      int lowstatsbins = bins;

      // If we are dealing with NPV or N_valid_jets then don't change the binning
      if(varname.Contains("N")) lowstatsbins = bins;

      // Keep track of which histogram (sample and variable) to fill in the category
      TString hkey;

      // Different categories for the analysis
      for(auto &c : categorySelection.categoryMap)
      {
          //number of bins for the histogram
          int hbins;

          // c.second is the category object, c.first is the category name
          TString hname = varname+"_"+c.first+"_"+s->name;
          hkey = varname+"_"+s->name;

          // The VBF categories have low stats so we use fewer bins
          // Same goes for 01jet categories with dijet variables
          if(c.first.Contains("VBF") || c.first.Contains("GGF") || (c.first.Contains("01_Jet") && varname.Contains("jj"))) hbins = lowstatsbins; 
          else hbins = bins;

          // Set up the histogram for the category and variable to plot
          c.second.histoMap[hkey] = new TH1F(hname, hname, hbins, min, max);
          c.second.histoList->Add(c.second.histoMap[hkey]); // need them ordered by xsec for the stack and ratio plot
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

        if(!s->vars.recoMuons.isTightMuon[0] || !s->vars.recoMuons.isTightMuon[1])
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
            // dimuCand.recoCandMass
            if(varname.EqualTo("dimu_mass")) 
            {
                float varvalue = s->vars.dimuCand.recoCandMassPF;
                // blind the signal region for data but not for MC
                if(!(s->sampleType.Contains("data") && varvalue >= 110 && varvalue < 140))
                    // if the event is in the current category then fill the category's histogram for the given sample and variable
                    if(c.second.inCategory) c.second.histoMap[hkey]->Fill(varvalue, s->getWeight());
                    //std::cout << "    " << c.first << ": " << varvalue;
            }

            if(varname.EqualTo("dimu_pt"))
            {
                // if the event is in the current category then fill the category's histogram for the given sample and variable
                if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.dimuCand.recoCandPtPF, s->getWeight());
            }

            if(varname.EqualTo("mu_pt"))
            {
                if(c.second.inCategory)
                {
                    c.second.histoMap[hkey]->Fill(s->vars.recoMuons.pt[0], s->getWeight());
                    c.second.histoMap[hkey]->Fill(s->vars.recoMuons.pt[1], s->getWeight());
                }
            }

            // recoMu_Eta
            if(varname.EqualTo("mu_eta"))
            {
                if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.recoMuons.eta[0], s->getWeight());
                if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.recoMuons.eta[1], s->getWeight());
            }

            // NPV
            if(varname.EqualTo("NPV"))
            {
                 if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.vertices.nVertices, s->getWeight());
            }

            // jet_pt
            if(varname.EqualTo("jet_pt"))
            {
                if(c.second.inCategory) 
                {
                    for(unsigned int j=0; j<s->vars.validJets.size(); j++)
                        c.second.histoMap[hkey]->Fill(s->vars.validJets[j].Pt(), s->getWeight());
                }
            }

            // jet_eta
            if(varname.EqualTo("jet_eta"))
            {
                if(c.second.inCategory) 
                {
                    for(unsigned int j=0; j<s->vars.validJets.size(); j++)
                        c.second.histoMap[hkey]->Fill(s->vars.validJets[j].Eta(), s->getWeight());
                }
            }

            // N_valid_jets
            if(varname.EqualTo("N_valid_jets"))
            {
                 if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.validJets.size(), s->getWeight());
            }

            // m_jj
            if(varname.EqualTo("m_jj"))
            {
                 if(s->vars.validJets.size() >= 2)
                 {
                     TLorentzVector dijet = s->vars.validJets[0] + s->vars.validJets[1];
                     if(c.second.inCategory) c.second.histoMap[hkey]->Fill(dijet.M(), s->getWeight());
                 }
            }

            // dEta_jj
            if(varname.EqualTo("dEta_jj"))
            {
                 if(s->vars.validJets.size() >= 2)
                 {
                     float dEta = s->vars.validJets[0].Eta() - s->vars.validJets[1].Eta();
                     if(c.second.inCategory) c.second.histoMap[hkey]->Fill(dEta, s->getWeight());
                 }
            }

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

    } // end sample loop

    TList* varstacklist = new TList();  // list to save all of the stacks
    TList* histolist = new TList();     // list to save all of the histos
    TList* netlist = new TList();       // list to save all of the histos

    for(auto &c : categorySelection.categoryMap)
    {
        // Create the stack and ratio plot    
        TString cname = c.first+"_"+varname+"_stack";
        //stackedHistogramsAndRatio(TList* list, TString name, TString title, TString xaxistitle, TString yaxistitle, bool rebin = false, bool fit = true,
                                  //TString ratiotitle = "Data/MC", bool log = true, bool stats = false, int legend = 0);
        TCanvas* stack = dps->stackedHistogramsAndRatio(c.second.histoList, cname, cname, varname, "Num Entries", rebin);
        varstacklist->Add(stack);
        histolist->Add(c.second.histoList);
       
        // we need the data histo, the net signal, and the net bkg for the datacards
        // so we make these histos if we are plotting the dimu_mass spectrum
        if(varname.Contains("dimu_mass"))
        {
            TH1F* hNetSignal = dps->addHists(c.second.signalList, c.first+"_Net_Signal", c.first+"Net_Signal");
            TH1F* hNetBkg    = dps->addHists(c.second.bkgList, c.first+"_Net_Bkg", c.first+"_Net_Bkg");
            TH1F* hNetData   = dps->addHists(c.second.dataList, c.first+"_Net_Data", c.first+"_Net_Data");

            netlist->Add(hNetSignal);
            netlist->Add(hNetBkg);
            netlist->Add(hNetData);
        }

        stack->SaveAs("imgs/"+cname+".png");
    }
    std::cout << std::endl;

    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile("rootfiles/validate_"+varname+Form("_%d_%d", (int)min, (int)max)+
                                "_x69p2_8_0_X_MC_categories_"+Form("%d",(int)luminosity)+Form("_rebin%d.root", (int)rebin), "RECREATE");
    TDirectory* stacks = savefile->mkdir("stacks");
    TDirectory* histos = savefile->mkdir("histos");
    TDirectory* net_histos = savefile->mkdir("net_histos");

    // save the different histos and stacks in the appropriate directories in the tfile
    stacks->cd();
    varstacklist->Write();
    histos->cd();
    histolist->Write();
    net_histos->cd();
    netlist->Write();
    savefile->Close();

    return 0;
}
