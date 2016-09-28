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
    //float luminosity = 12900;      // pb-1 (ICHEP)
    //float luminosity = 15900;      // pb-1   (2016-08-03)
    float triggerSF = 0.913;       // no HLT trigger info available for the samples so we scale for the trigger efficiency instead
    float signalSF = 100;

    // ================================================================
    // Data -----------------------------------------------------------
    // ================================================================
/* 

    TString datafilename = 
    TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data/25ns/golden/CMSSW_8_0_X/old/stage_1_singleMuon_Run2016B_ALL.root");
    //TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data/25ns/golden/CMSSW_8_0_X/stage_1_singleMuon_Run2016BCD_ICHEP_ALL.root");
    //TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data/25ns/golden/CMSSW_8_0_X/stage_1_singleMuon_Run2016BCDE_ALL.root");

    Sample* datasample = new Sample(datafilename, "Data", "data");
    datasample->lumi = luminosity;
    datasample->xsec = 9999;
    datasample->pileupfile = "pu_reweight_trees/8_0_X/PU_2016B_xsec69mb_CMSSW_8_0_X.root";
    //datasample->pileupfile = "pu_reweight_trees/8_0_X/PU_2016BCD_ICHEP_xsec69p2mb_CMSSW_8_0_X.root";
    //datasample->pileupfile = "pu_reweight_trees/8_0_X/PU_2016BCDE_xsec69p2mb_CMSSW_8_0_X.root";
    samples["Data"] = datasample;
*/
    // ================================================================
    // DYJetsToLL -----------------------------------------------------
    // ================================================================

    TString dyfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/bg/dy/CMSSW_8_0_X/stage_1_dy_jetsToLL_ALL.root");
    samples["DYJetsToLL"] = new Sample(dyfilename, "DYJetsToLL", "background");
    samples["DYJetsToLL"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_DYJetsToLL.root"; //nPU
    samples["DYJetsToLL"]->xsec = 6025.2; // pb

    ///////////////////////////////////////////////////////////////////
    // PREPROCESSING---------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // Loop through all of the samples to do some pre-processing
    std::cout << std::endl;
    std::cout << "======== Preprocess the samples... " << std::endl;
    std::cout << std::endl;
    
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
/*
        if(!i.second->sampleType.Contains("data"))
        {
            // Pileup reweighting
            std::cout << "    +++ PU Reweighting " << i.second->name << "..."  << std::endl;
            std::cout << std::endl;

            i.second->lumiWeights = new reweight::LumiReWeighting(i.second->pileupfile.Data(), samples["Data"]->pileupfile.Data(), "pileup", "pileup");
            std::cout << "        " << i.first << "->lumiWeights: " << i.second->lumiWeights << std::endl;
            std::cout << std::endl;
        }
*/
        samplevec.push_back(i.second);
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 
    

    ///////////////////////////////////////////////////////////////////
    // Cut and Categorize ---------------------------------------------
    ///////////////////////////////////////////////////////////////////
    
    // Objects to help with the cuts and selections
    JetSelectionTools jetSelectionTools;
    FEWZCompareCuts FEWZselection;
    CategorySelectionFEWZ categorySelection;

    TString varname;
    int nbins, wbins;
    float nmin, wmin;
    float nmax, wmax;

    // dimu_mass
    if(input == 0)
    {
        nbins = 50;
        nmin = 110;
        nmax = 160;

        wbins = 100;
        wmin = 110;
        wmax = 310;

        varname = "dimu_mass";
    }

    // dimu_pt 
    if(input == 1)
    {
        nbins = 200;
        nmin = 0;
        nmax = 100;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "dimu_pt";
    }

    // recoMu_pt
    if(input == 2)
    {
        nbins = 200;
        nmin = 0;
        nmax = 150;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "recoMu_pt";
    }
 
    // recoMu_eta
    if(input == 3)
    {
        nbins = 100;
        nmin = -2.5;
        nmax = 2.5;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "recoMu_eta";
    }

    // NPV
    if(input == 4)
    {
        nbins = 50;
        nmin = 0;
        nmax = 50;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "NPV";
    }

    // jet_pt
    if(input == 5)
    {   
        nbins = 200;
        nmin = 0;
        nmax = 200;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "jet_pt";
    }   

    // jet_eta 
    if(input == 6)
    {   
        nbins = 100;
        nmin = -5; 
        nmax = 5;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "jet_eta";
    }   

    // N_valid_jets
    if(input == 7)
    {   
        nbins = 11;
        nmin = 0; 
        nmax = 11;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "N_valid_jets";
    }   

    // m_jj
    if(input == 8)
    {   
        nbins = 200;
        nmin = 0; 
        nmax = 2000;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "m_jj";
    }   

    // dEta_jj
    if(input == 9)
    {   
        nbins = 100;
        nmin = -10; 
        nmax = 10;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "dEta_jj";
    }   

    std::cout << std::endl;
    std::cout << "======== Plot Configs ========" << std::endl;
    std::cout << "var          : " << varname << std::endl;
    std::cout << "nmin         : " << nmin << std::endl;
    std::cout << "nmax         : " << nmax << std::endl;
    std::cout << "nbins        : " << nbins << std::endl;
    std::cout << std::endl;
    std::cout << "wmin         : " << wmin << std::endl;
    std::cout << "wmax         : " << wmax << std::endl;
    std::cout << "wbins        : " << wbins << std::endl;
    std::cout << std::endl;

    // Not sure how to deal with the scaling correctly when using a subset of events
    float reductionFactor = 1;

    for(auto &s : samplevec)
    {
      // Output some info about the current file
      std::cout << std::endl;
      std::cout << "  /// Looping over " << s->name << std::endl;

      ///////////////////////////////////////////////////////////////////
      // HISTOGRAMS TO FILL ---------------------------------------------
      ///////////////////////////////////////////////////////////////////

      // Keep track of which histogram (sample and variable) to fill in the category
      TString hkey;

      // We usually fill the same plots for each category, just with different events
      // so we loop through each category and fill the appropriate events
      for(auto &c : categorySelection.categoryMap)
      {
          //number of bins for the histogram, initialize to narrow category binning
          //overwrite later if it's a wide category
          int hbins = nbins;
          double hmin = nmin;
          double hmax = nmax;

          // c.second is the category object, c.first is the category name
          TString hname = varname+"_"+c.first+"_"+s->name;
          hkey = varname+"_"+s->name;

          // The VBF categories have low stats so we use fewer bins
          // Same goes for 01jet categories with dijet variables
          if(c.first.Contains("Wide")) 
          {
              hbins = wbins; 
              hmin = wmin;
              hmax = wmax;
          }

          // Set up the histogram for the category and variable to plot
          c.second.histoMap[hkey] = new TH1F(hname, hname, hbins, hmin, hmax);
          c.second.histoList->Add(c.second.histoMap[hkey]); // need them ordered by xsec for the stack and ratio plot
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
        
        // don't apply iso cuts
        FEWZselection.cutset.cuts[7].on = false;
        FEWZselection.cutset.cuts[8].on = false;

        if(!s->vars.reco1.isTightMuon || !s->vars.reco2.isTightMuon)
        {
            continue; 
        }
        if(!FEWZselection.evaluate(s->vars)) 
        {
            continue; 
        }

        // Figure out which category the event belongs to
        categorySelection.evaluate(s->vars);

        // Look at each category
        for(auto &c : categorySelection.categoryMap)
        {
            // recoCandMass
            if(varname.Contains("dimu_mass")) 
            {
                float varvalue = s->vars.recoCandMassPF;
                // blind the signal region for data but not for MC
                if(!(s->sampleType.Contains("data") && varvalue >= 110 && varvalue < 140))
                    // if the event is in the current category then fill the category's histogram for the given sample and variable
                    if(c.second.inCategory) c.second.histoMap[hkey]->Fill(varvalue, s->getWeight());
                    //std::cout << "    " << c.first << ": " << varvalue;
            }

            if(varname.Contains("dimu_pt"))
            {
                // if the event is in the current category then fill the category's histogram for the given sample and variable
                if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.recoCandPtPF, s->getWeight());
            }

            if(varname.Contains("recoMu_pt"))
            {
                if(c.second.inCategory)
                {
                    c.second.histoMap[hkey]->Fill(s->vars.reco1.pt, s->getWeight());
                    c.second.histoMap[hkey]->Fill(s->vars.reco2.pt, s->getWeight());
                }
            }

            // recoMu_Eta
            if(varname.Contains("recoMu_eta"))
            {
                if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.reco1.eta, s->getWeight());
                if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.reco2.eta, s->getWeight());
            }

            // NPV
            if(varname.Contains("NPV"))
            {
                 if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.vertices.nVertices, s->getWeight());
            }

            // jet_pt
            if(varname.Contains("jet_pt"))
            {
                if(c.second.inCategory) 
                {
                    for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                        c.second.histoMap[hkey]->Fill(s->vars.validJets[i].Pt(), s->getWeight());
                }
            }

            // jet_eta
            if(varname.Contains("jet_Eta"))
            {
                if(c.second.inCategory) 
                {
                    for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                        c.second.histoMap[hkey]->Fill(s->vars.validJets[i].Eta(), s->getWeight());
                }
            }

            // N_valid_jets
            if(varname.Contains("N_valid_jets"))
            {
                 if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.validJets.size(), s->getWeight());
            }

            // m_jj
            if(varname.Contains("m_jj"))
            {
                 if(s->vars.validJets.size() >= 2)
                 {
                     TLorentzVector dijet = s->vars.validJets[0] + s->vars.validJets[1];
                     if(c.second.inCategory) c.second.histoMap[hkey]->Fill(dijet.M(), s->getWeight());
                 }
            }

            // dEta_jj
            if(varname.Contains("dEta_jj"))
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
          outputEvent(s->vars, categorySelection);

        // Reset the flags in preparation for the next event
        categorySelection.reset();

      } // end event loop

      // Scale according to luminosity and sample xsec now that the histograms are done being filled for that sample
      for(auto &c : categorySelection.categoryMap)
          c.second.histoMap[hkey]->Scale(s->getScaleFactor(luminosity));

    } // end sample loop

    //TList* varstacklist = new TList();  // list to save all of the stacks
    TList* histolist = new TList();     // list to save all of the histos

    // Not ready for stack & ratio yet, need to make those plots with corresponding FEWZ plots instead of the usual
    // Get the histos from each category and add to net histo list for saving
    // For each category we want to compare the dimu_mass for DY with the dimu mass for FEWZ
    for(auto &c : categorySelection.categoryMap)
    {
    //    // Name the canvas 
    //    TString cname = c.first+"_"+varname+"_stack";
    //    // Get the dy histogram that we just made earlier
    //    TCanvas* stack = dps->stackedHistogramsAndRatio(c.second.histoList, cname, cname, varname, "Num Entries", "FEWZ/DY_MC");
    //    varstacklist->Add(stack);
        histolist->Add(c.second.histoList);
    //    stack->SaveAs("../python/fewz/img/"+cname+".png");
    }
   
    //std::cout << std::endl;

    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile("rootfiles/validate_"+varname+"_DY-FEWZ_MC_categories_"+Form("%d",(int)luminosity)+".root", "RECREATE");
    //TDirectory* stacks = savefile->mkdir("stacks");
    TDirectory* histos = savefile->mkdir("histos");

    // save the different histos and stacks in the appropriate directories in the tfile
    //stacks->cd();
    //varstacklist->Write();
    histos->cd();
    histolist->Write();
    savefile->Close();

    return 0;
}
