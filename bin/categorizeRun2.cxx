// plot different variables in the run 2 categories.
// may be used to look for discrepancies between data and mc or to make dimu_mass plots for limit setting.
// outputs mc stacks with data overlayed and a ratio plot underneath.
// also saves the histos needed to make the mc stack, data.
// Also saves net BKG histo and net signal histo for limit setting.

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

    int input = 0;       // the variable to plot, 0 is dimu_mass for instance
    bool rebin = false;  // rebin the histograms so that the ratio plots have small errors
    int binning = 0;     // binning = 1 -> plot dimu_mass from 110 to 160 for limit setting
    unsigned int nPartitions = 1; // break the samples up into nPartitions to run in parallel
    unsigned int partition = 0;   // which partition to make the plots for

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) ss >> input;
        if(i==2) ss >> rebin;
        if(i==3) ss >> binning;
        if(i==4) ss >> nPartitions;
        if(i==5) ss >> partition;
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

    // N_valid_muons
    if(input == 10)
    {   
        bins = 11;
        min = 0; 
        max = 11;
        varname = "N_valid_muons";
    }   

    // N_valid_extra_muons
    if(input == 11)
    {   
        bins = 11;
        min = 0; 
        max = 11;
        varname = "N_valid_extra_muons";
    }   

    // extra_muon_pt
    if(input == 12)
    {
        bins = 200;
        min = 0;
        max = 150;
        varname = "extra_muon_pt";
    }
 
    // extra_muon_eta
    if(input == 13)
    {
        bins = 100;
        min = -3;
        max = 3;
        varname = "extra_muon_eta";
    }

    // N_valid_electrons
    if(input == 14)
    {   
        bins = 11;
        min = 0; 
        max = 11;
        varname = "N_valid_electrons";
    }   

    // electron_pt
    if(input == 15)
    {
        bins = 200;
        min = 0;
        max = 150;
        varname = "electron_pt";
    }
 
    // electron_eta
    if(input == 16)
    {
        bins = 100;
        min = -3;
        max = 3;
        varname = "electron_eta";
    }

    // N_valid_extra_leptons
    if(input == 17)
    {   
        bins = 11;
        min = 0; 
        max = 11;
        varname = "N_valid_extra_leptons";
    }   

    // N_valid_bjets
    if(input == 18)
    {   
        bins = 11;
        min = 0; 
        max = 11;
        varname = "N_valid_bjets";
    }   

    // bjet_pt
    if(input == 19)
    {   
        bins = 200;
        min = 0;
        max = 200;
        varname = "bjet_pt";
    }   

    // bjet_eta 
    if(input == 20)
    {   
        bins = 100;
        min = -5; 
        max = 5;
        varname = "bjet_eta";
    }   

    // m_bb
    if(input == 21)
    {   
        bins = 200;
        min = 0; 
        max = 2000;
        varname = "m_bb";
    }   

    // mT_b_MET
    if(input == 22)
    {   
        bins = 200;
        min = 0; 
        max = 2000;
        varname = "mT_b_MET";
    }   

    // MET
    if(input == 23)
    {
        bins = 200;
        min = 0;
        max = 150;
        varname = "MET";
    }
 
    // dEta_jj_mumu
    if(input == 24)
    {
        bins = 100;
        min = -10; 
        max = 10;
        varname = "dEta_jj_mumu";
    }
 
    std::cout << std::endl;
    std::cout << "======== Plot Configs ========" << std::endl;
    std::cout << "var         : " << varname << std::endl;
    std::cout << "rebin       : " << rebin << std::endl;
    std::cout << "binning     : " << binning << std::endl;
    std::cout << "nPartitions : " << nPartitions << std::endl;
    std::cout << "partition   : " << partition << std::endl;
    std::cout << std::endl;
    std::cout << "min         : " << min << std::endl;
    std::cout << "max         : " << max << std::endl;
    std::cout << "bins        : " << bins << std::endl;
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

      ///////////////////////////////////////////////////////////////////
      // HISTOGRAMS TO FILL ---------------------------------------------
      ///////////////////////////////////////////////////////////////////

      // Fewer bins for lowstats categories if necessary
      int lowstatsbins = bins;
      if(!rebin)  lowstatsbins = bins/5;

      // If we are dealing with NPV or N_valid_jets then don't change the binning
      if(varname.Contains("N")) lowstatsbins = bins;

      // Keep track of which histogram (sample and variable) to fill in the category
      TString hkey;

      // Different categories for the analysis
      for(auto &c : categorySelection.categoryMap)
      {
          //number of bins for the histogram
          int hbins = bins;

          // c.second is the category object, c.first is the category name
          TString hname = c.first+"_"+s->name;
          hkey = s->name;

          // Not sure what the low stats categories are for run2 yet
          //if(c.first.Contains("VBF") || c.first.Contains("GGF") || (c.first.Contains("01_Jet") && varname.Contains("jj"))) hbins = lowstatsbins; 
          //else hbins = bins;

          // Set up the histogram for the category and variable to plot
          c.second.histoMap[hkey] = new TH1D(hname, hname, hbins, min, max);
          c.second.histoMap[hkey]->GetXaxis()->SetTitle(varname);
          c.second.histoList->Add(c.second.histoMap[hkey]); // need them ordered by xsec for the stack and ratio plot
          if(s->sampleType.Contains("data")) c.second.dataList->Add(c.second.histoMap[hkey]);      // data histo
          if(s->sampleType.Contains("signal")) c.second.signalList->Add(c.second.histoMap[hkey]);  // signal histos
          if(s->sampleType.Contains("background")) c.second.bkgList->Add(c.second.histoMap[hkey]); // bkg histos
      }

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

        //std::cout << "Apply tight mu..." << std::endl;
        if(!s->vars.recoMuons.isTightMuon[0] || !s->vars.recoMuons.isTightMuon[1])
        { 
            continue; 
        }
        //std::cout << "Apply run1 event selection..." << std::endl;
        if(!run1EventSelectionData.evaluate(s->vars) && s->sampleType.Contains("data"))
        { 
            continue; 
        }
        if(!run1EventSelectionMC.evaluate(s->vars) && !s->sampleType.Contains("data"))
        { 
            continue; 
        }
        //std::cout << "Apply run1 muon selection..." << std::endl;
        if(!run1MuonSelection.evaluate(s->vars)) 
        {
            continue; 
        }

        // Figure out which category the event belongs to
        //std::cout << "Evaluate categories..." << std::endl;
        categorySelection.evaluate(s->vars);

        // Look at each category
        //std::cout << "Fill each category histo..." << std::endl;
        for(auto &c : categorySelection.categoryMap)
        {
            //std::cout << "in category loop..." << std::endl;
            // dimuCand.recoCandMass
            if(varname.EqualTo("dimu_mass")) 
            {
                //std::cout << "in dimu_mass fill..." << std::endl;
                float varvalue = s->vars.dimuCand.recoCandMassPF;
                // blind the signal region for data but not for MC
                //if(c.second.inCategory) std::cout << "    " << c.first << ": " << varvalue << std::endl;
                if(!(s->sampleType.Contains("data") && varvalue >= 110 && varvalue < 140))
                {
                    // if the event is in the current category then fill the category's histogram for the given sample and variable
                    if(c.second.inCategory) c.second.histoMap[hkey]->Fill(varvalue, s->getWeight());
                }
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
            // N_valid_muons
            if(varname.EqualTo("N_valid_muons"))
            {
                 if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.validMuons.size(), s->getWeight());
            }

            // N_valid_extra_muons
            if(varname.EqualTo("N_valid_extra_muons"))
            {
                 if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.validExtraMuons.size(), s->getWeight());
            }
            
            // extra_muon_pt
            if(varname.EqualTo("extra_muon_pt"))
            {
                if(c.second.inCategory) 
                {
                    for(unsigned int j=0; j<s->vars.validExtraMuons.size(); j++)
                        c.second.histoMap[hkey]->Fill(s->vars.validExtraMuons[j].Pt(), s->getWeight());
                }
            }
 
            // extra_muon_eta
            if(varname.EqualTo("extra_muon_eta"))
            {
                if(c.second.inCategory) 
                {
                    for(unsigned int j=0; j<s->vars.validExtraMuons.size(); j++)
                        c.second.histoMap[hkey]->Fill(s->vars.validExtraMuons[j].Eta(), s->getWeight());
                }
            }
 
            // N_valid_electrons
            if(varname.EqualTo("N_valid_electrons"))
            {
                 if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.validElectrons.size(), s->getWeight());
            }
            
            // electron_pt
            if(varname.EqualTo("electron_pt"))
            {
                if(c.second.inCategory) 
                {
                    for(unsigned int j=0; j<s->vars.validElectrons.size(); j++)
                        c.second.histoMap[hkey]->Fill(s->vars.validElectrons[j].Pt(), s->getWeight());
                }
            }
 
            // electron_eta
            if(varname.EqualTo("electron_eta"))
            {
                if(c.second.inCategory) 
                {
                    for(unsigned int j=0; j<s->vars.validElectrons.size(); j++)
                        c.second.histoMap[hkey]->Fill(s->vars.validElectrons[j].Eta(), s->getWeight());
                }
            }

            // N_valid_extra_leptons
            if(varname.EqualTo("N_valid_extra_leptons"))
            {
                 if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.validElectrons.size() + s->vars.validExtraMuons.size(), s->getWeight());
            }

            // N_valid_bjets
            if(varname.EqualTo("N_valid_bjets"))
            {
                 if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.validBJets.size(), s->getWeight());
            }

            // bjet_pt
            if(varname.EqualTo("bjet_pt"))
            {
                if(c.second.inCategory) 
                {
                    for(unsigned int j=0; j<s->vars.validBJets.size(); j++)
                        c.second.histoMap[hkey]->Fill(s->vars.validBJets[j].Pt(), s->getWeight());
                }
            }
            
            // bjet_eta 
            if(varname.EqualTo("bjet_eta"))
            {
                if(c.second.inCategory) 
                {
                    for(unsigned int j=0; j<s->vars.validBJets.size(); j++)
                        c.second.histoMap[hkey]->Fill(s->vars.validBJets[j].Eta(), s->getWeight());
                }
            }

            // m_bb
            if(varname.EqualTo("m_bb"))
            {
                 if(s->vars.validBJets.size() >= 2)
                 {
                     TLorentzVector dijet = s->vars.validBJets[0] + s->vars.validBJets[1];
                     if(c.second.inCategory) c.second.histoMap[hkey]->Fill(dijet.M(), s->getWeight());
                 }
            }

            // mT_b_MET
            if(varname.EqualTo("mT_b_MET"))
            {
                 if(s->vars.validBJets.size() > 0)
                 {
                     TLorentzVector met(s->vars.met.px, s->vars.met.py, 0, s->vars.met.sumEt);
                     TLorentzVector bjet = s->vars.validBJets[0];
                     TLorentzVector bjet_t(bjet.Px(), bjet.Py(), 0, bjet.Et());
                     TLorentzVector bmet_t = met + bjet_t;
                     
                     if(c.second.inCategory) c.second.histoMap[hkey]->Fill(bmet_t.M(), s->getWeight());
                 }
            }

            // MET
            if(varname.EqualTo("MET"))
            {
                 if(c.second.inCategory) c.second.histoMap[hkey]->Fill(s->vars.met.pt, s->getWeight());
            }

            // dEta_jj_mumu
            if(varname.EqualTo("dEta_jj_mumu"))
            {
                 if(s->vars.validJets.size() >= 2)
                 {
                     TLorentzVector dijet = s->vars.validJets[0] + s->vars.validJets[1];
                     float dEta = dijet.Eta() - s->vars.dimuCand.recoCandEtaPF;
                     if(c.second.inCategory) c.second.histoMap[hkey]->Fill(dEta, s->getWeight());
                 }
            }
        } // end category loop

        if(false)
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

    TList* varstacklist = new TList();   // list to save all of the stacks
    TList* signallist = new TList();     // list to save all of the signal histos
    TList* bglist = new TList();         // list to save all of the background histos
    TList* datalist = new TList();       // list to save all of the data histos
    TList* netlist = new TList();        // list to save all of the net histos

    for(auto &c : categorySelection.categoryMap)
    {
        // some categories are intermediate and we don't want to save the plots for those
        if(c.second.hide) continue;

        // Create the stack and ratio plot    
        TString cname = c.first+"_stack";
        //stackedHistogramsAndRatio(TList* list, TString name, TString title, TString xaxistitle, TString yaxistitle, bool rebin = false, bool fit = true,
                                  //TString ratiotitle = "Data/MC", bool log = true, bool stats = false, int legend = 0);
        // stack signal, bkg, and data
        TCanvas* stack = dps->stackedHistogramsAndRatio(c.second.histoList, cname, cname, varname, "Num Entries", rebin);
        varstacklist->Add(stack);
       
        // lists will contain signal, bg, and data histos for every category
        signallist->Add(c.second.signalList);
        bglist->Add(c.second.bkgList);
        datalist->Add(c.second.dataList);

        // we need the data histo, the net signal, and the net bkg dimu mass histos for the datacards
        // so we make these histos. Might as well make them for every variable, not just dimu_mass.
        TH1D* hNetSignal = dps->addHists(c.second.signalList, c.first+"_Net_Signal", c.first+"_Net_Signal");
        TH1D* hNetBkg    = dps->addHists(c.second.bkgList,    c.first+"_Net_Bkg",    c.first+"_Net_Bkg");
        TH1D* hNetData   = dps->addHists(c.second.dataList,   c.first+"_Net_Data",   c.first+"_Net_Data");

        netlist->Add(hNetSignal);
        netlist->Add(hNetBkg);
        netlist->Add(hNetData);

        // fails here with seg-fault if the TCanvas 'stack' was created with bad info, aka 0 integral histograms for numerator or denominator
        stack->SaveAs("imgs/"+cname+".png");
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile("rootfiles/validate_"+varname+Form("_%d_%d", (int)min, (int)max)+
                                "_x69p2_8_0_X_MC_lotsOfCategoriesRun2_"+Form("%d",(int)luminosity)+
                                Form("_rebin%d_partition%d-%d.root", (int)rebin, partition, nPartitions), "RECREATE");

    TDirectory* stacks        = savefile->mkdir("stacks");
    TDirectory* signal_histos = savefile->mkdir("signal_histos");
    TDirectory* bg_histos     = savefile->mkdir("bg_histos");
    TDirectory* data_histos   = savefile->mkdir("data_histos");
    TDirectory* net_histos    = savefile->mkdir("net_histos");

    // save the different histos and stacks in the appropriate directories in the tfile
    stacks->cd();
    varstacklist->Write();

    signal_histos->cd();
    signallist->Write();

    bg_histos->cd();
    bglist->Write();

    data_histos->cd();
    datalist->Write();

    net_histos->cd();
    netlist->Write();

    savefile->Close();

    return 0;
}
