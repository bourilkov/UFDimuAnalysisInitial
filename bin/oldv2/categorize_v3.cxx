// plot different variables in the run 1 or run 2 categories.
// may be used to look for discrepancies between data and mc or to make dimu_mass plots for limit setting.
// outputs mc stacks with data overlayed and a ratio plot underneath.
// also saves the histos needed to make the mc stack, data.
// Also saves net BKG histo and net signal histo for limit setting, saves individuals as well

// Missing HLT trigger info in CMSSW_8_0_X MC so we have to compare Data and MC in a different manner.
// We apply triggers to data but not to MC. Then scale MC for trigger efficiency.

#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "EventSelection.h"
#include "MuonSelection.h"
#include "CategorySelection.h"
#include "JetCollectionCleaner.h"
#include "MuonCollectionCleaner.h"
#include "EleCollectionCleaner.h"

#include "EventTools.h"
#include "PUTools.h"

#include "SignificanceMetrics.hxx"
#include "SampleDatabase.cxx"
#include "ThreadPool.hxx"

#include <sstream>
#include <map>
#include <vector>
#include <utility>

#include "TLorentzVector.h"
#include "TSystem.h"
#include "TBranchElement.h"
//#include "TStreamerInfo.h"

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

TList* getSortedMC(TList* list, std::vector<Sample*>& sampleVec)
{
// Sorts according to the xsec of the sample

    // map sample name to sorted vector location
    std::map<TString, unsigned int> xsecMap;
    for(unsigned int i=0; i< sampleVec.size(); i++)
        xsecMap[sampleVec[i]->name] = i;

    // strip data and get a sorted vector of mc samples
    std::vector<TH1D*> mcVec(sampleVec.size());
    unsigned int ndata = 0;
    for(unsigned int i=0; i<list->GetSize(); i++)
    {
        TH1D* hist = (TH1D*)list->At(i); 
        TString name = hist->GetName();
        int uscore = name.Last('_');
        TString sampleName = name(uscore+1,name.Length());

        if(sampleName.Contains("Run") || sampleName.Contains("Data"))
            ndata++;
        else
            mcVec[xsecMap[sampleName]] = hist;    
    }

    // put the sorted mc into a list
    TList* sortedMCList = new TList();
    for(unsigned int i=0; i<mcVec.size()-ndata; i++)
    {
        sortedMCList->Add(mcVec[i]); 
    }

    return sortedMCList;
}


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

// set bins, min, max, and the var to plot based upon terminal varNumber: varNumber and binning
void initPlotSettings(int varNumber, int binning, int& bins, float& min, float& max, TString& varname)
{
    // dimu_mass
    if(varNumber == 0)
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
    if(varNumber == 1)
    {
        bins = 200;
        min = 0;
        max = 100;
        varname = "dimu_pt";
    }

    // mu_pt
    if(varNumber == 2)
    {
        bins = 200;
        min = 0;
        max = 150;
        varname = "mu_pt";
    }
 
    // mu_eta
    if(varNumber == 3)
    {
        bins = 100;
        min = -2.5;
        max = 2.5;
        varname = "mu_eta";
    }

    // NPV
    if(varNumber == 4)
    {
        bins = 50;
        min = 0;
        max = 50;
        varname = "NPV";
    }

    // jet_pt
    if(varNumber == 5)
    {   
        bins = 200;
        min = 0;
        max = 200;
        varname = "jet_pt";
    }   

    // jet_eta 
    if(varNumber == 6)
    {   
        bins = 100;
        min = -5; 
        max = 5;
        varname = "jet_eta";
    }   

    // nValJets
    if(varNumber == 7)
    {   
        bins = 11;
        min = 0; 
        max = 11;
        varname = "nValJets";
    }   

    // m_jj
    if(varNumber == 8)
    {   
        bins = 200;
        min = 0; 
        max = 2000;
        varname = "m_jj";
    }   

    // dEta_jj
    if(varNumber == 9)
    {   
        bins = 100;
        min = -10; 
        max = 10;
        varname = "dEta_jj";
    }   
    // N_valid_muons
    if(varNumber == 10) 
    {   
        bins = 11; 
        min = 0;  
        max = 11; 
        varname = "N_valid_muons";
    }   

    // nExtraMu
    if(varNumber == 11) 
    {   
        bins = 11; 
        min = 0;  
        max = 11; 
        varname = "nExtraMu";
    }   

    // extra_muon_pt
    if(varNumber == 12) 
    {   
        bins = 200;
        min = 0;
        max = 150;
        varname = "extra_muon_pt";
    }   
 
    // extra_muon_eta
    if(varNumber == 13)
    {
        bins = 100;
        min = -3;
        max = 3;
        varname = "extra_muon_eta";
    }

    // nEle
    if(varNumber == 14)
    {
        bins = 11;
        min = 0;
        max = 11;
        varname = "nEle";
    }
    // electron_pt
    if(varNumber == 15)
    {
        bins = 200;
        min = 0;
        max = 150;
        varname = "electron_pt";
    }

    // electron_eta
    if(varNumber == 16)
    {
        bins = 100;
        min = -3;
        max = 3;
        varname = "electron_eta";
    }

    // nExtraLep
    if(varNumber == 17)
    {
        bins = 11;
        min = 0;
        max = 11;
        varname = "nExtraLep";
    }

    // nValBTags
    if(varNumber == 18)
    {
        bins = 11;
        min = 0;
        max = 11;
        varname = "nValBTags";
    }

    // bjet_pt
    if(varNumber == 19)
    {
        bins = 200;
        min = 0;
        max = 200;
        varname = "bjet_pt";
    }
    // bjet_eta 
    if(varNumber == 20)
    {
        bins = 100;
        min = -5;
        max = 5;
        varname = "bjet_eta";
    }

    // m_bb
    if(varNumber == 21)
    {
        bins = 200;
        min = 0;
        max = 2000;
        varname = "m_bb";
    }

    // mT_b_MET
    if(varNumber == 22)
    {
        bins = 200;
        min = 0;
        max = 2000;
        varname = "mT_b_MET";
    }

    // MET
    if(varNumber == 23)
    {
        bins = 200;
        min = 0;
        max = 150;
        varname = "MET";
    }

    // dEta_jj_mumu
    if(varNumber == 24)
    {
        bins = 100;
        min = -10;
        max = 10;
        varname = "dEta_jj_mumu";
    }
}


//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    gROOT->SetBatch();

    // save the errors for the histogram correctly so they depend upon 
    // the number used to fill originally rather than the scaling
    TH1::SetDefaultSumw2();

    int whichCategories = 1;  // run2categories = 1, run2categories = 2
    int varNumber = 0;        // the variable to plot, 0 is dimu_mass for instance
    int nthreads = 5;         // number of threads to use in parallelization
    bool rebin = true;        // rebin the histograms so that the ratio plots have small errors
    int binning = 0;          // binning = 1 -> plot dimu_mass from 110 to 160 for limit setting
    bool fitratio = 0;        // fit the ratio plot (data/mc) under the stack w/ a straight line

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) ss >> whichCategories;
        if(i==2) ss >> varNumber;
        if(i==3) ss >> nthreads;
        if(i==4) ss >> rebin;
        if(i==5) ss >> binning;
        if(i==6) ss >> fitratio;
    }   
    std::cout << "@@@ nCPUs Available: " << getNumCPUs() << std::endl;
    std::cout << "@@@ nCPUs used     : " << nthreads << std::endl;

    // Not sure that we need a map if we have a vector
    // Should use this as the main database and choose from it to make the vector
    std::map<TString, Sample*> samples;

    // Second container so that we can have a copy sorted by cross section.
    std::vector<Sample*> samplevec;

    // Use this to plot some things if we wish
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

    float luminosity = 36814;      // pb-1
    float triggerSF = 0.913;       // no HLT trigger info available for the samples so we scale for the trigger efficiency instead
    float signalSF = 100;          // not using this at the moment, but scale the signal samples to see them better in the plots if you want
    float reductionFactor = 100;   // reduce the number of events you run over in case you want to debug or some such thing

    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // gather samples map from SamplesDatabase.cxx
    getSamples(luminosity, samples);

    ///////////////////////////////////////////////////////////////////
    // PREPROCESSING: Get PU weights ----------------------------------
    ///////////////////////////////////////////////////////////////////

    // Loop through all of the samples to do some pre-processing
    // Add to the vector we loop over to categorize
    
    //std::cout << std::endl;
    //std::cout << "======== Preprocess the samples... " << std::endl;
    //std::cout << std::endl;

    //makePUHistos(samples);
    
    for(auto &i : samples)
    {
        // Output some info about the current file
        std::cout << "  /// Using the following samples " << i.second->name << std::endl;
        std::cout << std::endl;
        std::cout << "    sample name:       " << i.second->name << std::endl;
        std::cout << "    sample file:       " << i.second->filename << std::endl;
        std::cout << "    pileup file:       " << i.second->pileupfile << std::endl;
        std::cout << "    nOriginal:         " << i.second->nOriginal << std::endl;
        std::cout << "    N:                 " << i.second->N << std::endl;
        std::cout << "    nOriginalWeighted: " << i.second->nOriginalWeighted << std::endl;
        std::cout << std::endl;

        samplevec.push_back(i.second);

        // Pileup Reweighting now available in the TTree info, no need to calculate
        //if(!i.second->sampleType.EqualTo("data"))
        //{
        //    // Pileup reweighting
        //    std::cout << "    +++ PU Reweighting " << i.second->name << "..."  << std::endl;
        //    std::cout << std::endl;

        //    // Every data sample has the pileup hist for the entire lumi
        //    i.second->lumiWeights = new reweight::LumiReWeighting(i.second->pileupfile.Data(), samples["RunB"]->pileupfile.Data(), "pileup", "pileup");
        //    std::cout << "        " << i.first << "->lumiWeights: " << i.second->lumiWeights << std::endl;
        //    std::cout << std::endl;
        //}
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 

    ///////////////////////////////////////////////////////////////////
    // Get Histo Information from Input -------------------------------
    ///////////////////////////////////////////////////////////////////

    TString varname;
    int bins;
    float min;
    float max;

    // set nbins, min, max, and var to plot based upon terminal varNumbers: varNumber and binning
    initPlotSettings(varNumber, binning, bins, min, max, varname);

    std::cout << std::endl;
    std::cout << "======== Plot Configs ========" << std::endl;
    std::cout << "categories  : " << whichCategories << std::endl;
    std::cout << "var         : " << varname << std::endl;
    std::cout << "min         : " << min << std::endl;
    std::cout << "max         : " << max << std::endl;
    std::cout << "bins        : " << bins << std::endl;
    std::cout << "rebin       : " << rebin << std::endl;
    std::cout << std::endl;

    std::mutex mtx;

    ///////////////////////////////////////////////////////////////////
    // Define Task for Parallelization -------------------------------
    ///////////////////////////////////////////////////////////////////

    auto makeHistoForSample = [&mtx, varname, bins, min, max, rebin, whichCategories, triggerSF, luminosity, reductionFactor](Sample* s)
    {
      // Output some info about the current file
      std::cout << std::endl;
      std::cout << "  /// Processing " << s->name << std::endl;

      ///////////////////////////////////////////////////////////////////
      // INIT Cuts and Categories ---------------------------------------
      ///////////////////////////////////////////////////////////////////
      
      // Objects to help with the cuts and selections
      JetCollectionCleaner      jetCollectionCleaner;
      MuonCollectionCleaner     muonCollectionCleaner;
      EleCollectionCleaner      eleCollectionCleaner;

      Run2MuonSelectionCuts     run2MuonSelection;
      Run2EventSelectionCuts80X run2EventSelectionData(true);
      Run2EventSelectionCuts80X run2EventSelectionMC;

      std::cout << "  /// Setting up categorizer " << s->name << std::endl;

      Categorizer* categorySelection = 0;
      if(whichCategories == 1) categorySelection = new CategorySelectionRun1();
      else categorySelection = new LotsOfCategoriesRun2();

      bool isData = s->sampleType.EqualTo("data");

      ///////////////////////////////////////////////////////////////////
      // INIT HISTOGRAMS TO FILL ----------------------------------------
      ///////////////////////////////////////////////////////////////////

      // Fewer bins for lowstats categories if necessary
      int lowstatsbins = bins/5;

      // If we are dealing with NPV or nValJets then don't change the binning
      if(varname.Contains("N")) lowstatsbins = bins;

      // Keep track of which histogram to fill in the category
      TString hkey = s->name;

      std::cout << "  /// Initializing histograms " << s->name << std::endl;
      // Different categories for the analysis
      for(auto &c : categorySelection->categoryMap)
      {
          //number of bins for the histogram
          int hbins;

          // c.second is the category object, c.first is the category name
          TString hname = c.first+"_"+s->name;

          // The VBF categories have low stats so we use fewer bins
          // Same goes for 01jet categories with dijet variables
          if(c.first.Contains("VBF") || c.first.Contains("GGF") || (c.first.Contains("01_Jet") && varname.Contains("jj"))) hbins = lowstatsbins; 
          else hbins = bins;

          // Set up the histogram for the category and variable to plot
          c.second.histoMap[hkey] = new TH1D(hname, hname, hbins, min, max);
          c.second.histoMap[hkey]->GetXaxis()->SetTitle(varname);
          c.second.histoList->Add(c.second.histoMap[hkey]);                                      // need them ordered by xsec for the stack and ratio plot
          if(s->sampleType.EqualTo("data")) c.second.dataList->Add(c.second.histoMap[hkey]);     // data histo
          if(s->sampleType.EqualTo("signal")) c.second.signalList->Add(c.second.histoMap[hkey]); // signal histos
          if(s->sampleType.EqualTo("background")) c.second.bkgList->Add(c.second.histoMap[hkey]);// bkg histos

      }

      std::cout << "  /// Getting branches " << s->name << std::endl;
      // Link only to the branches we need to save a lot of time
      // run 1 category info 
      TBranch* muPairsBranch  = s->tree->GetBranch("pairs");
      TBranch* muonsBranch      = s->tree->GetBranch("muons");
      TBranch* jetsBranch           = s->tree->GetBranch("jets");
      TBranch* mhtBranch            = s->tree->GetBranch("mht");
      TBranch* eventInfoBranch      = s->tree->GetBranch("event");
      TBranch* nVerticesBranch      = s->tree->GetBranch("nVertices");

      std::cout << "  /// Setting branch addresses " << s->name << std::endl;
      muPairsBranch->SetAddress(&s->vars.muPairs);
      muonsBranch->SetAddress(&s->vars.muons);
      jetsBranch->SetAddress(&s->vars.jets);
      mhtBranch->SetAddress(&s->vars.mht);
      eventInfoBranch->SetAddress(&s->vars.eventInfo);
      nVerticesBranch->SetAddress(&s->vars.nVertices);

      // extra branches needed for run 2 categories
      TBranch* electronsBranch  = 0;

      if(whichCategories == 2)
      {
          std::cout << "  /// Run2 category branches " << s->name << std::endl;
          electronsBranch = s->tree->GetBranch("eles");
          electronsBranch->SetAddress(&s->vars.electrons);
      }

      // extra branches needed for MC samples
      TBranch* nPUBranch = 0;
      TBranch* genWeightBranch = 0;
      TBranch* puWeightBranch = 0;
      TBranch* effWeightBranch = 0;

      if(!isData)
      { 
          std::cout << "  /// MC branches " << s->name << std::endl;
          nPUBranch       = s->tree->GetBranch("nPU");
          genWeightBranch = s->tree->GetBranch("GEN_wgt");
          puWeightBranch  = s->tree->GetBranch("PU_wgt");
          effWeightBranch = s->tree->GetBranch("IsoMu_eff_3");

          genWeightBranch->SetAddress(&s->vars.gen_wgt);
          nPUBranch->SetAddress(&s->vars.nPU);
          puWeightBranch->SetAddress(&s->vars.pu_wgt);
          effWeightBranch->SetAddress(&s->vars.eff_wgt);
      }

      ///////////////////////////////////////////////////////////////////
      // LOOP OVER EVENTS -----------------------------------------------
      ///////////////////////////////////////////////////////////////////

      std::cout << "  /// Loop over events " << s->name << std::endl;
      for(unsigned int i=0; i<s->N/reductionFactor; i++)
      {

        // only load essential information for the first set of cuts 
        muPairsBranch->GetEntry(i);
        muonsBranch->GetEntry(i);

        // insert loop over pairs
        if(s->vars.muPairs->size() >= 1) s->vars.dimuCand = &s->vars.muPairs->at(0);
        else continue;

        ///////////////////////////////////////////////////////////////////
        // CUTS  ----------------------------------------------------------
        ///////////////////////////////////////////////////////////////////

        if(!run2EventSelectionData.evaluate(s->vars) && isData)
        { 
            continue; 
        }
        if(!run2EventSelectionMC.evaluate(s->vars) && !isData)
        { 
            continue; 
        }
        if(!s->vars.muons->at(s->vars.dimuCand->iMu1).isTightID || !s->vars.muons->at(s->vars.dimuCand->iMu2).isTightID)
        { 
            continue; 
        }
        if(!run2MuonSelection.evaluate(s->vars)) 
        {
            continue; 
        }

        // Load the rest of the information needed for run2 categories
        jetsBranch->GetEntry(i);
        mhtBranch->GetEntry(i);
        nVerticesBranch->GetEntry(i);
        //eventInfoBranch->GetEntry(i);

        if(!isData)
        {
            genWeightBranch->GetEntry(i);
            nPUBranch->GetEntry(i);
            puWeightBranch->GetEntry(i);
            effWeightBranch->GetEntry(i);
        }

        if(whichCategories == 1) 
        {
            s->vars.validJets.clear();
            jetCollectionCleaner.getValidJetsdR(s->vars, s->vars.validJets);
        }
        //std::pair<int,int> e(s->vars.eventInfo.run, s->vars.eventInfo.event); // create a pair that identifies the event uniquely

        if(whichCategories == 2)
        {
            // load extra branches needed by run 2 categories
            //electronsBranch->GetEntry(i);

            // clear vectors for the valid collections
            s->vars.validMuons.clear();
            s->vars.validExtraMuons.clear();
            s->vars.validElectrons.clear();
            s->vars.validJets.clear();
            s->vars.validBJets.clear();

            // load valid collections from s->vars raw collections
            jetCollectionCleaner.getValidJets(s->vars, s->vars.validJets);
            jetCollectionCleaner.getValidBJets(s->vars, s->vars.validBJets);
            muonCollectionCleaner.getValidMuons(s->vars, s->vars.validMuons);
            eleCollectionCleaner.getValidElectrons(s->vars, s->vars.validElectrons);
          
            //for(unsigned int m=0; m<s->vars.validMuons.size(); m++)
            //    s->vars.validExtraMuons.push_back(s->vars.validMuons[m]);
        }

        // Figure out which category the event belongs to
        categorySelection->evaluate(s->vars);

        // Look at each category
        for(auto &c : categorySelection->categoryMap)
        {
            // skip categories
            if(c.second.hide) continue;
            if(!c.second.inCategory) continue;

            // dimuCand->recoCandMass
            if(varname.EqualTo("dimu_mass")) 
            {
                float varvalue = s->vars.dimuCand->mass_PF;
                // blind the signal region for data but not for MC
                if(!(isData && varvalue >= 110 && varvalue < 140))
                {
                    // if the event is in the current category then fill the category's histogram for the given sample and variable
                    c.second.histoMap[hkey]->Fill(varvalue, s->getWeight());
                    //std::cout << "    " << c.first << ": " << varvalue;
                }
               continue;
            }

            if(varname.EqualTo("dimu_pt"))
            {
                // if the event is in the current category then fill the category's histogram for the given sample and variable
                c.second.histoMap[hkey]->Fill(s->vars.dimuCand->pt_PF, s->getWeight());
                continue;
            }

            if(varname.EqualTo("mu_pt"))
            {
                c.second.histoMap[hkey]->Fill(s->vars.muons->at(0).pt, s->getWeight());
                c.second.histoMap[hkey]->Fill(s->vars.muons->at(1).pt, s->getWeight());
                continue;
            }

            // recoMu_Eta
            if(varname.EqualTo("mu_eta"))
            {
                c.second.histoMap[hkey]->Fill(s->vars.muons->at(0).eta, s->getWeight());
                c.second.histoMap[hkey]->Fill(s->vars.muons->at(1).eta, s->getWeight());
                continue;
            }

            // NPV
            if(varname.EqualTo("NPV"))
            {
                 c.second.histoMap[hkey]->Fill(s->vars.nVertices, s->getWeight());
                 continue;
            }

            // jet_pt
            if(varname.EqualTo("jet_pt"))
            {
                 for(unsigned int j=0; j<s->vars.validJets.size(); j++)
                     c.second.histoMap[hkey]->Fill(s->vars.validJets[j].Pt(), s->getWeight());
                 continue;
            }

            // jet_eta
            if(varname.EqualTo("jet_eta"))
            {
                 for(unsigned int j=0; j<s->vars.validJets.size(); j++)
                     c.second.histoMap[hkey]->Fill(s->vars.validJets[j].Eta(), s->getWeight());
                 continue;
            }

            // nValJets
            if(varname.EqualTo("nValJets"))
            {
                 c.second.histoMap[hkey]->Fill(s->vars.validJets.size(), s->getWeight());
                 continue;
            }

            // m_jj
            if(varname.EqualTo("m_jj"))
            {
                 if(s->vars.validJets.size() >= 2)
                 {
                     TLorentzVector dijet = s->vars.validJets[0] + s->vars.validJets[1];
                     c.second.histoMap[hkey]->Fill(dijet.M(), s->getWeight());
                 }
                 continue;
            }

            // dEta_jj
            if(varname.EqualTo("dEta_jj"))
            {
                 if(s->vars.validJets.size() >= 2)
                 {
                     float dEta = s->vars.validJets[0].Eta() - s->vars.validJets[1].Eta();
                     c.second.histoMap[hkey]->Fill(dEta, s->getWeight());
                 }
                 continue;
            }

            // N_valid_muons
            if(varname.EqualTo("N_valid_muons"))
            {
                 c.second.histoMap[hkey]->Fill(s->vars.validMuons.size(), s->getWeight());
                 continue;
            }

            // nExtraMu
            if(varname.EqualTo("nExtraMu"))
            {
                 c.second.histoMap[hkey]->Fill(s->vars.validExtraMuons.size(), s->getWeight());
                 continue;
            }

            // extra_muon_pt
            if(varname.EqualTo("extra_muon_pt"))
            {
                for(unsigned int j=0; j<s->vars.validExtraMuons.size(); j++)
                    c.second.histoMap[hkey]->Fill(s->vars.validExtraMuons[j].Pt(), s->getWeight());
                continue;
            }

            // extra_muon_eta
            if(varname.EqualTo("extra_muon_eta"))
            {
                for(unsigned int j=0; j<s->vars.validExtraMuons.size(); j++)
                    c.second.histoMap[hkey]->Fill(s->vars.validExtraMuons[j].Eta(), s->getWeight());
                continue;
            }

            // nEle
            if(varname.EqualTo("nEle"))
            {
                 c.second.histoMap[hkey]->Fill(s->vars.validElectrons.size(), s->getWeight());
                 continue;
            }

            // electron_pt
            if(varname.EqualTo("electron_pt"))
            {
                for(unsigned int j=0; j<s->vars.validElectrons.size(); j++)
                    c.second.histoMap[hkey]->Fill(s->vars.validElectrons[j].Pt(), s->getWeight());
                continue;
            }
            // electron_eta
            if(varname.EqualTo("electron_eta"))
            {
                for(unsigned int j=0; j<s->vars.validElectrons.size(); j++)
                    c.second.histoMap[hkey]->Fill(s->vars.validElectrons[j].Eta(), s->getWeight());
                continue;
            }

            // nExtraLep
            if(varname.EqualTo("nExtraLep"))
            {
                 c.second.histoMap[hkey]->Fill(s->vars.validElectrons.size() + s->vars.validExtraMuons.size(), s->getWeight());
                 continue;
            }

            // nValBTags
            if(varname.EqualTo("nValBTags"))
            {
                 c.second.histoMap[hkey]->Fill(s->vars.validBJets.size(), s->getWeight());
                 continue;
            }

            // bjet_pt
            if(varname.EqualTo("bjet_pt"))
            {
                for(unsigned int j=0; j<s->vars.validBJets.size(); j++)
                    c.second.histoMap[hkey]->Fill(s->vars.validBJets[j].Pt(), s->getWeight());
                continue;
            }

            // bjet_eta 
            if(varname.EqualTo("bjet_eta"))
            {
                for(unsigned int j=0; j<s->vars.validBJets.size(); j++)
                    c.second.histoMap[hkey]->Fill(s->vars.validBJets[j].Eta(), s->getWeight());
                continue;
            }
            // m_bb
            if(varname.EqualTo("m_bb"))
            {
                 if(s->vars.validBJets.size() >= 2)
                 {
                     TLorentzVector dijet = s->vars.validBJets[0] + s->vars.validBJets[1];
                     c.second.histoMap[hkey]->Fill(dijet.M(), s->getWeight());
                 }
                 continue;
            }

            // mT_b_MET
            //if(varname.EqualTo("mT_b_MET"))
            //{
            //     if(s->vars.validBJets.size() > 0)
            //     {
            //         TLorentzVector mht(s->vars.mht.px, s->vars.mht.py, 0, s->vars.met.sumEt);
            //         TLorentzVector bjet = s->vars.validBJets[0];
            //         TLorentzVector bjet_t(bjet.Px(), bjet.Py(), 0, bjet.Et());
            //         TLorentzVector bmet_t = met + bjet_t;

            //         c.second.histoMap[hkey]->Fill(bmet_t.M(), s->getWeight());
            //     }
            //     continue;
            //}

            // MHT
            if(varname.EqualTo("MHT"))
            {
                c.second.histoMap[hkey]->Fill(s->vars.mht->pt, s->getWeight());
            }

            // dEta_jj_mumu
            if(varname.EqualTo("dEta_jj_mumu"))
            {
                 if(s->vars.validJets.size() >= 2)
                 {
                     TLorentzVector dijet = s->vars.validJets[0] + s->vars.validJets[1];
                     float dEta = dijet.Eta() - s->vars.dimuCand->eta;
                     c.second.histoMap[hkey]->Fill(dEta, s->getWeight());
                 }
                 continue;
            }

        } // end category loop

        if(false)
          // ouput pt, mass info etc for the event
          EventTools::outputEvent(s->vars, *categorySelection);

        // Reset the flags in preparation for the next event
        categorySelection->reset();

      } // end event loop

      // Scale according to luminosity and sample xsec now that the histograms are done being filled for that sample
      for(auto &c : categorySelection->categoryMap)
      {
          c.second.histoMap[hkey]->Scale(s->getScaleFactor(luminosity));
          //if(!isData) c.second.histoMap[hkey]->Scale(triggerSF);
      }

      return categorySelection;

    }; // done defining makeHistoForSample

   ///////////////////////////////////////////////////////////////////
   // SAMPLE PARALLELIZATION------ ----------------------------------
   ///////////////////////////////////////////////////////////////////

    ThreadPool pool(nthreads);
    std::vector< std::future<Categorizer*> > results;

    TStopwatch timerWatch;
    timerWatch.Start();

    for(auto &s : samplevec)
        results.push_back(pool.enqueue(makeHistoForSample, s));

   ///////////////////////////////////////////////////////////////////
   // Gather all the Histos into one Categorizer----------------------
   ///////////////////////////////////////////////////////////////////

    Categorizer* cAll = 0;
    if(whichCategories == 1) cAll = new CategorySelectionRun1(); 
    else cAll = new LotsOfCategoriesRun2();

    // get histos from all categorizers and put them into one
    for(auto && categorizer: results)  // loop through each Categorizer object, one per sample
    {
        for(auto& category: categorizer.get()->categoryMap) // loop through each category
        {
            // category.first is the category name, category.second is the Category object
            if(category.second.hide) continue;
            for(auto& h: category.second.histoMap) // loop through each histogram in the category
            {
                // h.first is the sample name, category.histoMap<samplename, TH1D*>
                Sample* s = samples[h.first];

                cAll->categoryMap[category.first].histoMap[h.first] = h.second;
                cAll->categoryMap[category.first].histoList->Add(h.second);

                if(s->sampleType.EqualTo("signal"))          cAll->categoryMap[category.first].signalList->Add(h.second);
                else if(s->sampleType.EqualTo("background")) cAll->categoryMap[category.first].bkgList->Add(h.second);
                else                                         cAll->categoryMap[category.first].dataList->Add(h.second);
            }
        }
    }

   ///////////////////////////////////////////////////////////////////
   // Save All of the Histos-----------------------------------------
   ///////////////////////////////////////////////////////////////////

    TList* varstacklist = new TList();   // list to save all of the stacks
    TList* signallist = new TList();     // list to save all of the signal histos
    TList* bglist = new TList();         // list to save all of the background histos
    TList* datalist = new TList();       // list to save all of the data histos
    TList* netlist = new TList();        // list to save all of the net histos

    for(auto &c : cAll->categoryMap)
    {
        // some categories are intermediate and we don't want to save the plots for those
        if(c.second.hide) continue;

        // we need the data histo, the net signal, and the net bkg dimu mass histos for the datacards
        TH1D* hNetSignal = dps->addHists(c.second.signalList, c.first+"_Net_Signal", c.first+"_Net_Signal");
        TH1D* hNetBkg    = dps->addHists(c.second.bkgList,    c.first+"_Net_Bkg",    c.first+"_Net_Bkg");
        TH1D* hNetData   = dps->addHists(c.second.dataList,   c.first+"_Net_Data",   c.first+"_Net_Data");

        netlist->Add(hNetSignal);
        netlist->Add(hNetBkg);
        netlist->Add(hNetData);

        TList* sortedlist = getSortedMC(c.second.histoList, samplevec);
        sortedlist->Add(hNetData);

        // Create the stack and ratio plot    
        TString cname = c.first+"_stack";
        //stackedHistogramsAndRatio(TList* list, TString name, TString title, TString xaxistitle, TString yaxistitle, bool rebin = false, bool fit = true,
                                  //TString ratiotitle = "Data/MC", bool log = true, bool stats = false, int legend = 0);
        // stack signal, bkg, and data
        TCanvas* stack = dps->stackedHistogramsAndRatio(sortedlist, cname, cname, varname, "Num Entries", rebin, fitratio);
        varstacklist->Add(stack);

        // lists will contain signal, bg, and data histos for every category
        signallist->Add(c.second.signalList);
        bglist->Add(c.second.bkgList);
        datalist->Add(c.second.dataList);
       
        stack->SaveAs("imgs/"+cname+".png");
    }
    std::cout << std::endl;
    TString savename = Form("validate_%s_%d_%d_x69p2_8_0_X_MC_run%dcategories_%d.root", varname.Data(), (int)min, (int)max, whichCategories, (int)luminosity);

    std::cout << "  /// Saving plots to " << savename << " ..." << std::endl;
    std::cout << std::endl;

    TFile* savefile = new TFile(savename, "RECREATE");

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
 
    timerWatch.Stop();
    std::cout << "### DONE " << timerWatch.RealTime() << " seconds" << std::endl;

    return 0;
}
