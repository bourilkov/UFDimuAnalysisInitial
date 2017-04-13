/////////////////////////////////////////////////////////////////////////////
//                           categorize.cxx                                //
//=========================================================================//
//                                                                         //
// plot different variables in the run 1 or run 2 categories.              //
// may be used to look for discrepancies between data and mc or            //
// to make dimu_mass plots for limit setting or fitting.                   //
// outputs mc stacks with data overlayed and a ratio plot underneath.      //
// also saves the histos needed to make the mc stack, data.                //
// Also saves net BKG histo and net signal histo for limit setting,        //
// saves individuals as well.Missing HLT trigger info in CMSSW_8_0_X MC    //
// so we have to compare Data and MC in a different manner.                //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "EventSelection.h"
#include "MuonSelection.h"
#include "CategorySelection.h"
#include "JetCollectionCleaner.h"
#include "MuonCollectionCleaner.h"
#include "EleCollectionCleaner.h"

#include "EventTools.h"
#include "TMVATools.h"
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
#include "TROOT.h"
//#include "TStreamerInfo.h"

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

TList* groupMC(TList* list, TString categoryName)
{
// Group backgrounds into categories so there aren't a million of them in the legend
// drell_yan, ttbar + single top, diboson + rest, group VH together also

    categoryName+="_";

    TList* grouped_list = new TList();
    TList* drell_yan_list = new TList();
    TList* ttbar_list = new TList();
    TList* diboson_list = new TList();
    TList* vh_list = new TList();

    // strip data and get a sorted vector of mc samples
    for(unsigned int i=0; i<list->GetSize(); i++)
    {
        TH1D* hist = (TH1D*)list->At(i); 
        TString name = hist->GetName();
        // the sampleName is after the last underscore: categoryName_SampleName
        TString sampleName = name.ReplaceAll(categoryName, "");

        // filter out the data, since we add that to the end of this list later
        if(sampleName.Contains("Run") || sampleName.Contains("Data")) continue;
        else // group the monte carlo 
        {            
            if(sampleName.Contains("H2Mu"))
            {
                if(sampleName.Contains("gg")) 
                {
                    hist->SetTitle("GGF");
                    grouped_list->Add(hist); // go ahead and add the signal to the final list
                }
                if(sampleName.Contains("VBF")) 
                {
                    hist->SetTitle("VBF");
                    grouped_list->Add(hist); // go ahead and add the signal to the final list
                }
                if(sampleName.Contains("WH") || sampleName.Contains("ZH")) vh_list->Add(hist);
            }
            else if(sampleName.Contains("ZJets")) drell_yan_list->Add(hist);
            else if(sampleName.Contains("tt") || sampleName.Contains("tW") || sampleName.Contains("tZ")) ttbar_list->Add(hist);
            else diboson_list->Add(hist);
        }
    }
    // Don't forget to group VH together and retitle other signal samples
    TH1D* vh_histo = DiMuPlottingSystem::addHists(vh_list, categoryName+"VH", "VH");
    TH1D* drell_yan_histo = DiMuPlottingSystem::addHists(drell_yan_list, categoryName+"Drell_Yan", "Drell Yan");
    TH1D* ttbar_histo = DiMuPlottingSystem::addHists(ttbar_list, categoryName+"TTbar_Plus_SingleTop", "TTbar + SingleTop");
    TH1D* diboson_histo = DiMuPlottingSystem::addHists(diboson_list, categoryName+"Diboson_plus", "Diboson +");

    grouped_list->AddFirst(vh_histo);
    grouped_list->Add(diboson_histo);
    grouped_list->Add(ttbar_histo);
    grouped_list->Add(drell_yan_histo);

    return grouped_list;
}


//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

TList* getSortedMC(TList* list, std::vector<Sample*>& sampleVec, TString categoryName)
{
// Sorts "list" according to the xsec of the sample, using sampleVec
// where sampleVec was already sorted according to xsec via std::sort(sampleVec, byXsec[]{...}; )

    categoryName+="_";

    // map sample name to sorted vector location, since we already know the order
    // via sampleVec
    std::map<TString, unsigned int> xsecMap;
    for(unsigned int i=0; i< sampleVec.size(); i++)
    {
        xsecMap[sampleVec[i]->name] = i;
    }

    // strip data and get a sorted vector of mc samples
    std::vector<TH1D*> mcVec(sampleVec.size());
    unsigned int ndata = 0;
    for(unsigned int i=0; i<list->GetSize(); i++)
    {
        TH1D* hist = (TH1D*)list->At(i); 
        TString name = hist->GetName();
        // the sampleName is after the last underscore: categoryName_SampleName
        TString sampleName = name.ReplaceAll(categoryName, "");

        // Just count the number of data samples
        if(sampleName.Contains("Run") || sampleName.Contains("Data"))
            ndata++;
        else // put the mc histogram in its sorted location
        {            
            mcVec[xsecMap[sampleName]] = hist;    
        }
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

// set bins, min, max, and the var to plot based upon varNumber and binning
void initPlotSettings(TString varname, int binning, int& bins, float& min, float& max)
{
    bins = 100;
    min = -10;
    max = 10;

    // dimu_mass
    if(varname.Contains("dimu_mass"))
    {
        if(binning == 0)       // blind data in 120-130 GeV, plot in 50-200 GeV to include z-peak
        {                      // 1 GeV bins, used mostly for validation plots
            bins = 150;
            min = 50;
            max = 200;
        }
        else if(binning == -1) // unblind data in 120-130 GeV for limit setting, 110-160 window
        {                      // 1 GeV bins
            bins = 50;
            min = 110;
            max = 160;
        }
        else if(binning == 1)  // blind data in 120-130 GeV, 110-160 window
        {                      // 1 GeV bins, study background fits
            bins = 50;
            min = 110;
            max = 160;
        }
        else if(binning == -2) // unblind data in 120-130 GeV for limit setting, 110-310 window
        {                      // 1 GeV bins
            bins = 200;
            min = 110;
            max = 310;
        }
        else if(binning == 2)  // blind data in 120-130 GeV, plot in 110 to 310 GeV range
        {                      // 2 GeV bins, study background fits
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
    }

    // mu_pt
    if(varname.Contains("mu") && varname.Contains("pt"))
    {
        bins = 50;
        min = 0;
        max = 200;
        if(varname.Contains("dimu")) max = 300;
    }
 
    // mu_eta
    if(varname.Contains("mu") && varname.Contains("eta") && !varname.Contains("dEta"))
    {
        bins = 25;
        min = -2.5;
        max = 2.5;
    }

    // NPV
    if(varname == "NPV")
    {
        bins = 50;
        min = 0;
        max = 50;
    }

    // jet_pt
    if(varname.Contains("jet") && varname.Contains("pt"))
    {   
        bins = 50;
        min = 0;
        max = 200;
    }   

    // jet_eta 
    if(varname.Contains("jet") && varname.Contains("eta"))
    {   
        bins = 50;
        min = -5; 
        max = 5;
    }   

    // # of jets, ele, mu, etc
    if(varname.Contains("nJets") || varname.Contains("nVal") || varname.Contains("nExtra") || varname.Contains("nB") || varname.Contains("nEle"))
    {   
        bins = 11;
        min = 0; 
        max = 11;
    }   

    // m_jj
    if(varname.Contains("m_jj") || varname.Contains("m_bb") || (varname.Contains("dijet") && varname.Contains("mass")))
    {   
        bins = 50;
        min = 0; 
        max = 2000;
    }   

    // dEta
    if(varname.Contains("dEta"))
    {   
        bins = 50;
        min = -10; 
        max = 10;
        if(varname.Contains("mu"))
        {
            bins = 25;
            min = -5;
            max = 5;
        }
    }   

    // dPhi
    if(varname.Contains("dPhi"))
    {   
        bins = 50;
        min = -3.2; 
        max = 3.2;

        if(varname.Contains("Star"))
        {
            min = 0;
            max = 10;
        }
    }   

    // mT_b_MET/MT_had/mas_had
    if(varname.Contains("mT") || varname.Contains("MT") || varname.Contains("mass_had"))
    {
        bins = 50;
        min = 0;
        max = 2000;
    }

    // MET
    if(varname == "MHT" || varname == "MET")
    {
        bins = 50;
        min = 0;
        max = 150;
    }

    // bdt_score
    if(varname == "bdt_score")
    {
        bins = 50;
        min = -1;
        max = 1;
    }

    if(varname.Contains("abs"))
    {
        bins = bins/2;
        min = 0;
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

    int whichCategories = 1;         // run2categories = 1, run2categories = 2, "categories.xml" = 3 -> xmlcategories
    TString varname = "dimu_mass";   // the variable to plot, 0 is dimu_mass for instance
    int binning = 0;                 // binning = 1 -> plot dimu_mass from 110 to 160 for limit setting
                                     //  negative numbers unblind the data for limit setting
                                     //  see initPlotSettings above for more information

    bool rebin = true;        // rebin the ratio plots so that each point has small errors
    int nthreads = 10;        // number of threads to use in parallelization
    bool fitratio = 0;        // fit the ratio plot (data/mc) under the stack w/ a straight line

    TString xmlfile;          // filename for the xmlcategorizer, if you chose to use one 

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) 
        {
            size_t found = ss.str().find(".xml");
            if(found!=std::string::npos)
            {
                xmlfile = TString(ss.str().c_str());
                whichCategories = 3;
            }
            else
                ss >> whichCategories;
        }
        if(i==2) varname = TString(ss.str().c_str());
        if(i==3) ss >> binning;
        if(i==4) ss >> rebin;
        if(i==5) ss >> nthreads;
        if(i==6) ss >> fitratio;
    }   
    // Use this as the main database and choose from it to make the vector
    std::map<TString, Sample*> samples;

    // Second container so that we can have a copy sorted by cross section.
    std::vector<Sample*> samplevec;

    // Use this to plot some things if we wish
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

    float luminosity = 36814;       // pb-1
    float reductionFactor = 1;      // reduce the number of events you run over in case you want to debug or some such thing

    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // gather samples map from SamplesDatabase.cxx
    TString whichDY = "dyAMC";
    //TString whichDY = "dyMG";
    GetSamples(samples, "UF", "ALL_"+whichDY);

    ///////////////////////////////////////////////////////////////////
    // PREPROCESSING: SetBranchAddresses-------------------------------
    ///////////////////////////////////////////////////////////////////

    // Loop through all of the samples to do some pre-processing: filter out samples we don't want, set some branch addresses, etc
    // Add to the vector we loop over to categorize
    
    std::cout << std::endl;
    std::cout << "======== Preprocess the samples... " << std::endl;
    std::cout << std::endl;
    
    for(auto &i : samples)
    {
        //if(i.second->sampleType != "signal" && i.second->name != "ZJets_AMC" && i.second->name != "tt_ll_AMC" && i.second->sampleType != "data") continue;
        // Output some info about the current file
        std::cout << "  /// Using sample " << i.second->name << std::endl;
        std::cout << std::endl;
        std::cout << "    sample name:       " << i.second->name << std::endl;
        std::cout << "    sample file:       " << i.second->filenames[0] << std::endl;
        std::cout << "    pileup file:       " << i.second->pileupfile << std::endl;
        std::cout << "    nOriginal:         " << i.second->nOriginal << std::endl;
        std::cout << "    N:                 " << i.second->N << std::endl;
        std::cout << "    nOriginalWeighted: " << i.second->nOriginalWeighted << std::endl;
        std::cout << std::endl;

        i.second->setBranchAddresses();  // tell the ttree to load the variable values into sample->vars
        samplevec.push_back(i.second);   // add the sample to the vector of samples which we will run over
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 

    ///////////////////////////////////////////////////////////////////
    // Get Histo Information from Input -------------------------------
    ///////////////////////////////////////////////////////////////////

    // histo settings
    int bins;
    float min;
    float max;

    // set nbins, min, max, and var to plot based upon input from the terminal: varNumber and binning
    initPlotSettings(varname, binning, bins, min, max);

    std::cout << "@@@ nCPUs Available: " << getNumCPUs() << std::endl;
    std::cout << "@@@ nCPUs used     : " << nthreads << std::endl;
    std::cout << "@@@ nSamples used  : " << samplevec.size() << std::endl;

    // print out xmlfilename if using xmlcategorizer otherwise print out 1 or 2
    TString categoryString = Form("%d", whichCategories);
    if(whichCategories == 3) categoryString = xmlfile;

    std::cout << std::endl;
    std::cout << "======== Plot Configs ========" << std::endl;
    std::cout << "categories  : " << categoryString << std::endl;
    std::cout << "var         : " << varname << std::endl;
    std::cout << "min         : " << min << std::endl;
    std::cout << "max         : " << max << std::endl;
    std::cout << "bins        : " << bins << std::endl;
    std::cout << "binning     : " << binning << std::endl;
    std::cout << "rebin       : " << rebin << std::endl;
    std::cout << std::endl;

    ///////////////////////////////////////////////////////////////////
    // Define Task for Parallelization -------------------------------
    ///////////////////////////////////////////////////////////////////

    auto makeHistoForSample = [varname, binning, bins, min, max, rebin, whichCategories, 
                               whichDY, xmlfile, luminosity, reductionFactor](Sample* s)
    {

      // info to check that this event is different than the last event
      long long int last_run = -999;
      long long int last_event = -999;

      long long int this_run = -999;
      long long int this_event = -999;

      bool isblinded = true;              // negative binning: unblind data in 120-130 GeV 
      if(binning < 0) isblinded = false;

      // Output some info about the current file
      std::cout << Form("  /// Processing %s \n", s->name.Data());
      if(!s->vars.checkForVar(varname.Data()))
      {
          std::cout << Form("  !!! %s is not a valid variable %s \n", varname.Data(), s->name.Data());
      }


      /////////////////////////////////////////////////////
      // Load TMVA classifiers
     
      TString dir    = "classification/";
      //TString methodName = "BDTG_default";
      TString methodName = "BDTG_UF_v1";


      // sig vs bkg and multiclass (ggf, vbf, ... drell yan, ttbar) weight files
      TString weightfile = dir+"f_Opt_v1_all_sig_all_bkg_ge0j_BDTG_UF_v1.weights.xml";
      TString weightfile_multi = dir+"f_Opt_v1_multi_all_sig_all_bkg_ge0j_BDTG_UF_v1.weights.xml";

      /////////////////////////////////////////////////////
      // Book training and spectator vars into reader

      TMVA::Reader* reader = 0;
      std::map<TString, Float_t> tmap;
      std::map<TString, Float_t> smap;

      //TMVA::Reader* reader_multi = 0;
      //std::map<TString, Float_t> tmap_multi;
      //std::map<TString, Float_t> smap_multi;

      // load tmva binary classification and multiclass classifiers
      if(whichCategories == 3)
      {
          reader       = TMVATools::bookVars(methodName, weightfile, tmap, smap);
          //reader_multi = TMVATools::bookVars(methodName, weightfile_multi, tmap_multi, smap_multi);
      }

      ///////////////////////////////////////////////////////////////////
      // INIT Cuts and Categories ---------------------------------------
      ///////////////////////////////////////////////////////////////////
      
      // Objects to help with the cuts and selections
      JetCollectionCleaner      jetCollectionCleaner;
      MuonCollectionCleaner     muonCollectionCleaner;
      EleCollectionCleaner      eleCollectionCleaner;

      Run2MuonSelectionCuts  run2MuonSelection;
      Run2EventSelectionCuts run2EventSelection;

      Categorizer* categorySelection = 0;

      if(whichCategories == 1) categorySelection = new CategorySelectionRun1();         // run1 categories
      else if(whichCategories == 2) categorySelection = new LotsOfCategoriesRun2();     // Adrian's proposed run2 categories
      else if(whichCategories == 3 && xmlfile.Contains("hybrid")) categorySelection = new CategorySelectionBDT(xmlfile); // BDT based categories XML + object cuts
      else if(whichCategories == 3) categorySelection = new XMLCategorizer(xmlfile);    // BDT based categories XML only

      // set some flags
      bool isData = s->sampleType.EqualTo("data");
      bool isSignal = s->sampleType.EqualTo("signal");
      bool isMass = varname.Contains("dimu_mass");

      // use pf, roch, or kamu values for selections, categories, and fill?
      TString pf_roch_or_kamu = "PF";
      if(varname.Contains("PF")) pf_roch_or_kamu = "PF";
      else if(varname.Contains("Roch")) pf_roch_or_kamu = "Roch";
      else if(varname.Contains("KaMu")) pf_roch_or_kamu = "KaMu";

      ///////////////////////////////////////////////////////////////////
      // INIT HISTOGRAMS TO FILL ----------------------------------------
      ///////////////////////////////////////////////////////////////////

      // Keep track of which histogram to fill in the category
      TString hkey = s->name;

      // Different categories for the analysis
      // categorySelection has a categoryMap<TString, CategorizerObject>
      for(auto &c : categorySelection->categoryMap)
      {
          //number of bins for the histogram
          int hbins = bins;

          // c.second is the category object, c.first is the category name
          TString hname = c.first+"_"+s->name;

          // Set up the histogram for the category and variable to plot
          // Each category has a map and some lists to keep track of different histograms
          c.second.histoMap[hkey] = new TH1D(hname, hname, hbins, min, max);
          c.second.histoMap[hkey]->GetXaxis()->SetTitle(varname);
          c.second.histoList->Add(c.second.histoMap[hkey]);                                      // need them ordered by xsec for the stack and ratio plot
          if(s->sampleType.EqualTo("data")) c.second.dataList->Add(c.second.histoMap[hkey]);     // data histo
          if(s->sampleType.EqualTo("signal")) c.second.signalList->Add(c.second.histoMap[hkey]); // signal histos
          if(s->sampleType.EqualTo("background")) c.second.bkgList->Add(c.second.histoMap[hkey]);// bkg histos

      }


      ///////////////////////////////////////////////////////////////////
      // LOOP OVER EVENTS -----------------------------------------------
      ///////////////////////////////////////////////////////////////////

      // Sift the events into the different categories and fill the histograms for each sample x category
      for(unsigned int i=0; i<s->N/reductionFactor; i++)
      {
        // We are stitching together zjets_ht from 70-inf. We use the inclusive for
        // ht from 0-70.
        if(!isData)
        {
            s->branches.lhe_ht->GetEntry(i);
            if(s->name == "ZJets_MG" && s->vars.lhe_ht >= 70) continue;
        }

        // only load essential information for the first set of cuts 
        s->branches.muPairs->GetEntry(i);
        s->branches.muons->GetEntry(i);
        s->branches.eventInfo->GetEntry(i);

        // loop and find a good dimuon candidate
        if(s->vars.muPairs->size() < 1) continue;
        bool found_good_dimuon = false;

        // find the first good dimuon candidate and fill info
        for(auto& dimu: (*s->vars.muPairs))
        {
          // Reset the categorizer in preparation for the next event
          categorySelection->reset();

          // the dimuon candidate and the muons that make up the pair
          s->vars.dimuCand = &dimu; 
          MuonInfo& mu1 = s->vars.muons->at(s->vars.dimuCand->iMu1);
          MuonInfo& mu2 = s->vars.muons->at(s->vars.dimuCand->iMu2);

          // Selection cuts and categories use standard values e.g. mu.pt
          // to use PF, Roch, or KaMu for cuts and categories set these values 
          // to PF, Roch, or KaMu 
          s->vars.setCalibrationType(pf_roch_or_kamu);

          ///////////////////////////////////////////////////////////////////
          // CUTS  ----------------------------------------------------------
          ///////////////////////////////////////////////////////////////////

          // only use even signal events for limit setting 
          // as to separate training and evaluation events
          if(isSignal && whichCategories==3 && binning<0 && (s->vars.eventInfo->event % 2 == 1))
          {
              continue;
          }

          // only consider events w/ dimu_mass in the histogram range,
          // so the executable doesn't take forever, especially w/ tmva evaluation.
          if(isMass && (dimu.mass < min || dimu.mass > max))
          {
              continue;
          }

          // normal selections
          if(!run2EventSelection.evaluate(s->vars))
          { 
              continue; 
          }
          if(!mu1.isMediumID || !mu2.isMediumID)
          { 
              continue; 
          }
          if(!run2MuonSelection.evaluate(s->vars)) 
          {
              continue; 
          }

          // avoid double counting in RunF
          if(s->name == "RunF_1" && s->vars.eventInfo->run > 278801)
          {
              continue;
          }
          if(s->name == "RunF_2" && s->vars.eventInfo->run < 278802)
          {
              continue;
          }

          // Check that this event is different than the last event --------
          this_run = s->vars.eventInfo->run;
          this_event = s->vars.eventInfo->event;


          if(this_run == last_run && this_event == last_event)
              std::cout << Form("  !!! last run & event ==  this run & event: %d, %d, %d, %d, %s \n", last_run, last_event, this_run, this_event, s->name.Data());

          last_run = this_run;
          last_event = this_event;
          // ----------------------------------------------------------------

          // dimuon event passes selections, set flag to true so that we only fill info for
          // the first good dimu candidate
          found_good_dimuon = true; 

          ///////////////////////////////////////////////////////////////////
          // LOAD ALL BRANCHES, CLEAN COLLECTIONS ---------------------------
          ///////////////////////////////////////////////////////////////////

          // Load the rest of the information needed for run2 categories
          s->branches.getEntry(i);
          s->vars.setCalibrationType(pf_roch_or_kamu); // reloaded the branches, need to set mass,pt to correct calibrations again

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
          CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validMuons, 0.4);
          CollectionCleaner::cleanByDR(s->vars.validElectrons, s->vars.validMuons, 0.4);
          CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validElectrons, 0.4);

          //std::pair<int,int> e(s->vars.eventInfo.run, s->vars.eventInfo.event); // create a pair that identifies the event uniquely

          ///////////////////////////////////////////////////////////////////
          // CATEGORIZE -----------------------------------------------------
          ///////////////////////////////////////////////////////////////////

          // XML categories require the classifier score from TMVA
          if(whichCategories == 3)
          {
              //std::cout << i << " !!! SETTING JETS " << std::endl;
              //s->vars.setJets();    // jets sorted and paired by mjj, turn this off to simply take the leading two jets
              s->vars.bdt_out = TMVATools::getClassifierScore(reader, methodName, tmap, s->vars); // set tmva's bdt score

              // load multi results into varset
              //std::vector<float> bdt_multi_scores = TMVATools::getMulticlassScores(reader_multi, methodName, tmap_multi, s->vars);
              //s->vars.bdt_ggh_out = bdt_multi_scores[0];
              //s->vars.bdt_vbf_out = bdt_multi_scores[1];
              //s->vars.bdt_vh_out  = bdt_multi_scores[2];
              //s->vars.bdt_ewk_out = bdt_multi_scores[3];
              //s->vars.bdt_top_out = bdt_multi_scores[4];

              s->vars.setVBFjets();   // jets sorted and paired by vbf criteria
          }

          // Figure out which category the event belongs to
          categorySelection->evaluate(s->vars);

          // Look at each category, if the event belongs to that category fill the histogram for the sample x category
          for(auto &c : categorySelection->categoryMap)
          {
              // skip categories that we decided to hide (usually some intermediate categories)
              if(c.second.hide) continue;
              if(!c.second.inCategory) continue;

              ///////////////////////////////////////////////////////////////////
              // FILL HISTOGRAMS FOR REQUESTED VARIABLE -------------------------
              ///////////////////////////////////////////////////////////////////

              // now that we have a lot of these in the VarSet map we should just fill the s->vars.getValue("varName")
              // would save a lot of explicit code here ...

              if(varname.Contains("dimu_mass") || varname == "bdt_score") 
              {
                  if(isData && dimu.mass > 120 && dimu.mass < 130 && isblinded) continue; // blind signal region

                  // if the event is in the current category then fill the category's histogram for the given sample and variable
                  c.second.histoMap[hkey]->Fill(s->vars.getValue(varname.Data()), s->getWeight());
                  //std::cout << "    " << c.first << ": " << varvalue;
                  continue;
              }
              else
              {
                  c.second.histoMap[hkey]->Fill(s->vars.getValue(varname.Data()), s->getWeight());
              }

          } // end category loop

          ////////////////////////////////////////////////////////////////////
          // DEBUG ----------------------------------------------------------
          
          // ouput pt, mass info etc for the event
          if(false)
            EventTools::outputEvent(s->vars, *categorySelection);
           
          //------------------------------------------------------------------
          ////////////////////////////////////////////////////////////////////

          if(found_good_dimuon) break; // only fill one dimuon, break from dimu cand loop

        } // end dimu cand loop //
      } // end event loop //

      if(whichCategories == 3) delete reader;
      //if(whichCategories == 3) delete reader_multi;

      // Scale according to luminosity and sample xsec now that the histograms are done being filled for that sample
      for(auto &c : categorySelection->categoryMap)
      {
          // Only used half of the signal events in this case, need to boost the normalization by 2 to make up for that
          if(whichCategories == 3 && isSignal && binning < 0) c.second.histoMap[hkey]->Scale(2*s->getLumiScaleFactor(luminosity));
          else c.second.histoMap[hkey]->Scale(s->getLumiScaleFactor(luminosity));
      }

      std::cout << Form("  /// Done processing %s \n", s->name.Data());
      return categorySelection;

    }; // done defining makeHistoForSample

   ///////////////////////////////////////////////////////////////////
   // PARALLELIZE BY SAMPLE -----------------------------------------
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
    if(whichCategories == 1) cAll = new CategorySelectionRun1();      // run1 categories 
    else if(whichCategories == 2) cAll = new LotsOfCategoriesRun2();  // Adrian's proposed run2 categories
    else if(whichCategories == 3 && xmlfile.Contains("hybrid")) cAll = new CategorySelectionBDT(xmlfile); // BDT categories XML + Object cuts
    else if(whichCategories == 3) cAll = new XMLCategorizer(xmlfile);                                     // BDT categories using XML only

    // get histos from all categorizers and put them into one categorizer
    for(auto && categorizer: results)  // loop through each Categorizer object, one per sample
    {
        for(auto& category: categorizer.get()->categoryMap) // loop through each category for the given sample
        {
            // category.first is the category name, category.second is the Category object
            if(category.second.hide) continue;
            for(auto& h: category.second.histoMap) // loop through each histogram in this category's histo map
            {
                // std::cout << Form("%s: %f", h.first.Data(), h.second->Integral()) << std::endl;
                // we defined hkey = h.first as the sample name earlier so we have 
                // our histomap : category.histoMap<samplename, TH1D*>
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

    std::cout << "About to loop through cAll" << std::endl;
    for(auto &c : cAll->categoryMap)
    {
        // some categories are intermediate and we don't want to save the plots for those
        if(c.second.hide) continue;

        // we need the data histo, the net signal, and the net bkg dimu mass histos for the datacards
        TH1D* hNetSignal = dps->addHists(c.second.signalList, c.first+"_Net_Signal", "Net Signal");
        TH1D* hNetBkg    = dps->addHists(c.second.bkgList,    c.first+"_Net_Bkg",    "Net Background");
        TH1D* hNetData   = dps->addHists(c.second.dataList,   c.first+"_Net_Data",   "Data");

        TList* groupedlist = groupMC(c.second.histoList, c.first);
        groupedlist->Add(hNetData);

        //TIter next(groupedlist);
        //TObject* object = 0;

        //while ((object = next()))
        //{
        //    TH1D* h = (TH1D*) object;
        //    std::cout << Form("aboutToStack:: %s in sortedlist \n", h->GetName());
        //}   

        // Create the stack and ratio plot    
        TString cname = c.first+"_stack";
        //stackedHistogramsAndRatio(TList* list, TString name, TString title, TString xaxistitle, TString yaxistitle, bool rebin = false, bool fit = true,
                                  //TString ratiotitle = "Data/MC", bool log = true, bool stats = false, int legend = 0);
        // stack signal, bkg, and data
        TCanvas* stack = dps->stackedHistogramsAndRatio(groupedlist, cname, cname, varname, "Num Entries", rebin, fitratio);
        varstacklist->Add(stack);

        // lists will contain signal, bg, and data histos for every category
        signallist->Add(c.second.signalList);
        bglist->Add(c.second.bkgList);
        datalist->Add(c.second.dataList);

        netlist->Add(hNetSignal);
        netlist->Add(hNetBkg);
        netlist->Add(groupedlist);
       
        stack->SaveAs("imgs/"+varname+"_"+cname+"_"+whichDY+".png");
    }
    std::cout << std::endl;
    bool isblinded = binning >= 0; // binning == -1 means that we want dimu_mass from 110-160 unblinded
    TString blinded = "blinded";
    if(!isblinded) blinded = "UNBLINDED";

    TString xcategoryString = "";
    if(whichCategories==3) 
    {
        Ssiz_t i = xmlfile.Last('/');
        xcategoryString = xmlfile(i+1, xmlfile.Length()); // get the name of the xmlfile without all the /path/to/dir/ business 
        xcategoryString = xcategoryString.ReplaceAll(".xml", "");
        xcategoryString = "_"+xcategoryString;
    }

    if(varname.Contains("dimu_mass")) varname=blinded+"_"+varname;
    TString savename = Form("rootfiles/validate_%s_%d_%d_categories%d%s_%d_%s.root", varname.Data(), (int)min, 
                            (int)max, whichCategories, xcategoryString.Data(), (int)luminosity, whichDY.Data());

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
