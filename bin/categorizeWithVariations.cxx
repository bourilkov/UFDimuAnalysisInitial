/////////////////////////////////////////////////////////////////////////////
//                           categorizeWithVariations.cxx                  //
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

struct Settings
{
// default settings here, may be overwritten by terminal input, see main() below

    int whichCategories = 1;              // run2categories = 1, run2categories = 2, "categories.xml" = 3 -> xmlcategories
    TString varname = "dimu_mass_Roch";   // the variable to plot, must match name in ../lib/VarSet
    int binning = 0;                      // binning = 1 -> plot dimu_mass from 110 to 160 for limit setting
                                          //  negative numbers unblind the data for limit setting
                                          //  see initPlotSettings above for more information
    
    bool rebin = true;        // rebin the ratio plots so that each point has small errors
    int nthreads = 10;        // number of threads to use in parallelization
    bool fitratio = 0;        // fit the ratio plot (data/mc) under the stack w/ a straight line
    
    TString xmlfile;          // filename for the xmlcategorizer, if you chose to use one 

    float luminosity = 36814;                // pb-1
    float reductionFactor = 50;              // reduce the number of events you run over in case you want to debug or some such thing

    TString whichDY = "dyAMC";               // use amc@nlo or madgraph for Drell Yan : {"dyAMC", "dyMG"}
    TString systematics;                     // which systematics to plot

    float min;
    float max;
    int bins;
};

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

TList* groupMC(TList* list, TString categoryName, TString suffix)
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
    TH1D* vh_histo = DiMuPlottingSystem::addHists(vh_list, categoryName+"VH"+suffix, "VH");
    TH1D* drell_yan_histo = DiMuPlottingSystem::addHists(drell_yan_list, categoryName+"Drell_Yan_"+suffix, "Drell Yan");
    TH1D* ttbar_histo = DiMuPlottingSystem::addHists(ttbar_list, categoryName+"TTbar_Plus_SingleTop"+suffix, "TTbar + SingleTop");
    TH1D* diboson_histo = DiMuPlottingSystem::addHists(diboson_list, categoryName+"Diboson_plus"+suffix, "Diboson +");

    grouped_list->AddFirst(vh_histo);
    grouped_list->Add(diboson_histo);
    grouped_list->Add(ttbar_histo);
    grouped_list->Add(drell_yan_histo);

    return grouped_list;
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

// set bins, min, max, and the var to plot based upon varNumber and settings.binning
void initPlotSettings(Settings& settings)
{
    settings.bins = 100;
    settings.min = -10;
    settings.max = 10;

    // dimu_mass
    if(settings.varname.Contains("dimu_mass"))
    {
        if(settings.binning == 0)       // blind data in 120-130 GeV, plot in 50-200 GeV to include z-peak
        {                      // 1 GeV bins, used mostly for validation plots
            settings.bins = 150;
            settings.min = 50;
            settings.max = 200;
        }
        else if(settings.binning == -1) // unblind data in 120-130 GeV for limit setting, 110-160 window
        {                      // 1 GeV bins
            settings.bins = 50;
            settings.min = 110;
            settings.max = 160;
        }
        else if(settings.binning == 1)  // blind data in 120-130 GeV, 110-160 window
        {                      // 1 GeV bins, study background fits
            settings.bins = 50;
            settings.min = 110;
            settings.max = 160;
        }
        else if(settings.binning == -2) // unblind data in 120-130 GeV for limit setting, 110-310 window
        {                      // 1 GeV bins
            settings.bins = 200;
            settings.min = 110;
            settings.max = 310;
        }
        else if(settings.binning == 2)  // blind data in 120-130 GeV, plot in 110 to 310 GeV range
        {                      // 2 GeV bins, study background fits
            settings.bins = 100;
            settings.min = 110;
            settings.max = 310;
        }
        else
        {
            settings.bins = 150;
            settings.min = 50;
            settings.max = 200;
        }
    }

    // mu_pt
    if(settings.varname.Contains("mu") && settings.varname.Contains("pt"))
    {
        settings.bins = 50;
        settings.min = 0;
        settings.max = 200;
        if(settings.varname.Contains("dimu")) settings.max = 300;
    }
 
    // mu_eta
    if(settings.varname.Contains("mu") && settings.varname.Contains("eta") && !settings.varname.Contains("dEta"))
    {
        settings.bins = 25;
        settings.min = -2.5;
        settings.max = 2.5;
    }

    // NPV
    if(settings.varname == "NPV")
    {
        settings.bins = 50;
        settings.min = 0;
        settings.max = 50;
    }

    // jet_pt
    if(settings.varname.Contains("jet") && settings.varname.Contains("pt"))
    {   
        settings.bins = 50;
        settings.min = 0;
        settings.max = 200;
    }   

    // jet_eta 
    if(settings.varname.Contains("jet") && settings.varname.Contains("eta"))
    {   
        settings.bins = 50;
        settings.min = -5; 
        settings.max = 5;
    }   

    // # of jets, ele, mu, etc
    if(settings.varname.Contains("nJets") || settings.varname.Contains("nVal") || settings.varname.Contains("nExtra") || settings.varname.Contains("nB") || settings.varname.Contains("nEle"))
    {   
        settings.bins = 11;
        settings.min = 0; 
        settings.max = 11;
    }   

    // m_jj
    if(settings.varname.Contains("m_jj") || settings.varname.Contains("m_bb") || (settings.varname.Contains("dijet") && settings.varname.Contains("mass")))
    {   
        settings.bins = 50;
        settings.min = 0; 
        settings.max = 2000;
    }   

    // dEta
    if(settings.varname.Contains("dEta"))
    {   
        settings.bins = 50;
        settings.min = -10; 
        settings.max = 10;
        if(settings.varname.Contains("mu"))
        {
            settings.bins = 25;
            settings.min = -5;
            settings.max = 5;
        }
    }   

    // dPhi
    if(settings.varname.Contains("dPhi"))
    {   
        settings.bins = 50;
        settings.min = -3.2; 
        settings.max = 3.2;

        if(settings.varname.Contains("Star"))
        {
            settings.min = 0;
            settings.max = 10;
        }
    }   

    // mT_b_MET/MT_had/mas_had
    if(settings.varname.Contains("mT") || settings.varname.Contains("MT") || settings.varname.Contains("mass_had"))
    {
        settings.bins = 50;
        settings.min = 0;
        settings.max = 2000;
    }

    // MET
    if(settings.varname == "MHT" || settings.varname == "MET")
    {
        settings.bins = 50;
        settings.min = 0;
        settings.max = 150;
    }

    // bdt_score
    if(settings.varname == "bdt_score")
    {
        settings.bins = 50;
        settings.min = -1;
        settings.max = 1;
    }

    if(settings.varname.Contains("abs"))
    {
        settings.bins = settings.bins/2;
        settings.min = 0;
    }
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

Categorizer* plotWithSystematic(TString systematic, Settings& settings)
{
    gROOT->SetBatch();
    // save the errors for the histogram correctly so they depend upon 
    // the number used to fill originally rather than the scaling
    TH1::SetDefaultSumw2();

    // Use this as the main database and choose from it to make the vector
    std::map<TString, Sample*> samples;

    // Second container so that we can have a copy sorted by cross section.
    std::vector<Sample*> samplevec;

    // Use this to plot some things if we wish
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    GetSamples(samples, "UF", "ALL_"+settings.whichDY);

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
        if(i.second->sampleType == "data" && systematic!="") continue;
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

        i.second->setBranchAddresses(systematic);          // tell the ttree to load the variable values into sample->vars
        samplevec.push_back(i.second);                              // add the sample to the vector of samples which we will run over
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 

    ///////////////////////////////////////////////////////////////////
    // Get Histo Information from Input -------------------------------
    ///////////////////////////////////////////////////////////////////

    // set nbins, min, max, and var to plot based upon input from the terminal: varNumber and settings.binning
    initPlotSettings(settings);

    std::cout << "@@@ nCPUs Available: " << getNumCPUs() << std::endl;
    std::cout << "@@@ nCPUs used     : " << settings.nthreads << std::endl;
    std::cout << "@@@ nSamples used  : " << samplevec.size() << std::endl;

    // print out settings.xmlfilename if using xmlcategorizer otherwise print out 1 or 2
    TString categoryString = Form("%d", settings.whichCategories);
    if(settings.whichCategories == 3) categoryString = settings.xmlfile;

    std::cout << std::endl;
    std::cout << "======== Plot Configs ========" << std::endl;
    std::cout << "categories     : " << categoryString << std::endl;
    std::cout << "systematic     : " << systematic << std::endl;
    std::cout << "var            : " << settings.varname << std::endl;
    std::cout << "min            : " << settings.min << std::endl;
    std::cout << "max            : " << settings.max << std::endl;
    std::cout << "bins           : " << settings.bins << std::endl;
    std::cout << "binning        : " << settings.binning << std::endl;
    std::cout << "reductionFactor: " << settings.reductionFactor << std::endl;
    std::cout << "whichDY        : " << settings.whichDY << std::endl;
    std::cout << "rebin ratio    : " << settings.rebin << std::endl;
    std::cout << "fit ratio      : " << settings.fitratio << std::endl;
    std::cout << std::endl;

    ///////////////////////////////////////////////////////////////////
    // Define Task for Parallelization -------------------------------
    ///////////////////////////////////////////////////////////////////

    auto makeHistoForSample = [settings, systematic](Sample* s)
    {

      // info to check that this event is different than the last event
      long long int last_run = -999;
      long long int last_event = -999;

      long long int this_run = -999;
      long long int this_event = -999;

      bool isblinded = true;              // negative settings.binning: unblind data in 120-130 GeV 
      if(settings.binning < 0) isblinded = false;

      // Output some info about the current file
      std::cout << Form("  /// Processing %s \n", s->name.Data());
      if(!s->vars.checkForVar(settings.varname.Data()))
      {
          std::cout << Form("  !!! %s is not a valid variable %s \n", settings.varname.Data(), s->name.Data());
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
      if(settings.whichCategories == 3)
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

      if(settings.whichCategories == 1) categorySelection = new CategorySelectionRun1();                  // run1 categories
      else if(settings.whichCategories == 2) categorySelection = new LotsOfCategoriesRun2();              // Adrian's proposed run2 categories
      else if(settings.whichCategories == 3 && settings.xmlfile.Contains("hybrid")) 
          categorySelection = new CategorySelectionBDT(settings.xmlfile);                                 // BDT based categories XML + object cuts
      else if(settings.whichCategories == 3) categorySelection = new XMLCategorizer(settings.xmlfile);    // BDT based categories XML only

      // set some flags
      bool isData = s->sampleType.EqualTo("data");
      bool isSignal = s->sampleType.EqualTo("signal");
      bool isMass = settings.varname.Contains("dimu_mass");

      // use pf, roch, or kamu values for selections, categories, and fill?
      TString pf_roch_or_kamu = "PF";
      if(settings.varname.Contains("PF")) pf_roch_or_kamu = "PF";
      else if(settings.varname.Contains("Roch")) pf_roch_or_kamu = "Roch";
      else if(settings.varname.Contains("KaMu")) pf_roch_or_kamu = "KaMu";

      ///////////////////////////////////////////////////////////////////
      // INIT HISTOGRAMS TO FILL ----------------------------------------
      ///////////////////////////////////////////////////////////////////

      // Keep track of which histogram to fill in the category
      TString hkey = s->name;

      // Different categories for the analysis
      // categorySelection has a categoryMap<TString, CategorizerObject>
      for(auto &c : categorySelection->categoryMap)
      {
          // c.second is the category object, c.first is the category name
          TString hname = c.second.name+"_"+s->name;
          if(systematic != "" && !isData) hname+="_"+systematic;

          // Set up the histogram for the category and variable to plot
          // Each category has a map and some lists to keep track of different histograms
          c.second.histoMap[hkey] = new TH1D(hname, hname, settings.bins, settings.min, settings.max);
          c.second.histoMap[hkey]->GetXaxis()->SetTitle(settings.varname);
          c.second.histoList->Add(c.second.histoMap[hkey]);                                      // need them ordered by xsec for the stack and ratio plot
          if(s->sampleType.EqualTo("data")) c.second.dataList->Add(c.second.histoMap[hkey]);     // data histo
          if(s->sampleType.EqualTo("signal")) c.second.signalList->Add(c.second.histoMap[hkey]); // signal histos
          if(s->sampleType.EqualTo("background")) c.second.bkgList->Add(c.second.histoMap[hkey]);// bkg histos

      }


      ///////////////////////////////////////////////////////////////////
      // LOOP OVER EVENTS -----------------------------------------------
      ///////////////////////////////////////////////////////////////////

      // Sift the events into the different categories and fill the histograms for each sample x category
      for(unsigned int i=0; i<s->N/settings.reductionFactor; i++)
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
          if(isSignal && settings.whichCategories==3 && settings.binning<0 && (s->vars.eventInfo->event % 2 == 1))
          {
              continue;
          }

          // only consider events w/ dimu_mass in the histogram range,
          // so the executable doesn't take forever, especially w/ tmva evaluation.
          if(isMass && (dimu.mass < settings.min || dimu.mass > settings.max))
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
          if(settings.whichCategories == 3)
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

              if(settings.varname.Contains("dimu_mass") || settings.varname == "bdt_score") 
              {
                  if(isData && dimu.mass > 120 && dimu.mass < 130 && isblinded) continue; // blind signal region

                  // if the event is in the current category then fill the category's histogram for the given sample and variable
                  c.second.histoMap[hkey]->Fill(s->vars.getValue(settings.varname.Data()), s->getWeight());
                  //std::cout << "    " << c.first << ": " << varvalue;
                  continue;
              }
              else
              {
                  c.second.histoMap[hkey]->Fill(s->vars.getValue(settings.varname.Data()), s->getWeight());
              }

          } // end category loop

          ////////////////////////////////////////////////////////////////////
          // DEBUG ----------------------------------------------------------
          
          // ouput pt, mass info etc for the event
          if(false)
          {
            EventTools::outputEvent(s->vars, *categorySelection);
            return categorySelection;
          }
           
          //------------------------------------------------------------------
          ////////////////////////////////////////////////////////////////////

          if(found_good_dimuon) break; // only fill one dimuon, break from dimu cand loop

        } // end dimu cand loop //
      } // end event loop //

      if(settings.whichCategories == 3) delete reader;
      //if(settings.whichCategories == 3) delete reader_multi;

      // Scale according to settings.luminosity and sample xsec now that the histograms are done being filled for that sample
      for(auto &c : categorySelection->categoryMap)
      {
          // Only used half of the signal events in this case, need to boost the normalization by 2 to make up for that
          if(settings.whichCategories == 3 && isSignal && settings.binning < 0) c.second.histoMap[hkey]->Scale(2*s->getLumiScaleFactor(settings.luminosity));
          else c.second.histoMap[hkey]->Scale(s->getLumiScaleFactor(settings.luminosity));
      }

      std::cout << Form("  /// Done processing %s \n", s->name.Data());
      delete s;
      return categorySelection;

    }; // done defining makeHistoForSample

   ///////////////////////////////////////////////////////////////////
   // PARALLELIZE BY SAMPLE -----------------------------------------
   ///////////////////////////////////////////////////////////////////

    ThreadPool pool(settings.nthreads);
    std::vector< std::future<Categorizer*> > results;

    for(auto &s : samplevec)
        results.push_back(pool.enqueue(makeHistoForSample, s));

   ///////////////////////////////////////////////////////////////////
   // Gather all the Histos into one Categorizer----------------------
   ///////////////////////////////////////////////////////////////////

    Categorizer* cAll = 0;
    if(settings.whichCategories == 1) cAll = new CategorySelectionRun1();                            // run1 categories 
    else if(settings.whichCategories == 2) cAll = new LotsOfCategoriesRun2();                        // Adrian's proposed run2 categories
    else if(settings.whichCategories == 3 && settings.xmlfile.Contains("hybrid")) 
        cAll = new CategorySelectionBDT(settings.xmlfile);                                           // BDT categories XML + Object cuts
    else if(settings.whichCategories == 3) cAll = new XMLCategorizer(settings.xmlfile);              // BDT categories using XML only

    // get histos from all categorizers and put them into one categorizer
    for(auto && categorizer: results)  // loop through each Categorizer object, one per sample
    {
        int i = 0;
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

    return cAll;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    Settings settings;


    ///////////////////////////////////////////////////////////////
    // Parse Arguments -------------------------------------------
    ///////////////////////////////////////////////////////////////

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        TString in = argv[i];
        TString option = in(0, in.First("="));
        option = option.ReplaceAll("--", "");
        TString value  = in(in.First("=")+1, in.Length());
        value = value.ReplaceAll("\"", "");
        ss << value.Data();

        if(option=="categories") 
        {
            if(value.Contains(".xml"))
            {
                settings.xmlfile = value;
                settings.whichCategories = 3;
            }
            else
                ss >> settings.whichCategories;
        }
        else if(option=="var")             settings.varname = value;
        else if(option=="binning")         ss >> settings.binning;
        else if(option=="nthreads")        ss >> settings.nthreads;
        else if(option=="rebinRatio")      ss >> settings.rebin;
        else if(option=="fitRatio")        ss >> settings.fitratio;
        else if(option=="reductionFactor") ss >> settings.reductionFactor;
        else if(option=="luminosity")      ss >> settings.luminosity;
        else if(option=="whichDY")         settings.whichDY = value;
        else
        {
            std::cout << Form("!!! %s is not a recognized option.", option) << std::endl;
        }
    }   

    std::vector<TString> systematics = {"", "JES_down"}; //, "JES_down", "PU_up", "PU_down"};
    TH1D* hNetData = 0; // keep the net data so that we can make DATA/MC stacks for other systematics
                        // without rerunning over the data samples

    TStopwatch timerWatch;
    timerWatch.Start();

    for(auto& systematic: systematics)
    {
        std::cout << std::endl;
        std::cout << "/////////////////////////////////////////////////////////////////////" << std::endl;
        std::cout << "/// SYSTEMATIC: " << systematic << std::endl;
        std::cout << "/////////////////////////////////////////////////////////////////////" << std::endl;
        std::cout << std::endl;

        Categorizer* cAll = plotWithSystematic(systematic, settings);

        ///////////////////////////////////////////////////////////////////
        // Save All of the Histos-----------------------------------------
        ///////////////////////////////////////////////////////////////////


        TList* varstacklist = new TList();   // list to save all of the stacks
        TList* signallist = new TList();     // list to save all of the signal histos
        TList* bglist = new TList();         // list to save all of the background histos
        TList* datalist = new TList();       // list to save all of the data histos
        TList* netlist = new TList();        // list to save all of the net histos

        TString suffix = "";
        if(systematic != "") suffix = "_"+systematic;

        std::cout << "About to loop through cAll" << std::endl;
        int i = 0;
        for(auto &c : cAll->categoryMap)
        {
            // some categories are intermediate and we don't want to save the plots for those
            if(c.second.hide) continue;

            // we need the data histo, the net signal, and the net bkg dimu mass histos for the datacards
            TH1D* hNetSignal = DiMuPlottingSystem::addHists(c.second.signalList, c.second.name+"_Net_Signal"+suffix, "Net Signal");
            TH1D* hNetBkg    = DiMuPlottingSystem::addHists(c.second.bkgList,    c.second.name+"_Net_Bkg"+suffix,    "Net Background");
            if(c.second.dataList->GetSize()!=0) 
                hNetData = DiMuPlottingSystem::addHists(c.second.dataList,   c.second.name+"_Net_Data"+suffix,   "Data");

            TList* groupedlist = groupMC(c.second.histoList, c.second.name, suffix);
            groupedlist->Add(hNetData);

            // Create the stack and ratio plot    
            TString stackname = c.second.name+"_stack"+suffix;

            //stackedHistogramsAndRatio(TList* list, TString name, TString title, TString xaxistitle, TString yaxistitle, bool settings.rebin = false, bool fit = true,
                                      //TString ratiotitle = "Data/MC", bool log = true, bool stats = false, int legend = 0);
            // stack signal, bkg, and data
            TCanvas* stack = DiMuPlottingSystem::stackedHistogramsAndRatio(groupedlist, stackname, stackname, settings.varname, "Num Entries", settings.rebin, settings.fitratio);
            varstacklist->Add(stack);

            // lists will contain signal, bg, and data histos for every category
            signallist->Add(c.second.signalList);
            bglist->Add(c.second.bkgList);
            datalist->Add(c.second.dataList);

            netlist->Add(hNetSignal);
            netlist->Add(hNetBkg);
            netlist->Add(groupedlist);
           
            stack->SaveAs("imgs/"+settings.varname+"_"+stackname+"_"+settings.whichDY+".png");
        }
        std::cout << std::endl;
        bool isblinded = settings.binning >= 0; // settings.binning == -1 means that we want dimu_mass from 110-160 unblinded
        TString blinded = "blinded";
        if(!isblinded) blinded = "UNBLINDED";

        TString xcategoryString = "";
        if(settings.whichCategories==3) 
        {
            Ssiz_t i = settings.xmlfile.Last('/');
            xcategoryString = settings.xmlfile(i+1, settings.xmlfile.Length()); // get the name of the settings.xmlfile without all the /path/to/dir/ business 
            xcategoryString = xcategoryString.ReplaceAll(".xml", "");
            xcategoryString = "_"+xcategoryString;
        }

        TString xvarname = settings.varname;
        if(settings.varname.Contains("dimu_mass")) xvarname=blinded+"_"+xvarname;
        TString savename = Form("rootfiles/validate_%s_%d_%d_categories%d%s_%d_%s.root", settings.varname.Data(), (int)settings.min, 
                                (int)settings.max, settings.whichCategories, xcategoryString.Data(), (int)settings.luminosity, settings.whichDY.Data());

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
    }
    return 0;
}
