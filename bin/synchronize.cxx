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
#include "TROOT.h"
//#include "TStreamerInfo.h"

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
    gROOT->SetBatch();
    // save the errors for the histogram correctly so they depend upon 
    // the number used to fill originally rather than the scaling
    TH1::SetDefaultSumw2();

    int nthreads = 4;        // number of threads to use in parallelization

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) ss >> nthreads;
    }   
    // Not sure that we need a map if we have a vector
    // Should use this as the main database and choose from it to make the vector
    std::map<TString, Sample*> samples;

    // Second container so that we can have a copy sorted by cross section.
    std::vector<Sample*> samplevec;

    // Use this to plot some things if we wish
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

    float luminosity = 36814;      // pb-1
    float reductionFactor = 1;     // reduce the number of events you run over in case you want to debug or some such thing

    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // gather samples map from SamplesDatabase.cxx
    GetSamples(samples, "UF");

    ///////////////////////////////////////////////////////////////////
    // PREPROCESSING: SetBranchAddresses-------------------------------
    ///////////////////////////////////////////////////////////////////

    // Loop through all of the samples to do some pre-processing
    // Add to the vector we loop over to categorize
    
    std::cout << std::endl;
    std::cout << "======== Preprocess the samples... " << std::endl;
    std::cout << std::endl;

    //makePUHistos(samples);
    
    for(auto &i : samples)
    {

        if(i.second->sampleType != "data" && i.second->name != "H2Mu_gg" && i.second->name != "ZZ_2l_2v") continue;
        if(i.second->sampleType == "data" && i.second->name != "RunD" && i.second->name != "RunH") continue;

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

        i.second->setBranchAddresses(2);
        samplevec.push_back(i.second);
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 

    ///////////////////////////////////////////////////////////////////
    // Get Histo Information from Input -------------------------------
    ///////////////////////////////////////////////////////////////////

    std::cout << "@@@ nCPUs Available: " << getNumCPUs() << std::endl;
    std::cout << "@@@ nCPUs used     : " << nthreads << std::endl;
    std::cout << "@@@ nSamples used  : " << samplevec.size() << std::endl;

    ///////////////////////////////////////////////////////////////////
    // Define Task for Parallelization -------------------------------
    ///////////////////////////////////////////////////////////////////

    auto synchSample = [reductionFactor](Sample* s)
    {
      // Output some info about the current file

      bool useMedium2016 = s->sampleType == "data" && s->name != "RunG" && s->name != "RunH";
      std::cout << Form("  /// Processing %s with useMedium2016 = %d \n", s->name.Data(), (int)useMedium2016);

      ///////////////////////////////////////////////////////////////////
      // INIT Cuts and Categories ---------------------------------------
      ///////////////////////////////////////////////////////////////////
      
      // Objects to help with the cuts and selections
      JetCollectionCleaner      jetCollectionCleaner;
      MuonCollectionCleaner     muonCollectionCleaner;
      EleCollectionCleaner      eleCollectionCleaner;
      muonCollectionCleaner.cUseMedium2016 = useMedium2016; // set the correct medium ID depending on the Run/MC

      Run2MuonSelectionCuts    synchMuonSelection;
      SynchEventSelectionCuts  synchEventSelection;

      Categorizer* categorySelection = new CategorySelectionSynch();

      ///////////////////////////////////////////////////////////////////
      // INIT HISTOGRAMS TO FILL ----------------------------------------
      ///////////////////////////////////////////////////////////////////

      // Keep track of which histogram to fill in the category
      TString hkey = s->name;

      // Different categories for the analysis
      for(auto &c : categorySelection->categoryMap)
      {
          // histogram settings
          float min = 0;
          float max = 2;
          int hbins = 2;

          // c.second is the category object, c.first is the category name
          TString hname = c.first+"_"+s->name;

          // Set up the histogram for the category and variable to plot
          c.second.histoMap[hkey] = new TH1D(hname, hname, hbins, min, max);
          c.second.histoMap[hkey]->GetXaxis()->SetTitle("count");
          c.second.histoList->Add(c.second.histoMap[hkey]);
      }


      ///////////////////////////////////////////////////////////////////
      // LOOP OVER EVENTS -----------------------------------------------
      ///////////////////////////////////////////////////////////////////

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
          // Reset the flags in preparation for the next event
          categorySelection->reset();

          // the dimuon candidate and the muons that make up the pair
          s->vars.dimuCand = &dimu; 
          MuonInfo& mu1 = s->vars.recoMuons->at(s->vars.dimuCand->iMu1);
          MuonInfo& mu2 = s->vars.recoMuons->at(s->vars.dimuCand->iMu2);

          dimu.mass = dimu.mass_PF;
          mu1.pt = mu1.pt_PF;
          mu2.pt = mu2.pt_PF;

          // Load the rest of the information needed for run2 categories
          s->branches.jets->GetEntry(i);
          s->branches.mht->GetEntry(i);
          s->branches.nVertices->GetEntry(i);
          s->branches.recoElectrons->GetEntry(i);
          //eventInfoBranch->GetEntry(i);

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
          CollectionCleaner::cleanByDR(s->vars.validBJets, s->vars.validMuons, 0.4);
          //CollectionCleaner::cleanByDR(s->vars.validElectrons, s->vars.validMuons, 0.4);
          //CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validElectrons, 0.4);
          
          ///////////////////////////////////////////////////////////////////
          // CUTS  ----------------------------------------------------------
          ///////////////////////////////////////////////////////////////////

          if(!synchMuonSelection.evaluate(s->vars)) 
          {
              continue; 
          }
          if(!synchEventSelection.evaluate(s->vars))
          { 
              continue; 
          }
          if( useMedium2016 && (!mu1.isMediumID2016 || !mu2.isMediumID2016) )
          { 
              continue; 
          }
          if( !useMedium2016 && (!mu1.isMediumID || !mu2.isMediumID) )
          { 
              continue; 
          }
          if(!mu1.isGlobal || !mu1.isTracker || !mu2.isGlobal || !mu2.isTracker)
          {
              continue;
          }


          // dimuon event passes selections, set flag to true so that we only fill info for
          // the first good dimu candidate
          found_good_dimuon = true; 

          // Figure out which category the event belongs to
          categorySelection->evaluate(s->vars);

          // Look at each category
          for(auto &c : categorySelection->categoryMap)
          {
              // skip categories
              if(c.second.hide) continue;
              if(!c.second.inCategory) continue;

              c.second.histoMap[hkey]->Fill(1, 1);

          } // end category loop

          if(false)
            // ouput pt, mass info etc for the event
            EventTools::outputEvent(s->vars, *categorySelection);

          if(found_good_dimuon) break; // only fill one dimuon, break from dimu cand loop

        } // end dimu cand loop
      } // end event loop

      std::cout << Form("  /// Done processing %s \n", s->name.Data());
      return categorySelection;

    }; // done defining synchSample

   ///////////////////////////////////////////////////////////////////
   // SAMPLE PARALLELIZATION------ ----------------------------------
   ///////////////////////////////////////////////////////////////////

    ThreadPool pool(nthreads);
    std::vector< std::future<Categorizer*> > results;

    TStopwatch timerWatch;
    timerWatch.Start();

    for(auto &s : samplevec)
        results.push_back(pool.enqueue(synchSample, s));

   ///////////////////////////////////////////////////////////////////
   // Gather all the Histos into one Categorizer----------------------
   ///////////////////////////////////////////////////////////////////

    Categorizer* cAll = new CategorySelectionRun1();

    // get histos from all categorizers and put them into one
    for(auto && categorizer: results)  // loop through each Categorizer object, one per sample
    {
        for(auto& category: categorizer.get()->categoryMap) // loop through each category
        {
            // category.first is the category name, category.second is the Category object
            if(category.second.hide) continue;
            for(auto& h: category.second.histoMap) // loop through each histogram in the category
            {
                // std::cout << Form("%s: %f", h.first.Data(), h.second->Integral()) << std::endl;
                // h.first is the sample name, category.histoMap<samplename, TH1D*>
                Sample* s = samples[h.first];

                cAll->categoryMap[category.first].histoMap[h.first] = h.second;
                cAll->categoryMap[category.first].histoList->Add(h.second);
                std::cout << Form("%s: %d\n", h.second->GetName(), (int)h.second->Integral());
            }
        }
        std::cout << Form("\n");
    }

   ///////////////////////////////////////////////////////////////////
   // Save All of the Histos-----------------------------------------
   ///////////////////////////////////////////////////////////////////
   TList* list = new TList();

    std::cout << "About to loop through cAll" << std::endl;
    for(auto &c : cAll->categoryMap)
    {
        // some categories are intermediate and we don't want to save the plots for those
        if(c.second.hide) continue;
        for(auto& item: c.second.histoMap)
        {
            list->Add(item.second);
        }
    }
    std::cout << std::endl;
    TString savename = Form("rootfiles/synchronize.root");

    std::cout << "  /// Saving plots to " << savename << " ..." << std::endl;
    std::cout << std::endl;

    TFile* savefile = new TFile(savename, "RECREATE");
    savefile->cd();
    list->Write();
    savefile->Close();
 
    timerWatch.Stop();
    std::cout << "### DONE " << timerWatch.RealTime() << " seconds" << std::endl;

    return 0;
}
