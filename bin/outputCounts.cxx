// output the signal, background and significance for the different categories 

// Missing HLT trigger info in CMSSW_8_0_X MC so we have to compare Data and MC in a different manner.
// We apply triggers to data but not to MC. Then scale MC for trigger efficiency.

#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "MuonSelection.h"
#include "EventSelection.h"
#include "CategorySelection.h"
#include "JetCollectionCleaner.h"
#include "MuonCollectionCleaner.h"
#include "EleCollectionCleaner.h"
#include "SampleDatabase.cxx"

#include "EventTools.h"
#include "PUTools.h"
#include "ThreadPool.hxx"

#include "TLorentzVector.h"
#include "TSystem.h"
#include <sstream>
#include <map>
#include <vector>
#include <utility>

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
    // save the errors for the histogram correctly so they depend upon 
    // the number used to fill originally rather than the scaling
    TH1::SetDefaultSumw2();
    int nthreads = 10;

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        ss >> nthreads;
    }   

    // Not sure that we need a map if we have a vector
    // Should use this as the main database and choose from it to make the vector
    std::map<TString, Sample*> samples;

    // Second container so that we can have a copy sorted by cross section.
    std::vector<Sample*> samplevec;

    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    float reductionFactor = 1;
    float luminosity = 36814;      // pb-1
    GetSamples(samples, "UF");

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
        if(i.second->sampleType == "data") continue;

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

        i.second->setBranchAddresses(2);
        samplevec.push_back(i.second);
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 

    std::cout << "@@@ nCPUs Available: " << getNumCPUs() << std::endl;
    std::cout << "@@@ nCPUs used     : " << nthreads << std::endl;
    std::cout << "@@@ nSamples used  : " << samplevec.size() << std::endl;

    auto outputSampleInfo = [luminosity, reductionFactor](Sample* s)
    {
      // Output some info about the current file
      std::cout << Form("  /// Processing %s \n", s->name.Data());

      int bin = -1;
      float interval = 0.5;
      int nbins = 20;
      float massmin = 120;
      float massmax = 130;

      // Objects to help with the cuts and selections
      JetCollectionCleaner      jetCollectionCleaner;
      MuonCollectionCleaner     muonCollectionCleaner;
      EleCollectionCleaner      eleCollectionCleaner;

      Run2MuonSelectionCuts  run2MuonSelection;
      Run2EventSelectionCuts run2EventSelectionMC;

      Categorizer* categorySelection = new CategorySelectionRun1();

      TString hkeyn;
      TString hkeyw;
      for(auto &c : categorySelection->categoryMap)
      {   
            //number of bins for the histogram
            int hbins = 21;

            // c.second is the category object, c.first is the category name
            TString hnamen = c.first+"_"+s->name+"_num";
            TString hnamew = c.first+"_"+s->name+"_weights";
            hkeyn = s->name+"_num";
            hkeyw = s->name+"_weights";

            // The number of events histo
            c.second.histoMap[hkeyn] = new TH1D(hnamen, hnamen, hbins, -1, 20);
            c.second.histoMap[hkeyn]->GetXaxis()->SetTitle("bin");
            c.second.histoList->Add(c.second.histoMap[hkeyn]);                                        
            if(s->sampleType.Contains("data")) c.second.dataList->Add(c.second.histoMap[hkeyn]);      
            if(s->sampleType.Contains("signal")) c.second.signalList->Add(c.second.histoMap[hkeyn]);  
            if(s->sampleType.Contains("background")) c.second.bkgList->Add(c.second.histoMap[hkeyn]); 

            // The sum of weights histo
            c.second.histoMap[hkeyw] = new TH1D(hnamew, hnamew, hbins, -1, 20);
            c.second.histoMap[hkeyw]->GetXaxis()->SetTitle("bin");
            c.second.histoList->Add(c.second.histoMap[hkeyw]);                                        
            if(s->sampleType.Contains("data")) c.second.dataList->Add(c.second.histoMap[hkeyw]);      
            if(s->sampleType.Contains("signal")) c.second.signalList->Add(c.second.histoMap[hkeyw]);  
            if(s->sampleType.Contains("background")) c.second.bkgList->Add(c.second.histoMap[hkeyw]); 
      }   

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
          MuonInfo& mu1 = s->vars.recoMuons->at(dimu.iMu1);
          MuonInfo& mu2 = s->vars.recoMuons->at(dimu.iMu2);

          // only want to train on events that are in the higgs mass window
          if(!(dimu.mass_PF > 110 && dimu.mass_PF < 160))
          {
              continue;
          }
          // usual cuts
          if(!mu1.isTightID || !mu2.isTightID)
          { 
              continue; 
          }
          if(!run2EventSelectionMC.evaluate(s->vars))
          { 
              continue; 
          }
          if(!run2MuonSelection.evaluate(s->vars)) 
          {
              continue; 
          }

          // Load the rest of the information needed
          s->branches.jets->GetEntry(i);
          s->branches.mht->GetEntry(i);
          s->branches.nVertices->GetEntry(i);
          s->branches.recoElectrons->GetEntry(i);

          s->branches.gen_wgt->GetEntry(i);
          s->branches.nPU->GetEntry(i);
          s->branches.pu_wgt->GetEntry(i);
          s->branches.eff_wgt->GetEntry(i);

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

          categorySelection->evaluate(s->vars);

          // dimuon event passes selections, set flag to true so that we only fill info for
          // the first good dimu candidate
          found_good_dimuon = true;

          // bin = -1 for events outside signal region, 
          // [0,nbins) for events inside signal region
          if(dimu.mass_PF < massmin || dimu.mass_PF >= massmax) bin = -1;
          else
          {
              float diff = dimu.mass_PF - massmin;
              bin =  diff/interval;
          }

          // Look at each category
          for(auto &c : categorySelection->categoryMap)
          {
              // if the event is in the current category then fill the category's histogram for the given sample and variable
              if(c.second.inCategory) c.second.histoMap[hkeyn]->Fill(bin);
              if(c.second.inCategory) c.second.histoMap[hkeyw]->Fill(bin, s->getWeight());
          } // end category loop

          if(found_good_dimuon) break; // only fill one dimuon, break from dimu cand loop

        } // end dimucand loop
      } // end event loop

      for(auto& c : categorySelection->categoryMap)
      {
          for(auto& h: c.second.histoMap)
          {
              if(h.first.Contains("weights"))
                  h.second->Scale(s->getScaleFactor(luminosity));
              //if(!c.second.hide) std::cout << Form("    %s, %s, %f\n", c.first.Data(), h.first.Data(), h.second->Integral());
          }
      }

      std::cout << Form("  /// Done processing %s \n", s->name.Data());
      return categorySelection;
    }; // end sample lambda function

    ThreadPool pool(nthreads);
    std::vector< std::future<Categorizer*> > results;

    TStopwatch timerWatch;
    timerWatch.Start();

    for(auto& s: samplevec)
        results.push_back(pool.enqueue(outputSampleInfo, s));

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
                cAll->categoryMap[category.first].histoMap[h.first] = h.second;
                cAll->categoryMap[category.first].histoList->Add(h.second);

                if(h.first.Contains("H2Mu")) cAll->categoryMap[category.first].signalList->Add(h.second);
                else                         cAll->categoryMap[category.first].bkgList->Add(h.second);
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

        TIter isig(c.second.signalList);
        TIter ibkg(c.second.bkgList);
        TObject* obj = 0;

        TList* signalListNum = new TList();
        TList* signalListW = new TList();
        TList* bkgListNum = new TList();
        TList* bkgListW = new TList();

        // filter num/weight histos to appropriate lists
        while(obj = isig())
        {
            TH1D* h = (TH1D*) obj;
            if(TString(h->GetName()).Contains("num")) signalListNum->Add(h);
            if(TString(h->GetName()).Contains("weights")) signalListW->Add(h);
        }
        while(obj = ibkg())
        {
            TH1D* h = (TH1D*) obj;
            if(TString(h->GetName()).Contains("num")) bkgListNum->Add(h);
            if(TString(h->GetName()).Contains("weights")) bkgListW->Add(h);
        }

        // num histos
        TH1D* hNetSignalNum = DiMuPlottingSystem::addHists(signalListNum, c.first+"_Num_Signal", "Num Signal");
        TH1D* hNetBkgNum    = DiMuPlottingSystem::addHists(bkgListNum,    c.first+"_Num_Bkg",    "Num Background");

        // weighted histos
        TH1D* hNetSignalW = DiMuPlottingSystem::addHists(signalListW, c.first+"_Net_Signal", "Net Signal");
        TH1D* hNetBkgW    = DiMuPlottingSystem::addHists(bkgListW,    c.first+"_Net_Bkg",    "Net Background");

        netlist->Add(hNetSignalW);
        netlist->Add(hNetBkgW);
        netlist->Add(hNetSignalNum);
        netlist->Add(hNetBkgNum);
    }

    TString savename = Form("rootfiles/significance_run1categories.root");

    std::cout << "  /// Saving plots to " << savename << " ..." << std::endl;
    std::cout << std::endl;

    TFile* savefile = new TFile(savename, "RECREATE");
    savefile->cd();
    netlist->Write();
    savefile->Close();

    timerWatch.Stop();
    std::cout << "### DONE " << timerWatch.RealTime() << " seconds" << std::endl;

    return 0;

}
