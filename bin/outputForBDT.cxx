/////////////////////////////////////////////////////////////////////////////
//                           outputForBDT.cxx                              //
//=========================================================================//
//                                                                         //
// output the events in a small window around the higgs mass = 125 GeV     //
// and their variables so that we can train a bdt to make our              //
// categories for us.                                                      //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////


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
#include "TNtuple.h"
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

    //TString whichDY = "dyAMC";
    TString whichDY = "dyMG";
    GetSamples(samples, "UF", "ALL_"+whichDY);

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
        //if(i.second->name != "H2Mu_gg") continue;

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

    auto outputSampleInfo = [whichDY, luminosity, reductionFactor](Sample* s)
    {
      // Output some info about the current file
      std::cout << Form("  /// Processing %s \n", s->name.Data());

      int bin = -1;
      float interval = 0.5;
      int nbins = 20;
      float massmin = 120;
      float massmax = 130;

      //TFile* tfile = new TFile("rootfiles/bdt/"+s->name+"_bdt_training.root", "RECREATE");
      //tfile->cd();
      //TNtuple* ntuple = new TNtuple("theNtuple", "theNtuple");

      std::map<TString, double> vars;
      vars["bin"] = -999;
      vars["is_signal"] = -999;
      vars["weight"] = -999;

      // loop through all features in the varset (Sample->vars)
      // and initialize values
      for(auto& item: s->vars.varMap)
          vars[item.first.c_str()] = -999;

      // !!!! output first line of csv to file
      std::ofstream file(Form("csv/bdtcsv/%s_bdt_training_%s.csv", s->name.Data(), whichDY.Data()), std::ofstream::out);
      file << EventTools::outputMapKeysCSV(vars).Data() << std::endl;

      // ntuple requires list of variables to be separated by ":" rather than ","
      TString ntuplevars = EventTools::outputMapKeysCSV(vars).ReplaceAll(",", ":");
      TNtuple* ntuple = new TNtuple(s->name.Data(), s->name.Data(), ntuplevars.Data());

      // Objects to help with the cuts and selections
      JetCollectionCleaner      jetCollectionCleaner;
      MuonCollectionCleaner     muonCollectionCleaner;
      EleCollectionCleaner      eleCollectionCleaner;

      Run2MuonSelectionCuts  run2MuonSelection;
      Run2EventSelectionCuts run2EventSelection;

      Categorizer* categorySelection = new CategorySelectionRun1();

      for(unsigned int i=0; i<s->N/reductionFactor; i++)
      {

        // We are stitching together zjets_ht from 70-inf. We use the inclusive for
        // ht from 0-70.
        s->branches.lhe_ht->GetEntry(i);
        if(s->name == "ZJets_MG" && s->vars.lhe_ht >= 70) continue;

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

          // set pt calibration type: PF, Roch, or KaMu
          s->vars.setCalibrationType("PF");

          MuonInfo& mu1 = s->vars.recoMuons->at(dimu.iMu1);
          MuonInfo& mu2 = s->vars.recoMuons->at(dimu.iMu2);

          // only want to train on events that are in the higgs mass window
          if(!(dimu.mass > 110 && dimu.mass < 160))
          {
              continue;
          }
          // usual cuts
          if(!mu1.isMediumID || !mu2.isMediumID)
          { 
              continue; 
          }
          if(!run2EventSelection.evaluate(s->vars))
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
          s->branches.met->GetEntry(i);
          s->branches.nVertices->GetEntry(i);
          s->branches.recoElectrons->GetEntry(i);

          s->branches.nPU->GetEntry(i);
          s->branches.getEntryWeightsMC(i); // gen_wgt, pu_wgt and hlt, iso, mu_id scale factors

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

          // dimuon event passes selections, set flag to true so that we only fill info for
          // the first good dimu candidate
          found_good_dimuon = true;

          // bin = -1 for events outside signal region, 
          // [0,nbins) for events inside signal region
          if(dimu.mass < massmin || dimu.mass >= massmax) bin = -1;
          else
          {
              float diff = dimu.mass - massmin;
              bin =  diff/interval;
          }

          // nonfeature variables
          vars["bin"] = bin;
          vars["is_signal"] = s->sampleType.Contains("signal")?1:0;
          vars["weight"] = s->getLumiScaleFactor(luminosity)*s->getWeight();

          //std::cout << i << " !!! SETTING JETS " << std::endl;
          //s->vars.setJets();    // jets sorted and paired by mjj, turn this off to simply take the leading two jets
          s->vars.setVBFjets();   // jets sorted and paired by vbf criteria

          // put the actual value for each feature into our map
          // then output the map to CSV later
          for(auto& item: s->vars.varMap)
          {
              // item.first is the name of the feature
              //std::cout << item.first.c_str() << std::endl;
              vars[item.first.c_str()] = s->vars.getValue(item.first);
          }

          if(false)
            EventTools::outputEvent(s->vars, (*categorySelection));

          // !!!! output event info to file
          file << EventTools::outputMapValuesCSV(vars).Data() << std::endl;

          std::vector<Float_t> varvalues;
          for(auto& item: vars)
              varvalues.push_back((Float_t)item.second);

          // fill the ntuple
          ntuple->Fill(&varvalues[0]);

          if(found_good_dimuon) break; // only fill one dimuon, break from dimu cand loop

        } // end dimucand loop
      } // end event loop

      file.close();

      std::cout << Form("  /// Done processing %s \n", s->name.Data());
      return ntuple;
    }; // end sample lambda function

    ThreadPool pool(nthreads);
    std::vector< std::future<TNtuple*> > results;

    for(auto& s: samplevec)
        results.push_back(pool.enqueue(outputSampleInfo, s));

    std::vector<TNtuple*> ntuples;
    for(auto && result: results)  // get each ntuple and save it to a tfile
    {
        TNtuple* ntuple = result.get();
        ntuples.push_back(ntuple);
    }

    for(auto& ntuple: ntuples)
    {
        TString sname = ntuple->GetName();
        ntuple->SetName("theNtuple");
        ntuple->SetTitle("theNtuple");
        TString filename = Form("rootfiles/bdt/%s_bdt_training_%s.root", sname.Data(), whichDY.Data());
        std::cout << Form("  /// Writing ntuple to %s \n", filename.Data());
        TFile* tfile = new TFile(filename,"RECREATE");
        tfile->cd();
        ntuple->Write();
        tfile->Write();
        tfile->Close();
    }
}
