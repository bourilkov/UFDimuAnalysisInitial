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

#include "TSystem.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TROOT.h"

#include "SampleDatabase.cxx"
#include "ThreadPool.hxx"

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
  std::map<TString, Sample*> sampleMap;
  std::vector<Sample*> sampleVec;
  GetSamples(sampleMap, "UF", "ZJets_MG");

  int nthreads = 9;        // number of threads to use in parallelization
  float luminosity = 36814;
  float reductionFactor=1;
  std::map<TString, TH1D*> histMap;

  //TString xname = "dimu_mass";
  //int xbins = 140;
  //float xmin = 60;
  //float xmax = 200;

  TString xname = "HT";
  int xbins = 3000;
  float xmin = 0;
  float xmax = 3000;

  for(auto &s : sampleMap)
  {    
      //if(s.second->sampleType != "signal" && s.second->name != "ZJets_AMC" && s.second->name != "tt_ll_AMC" && s.second->sampleType != "data") continue;
      // Output some info about the current file
      std::cout << "  /// Using sample " << s.second->name << std::endl;
      std::cout << std::endl;
      std::cout << "    sample name:       " << s.second->name << std::endl;
      std::cout << "    sample file:       " << s.second->filenames[0] << std::endl;
      std::cout << "    pileup file:       " << s.second->pileupfile << std::endl;
      std::cout << "    nOriginal:         " << s.second->nOriginal << std::endl;
      std::cout << "    N:                 " << s.second->N << std::endl;
      std::cout << "    nOriginalWeighted: " << s.second->nOriginalWeighted << std::endl;
      std::cout << std::endl;

      s.second->setBranchAddresses(2); // tell the ttree to load the variable values into sample->vars
      sampleVec.push_back(s.second);   // add the sample to the vector of samples which we will run over

      // maps sample to histogram
      histMap[s.second->name] = new TH1D(s.second->name, s.second->name+Form(";%s;", xname), xbins, xmin, xmax);
  }    

  histMap["ZJets_MG_HT_0_70"] = new TH1D("ZJets_MG_HT_0_70", TString("ZJets_MG_HT_0_70")+Form(";%s;", xname), xbins, xmin, xmax);


  std::cout << "@@@ nCPUs Available: " << getNumCPUs() << std::endl;
  std::cout << "@@@ nCPUs used     : " << nthreads << std::endl;
  std::cout << "@@@ nSamples used  : " << sampleVec.size() << std::endl;


  auto makeHistoForSample = [&histMap, xname, reductionFactor, luminosity](Sample* s)
  {
    std::cout << Form("  /// Processing %s \n", s->name.Data());
    Run2MuonSelectionCuts  run2MuonSelection;
    Run2EventSelectionCuts run2EventSelection;

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
          // the dimuon candidate and the muons that make up the pair
          s->vars.dimuCand = &dimu; 
          MuonInfo& mu1 = s->vars.recoMuons->at(s->vars.dimuCand->iMu1);
          MuonInfo& mu2 = s->vars.recoMuons->at(s->vars.dimuCand->iMu2);

          ///////////////////////////////////////////////////////////////////
          // CUTS  ----------------------------------------------------------
          ///////////////////////////////////////////////////////////////////

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

          // dimuon event passes selections, set flag to true so that we only fill info for
          // the first good dimu candidate
          found_good_dimuon = true;

          // Load MC branches if MC
          if(s->sampleType!="data")
          {
              s->branches.gen_wgt->GetEntry(i);
              s->branches.nPU->GetEntry(i);
              s->branches.pu_wgt->GetEntry(i);
              s->branches.eff_wgt->GetEntry(i);
              s->branches.lhe_ht->GetEntry(i);
          }

          double varvalue = -999;
          if(xname == "dimu_mass") varvalue = dimu.mass;
          if(xname == "HT") varvalue = s->vars.lhe_ht;

          histMap[s->name]->Fill(varvalue, s->getWeight());
          if(s->name == "ZJets_MG" && s->vars.lhe_ht < 70)
          {
              histMap["ZJets_MG_HT_0_70"]->Fill(varvalue, s->getWeight()); 
          }

          if(found_good_dimuon) break; // only fill one dimuon, break from dimu cand loop

        } // end dimu cand loop //
    } // end event loop //

    histMap[s->name]->Scale(s->getScaleFactor(luminosity));
    if(s->name == "ZJets_MG")
        histMap["ZJets_MG_HT_0_70"]->Scale(s->getScaleFactor(luminosity));

    std::cout << Form("  /// Done processing %s \n", s->name.Data());
    return 0;
  };

  ///////////////////////////////////////////////////////////////////
  // PARALLELIZE BY SAMPLE -----------------------------------------
  ///////////////////////////////////////////////////////////////////

  ThreadPool pool(nthreads);
  std::vector< std::future<int> > results;

  TStopwatch timerWatch;
  timerWatch.Start();

  for(auto &s : sampleVec)
      results.push_back(pool.enqueue(makeHistoForSample, s));

  for(auto && result: results)
      std::cout << result.get();

  std::cout << std::endl;

  TList* ht = new TList();
  TList* stacklist = new TList();
  for(auto& h: histMap)
  {
      std::cout << Form("%s: %f\n", h.first.Data(), h.second->Integral());
  }

  ht->Add(histMap["ZJets_MG_HT_2500_inf"]);
  ht->Add(histMap["ZJets_MG_HT_1200_2500"]);
  ht->Add(histMap["ZJets_MG_HT_800_1200"]);
  ht->Add(histMap["ZJets_MG_HT_600_800"]);
  ht->Add(histMap["ZJets_MG_HT_400_600"]);
  ht->Add(histMap["ZJets_MG_HT_200_400"]);
  ht->Add(histMap["ZJets_MG_HT_100_200"]);
  ht->Add(histMap["ZJets_MG_HT_70_100"]);
  ht->Add(histMap["ZJets_MG_HT_0_70"]);

  TH1D* hNetHT   = DiMuPlottingSystem::addHists(ht, "Net_DY_HT", "Net_DY_HT");

  ht->Add(histMap["ZJets_MG"]);
  std::cout << std::endl;

  TCanvas* stack = DiMuPlottingSystem::stackedHistogramsAndRatio(ht, "DY_Compare", "DY_Compare", xname, "Num Entries", true, false, "Inclusive/Stitched");
  stacklist->Add(stack);

  TString savename = "rootfiles/compare_drell_yan.root";
  std::cout << "  /// Saving plots to " << savename << " ..." << std::endl;

  TFile* f = new TFile(savename, "RECREATE");
  f->cd();

  stacklist->Write();
  hNetHT->Write();
  ht->Write();

  f->Close();

  timerWatch.Stop();
  std::cout << "### DONE " << timerWatch.RealTime() << " seconds" << std::endl;
  return 0;
}
