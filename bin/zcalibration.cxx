// zcalibration.cxx, plot zmass and resolution in data vs phi+, phi-, eta+, eta-

#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "EventSelection.h"
#include "MuonSelection.h"
#include "CategorySelection.h"
#include "SampleDatabase.cxx"
#include "ThreadPool.hxx"

#include "TLorentzVector.h"
#include "TSystem.h"
#include "TROOT.h"

#include <sstream>
#include <map>
#include <vector>
#include <utility>

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

void initPlotSettings(int varNumber, float& fitsig, float& massmin, float& massmax, int& massbins, 
                      float& xmin, float& xmax, int& xbins, TString& xname)
{
    massmin = 86.2;
    massmax = 96.2;
    massbins = 50;
    xbins = 25;
    fitsig = 1;

    if(varNumber == 0)
    {
        xname = "phi_plus";
        xmin = -3.14;
        xmax = 3.14;
    }

    if(varNumber == 1)
    {
        xname = "phi_minus";
        xmin = -3.14;
        xmax = 3.14;
    }

    if(varNumber == 2)
    {
        xname = "eta_plus";
        xmin = -2.4;
        xmax = 2.4;
    }

    if(varNumber == 3)
    {
        xname = "eta_minus";
        xmin = -2.4;
        xmax = 2.4;
    }
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

int main(int argc, char* argv[])
{
    gROOT->SetBatch();
    int varNumber = 0;
    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) ss >> varNumber;
    }   

    float nthreads = 1;
    float luminosity = 36814;
    float reductionFactor = 10;
    std::map<TString, Sample*> samples;
    std::vector<Sample*> samplevec;
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

    ///////////////////////////////////////////////////////////////////
    // Plot Settings Depending on Input--------------------------------
    ///////////////////////////////////////////////////////////////////
    
    TString xname;
    float fitsig, massmin, massmax, xmin, xmax;
    int massbins, xbins;
    initPlotSettings(varNumber, fitsig, massmin, massmax, massbins, xmin, xmax, xbins, xname);

    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // gather samples map from SamplesDatabase.cxx
    GetSamples(samples, "UF");

    ///////////////////////////////////////////////////////////////////
    // PREPROCESSING: SetBranchAddresses-------------------------------
    ///////////////////////////////////////////////////////////////////

    // Loop through all of the samples to do some pre-processing
    // Add only data samples to the vector we loop over to categorize

    std::cout << std::endl;
    std::cout << "======== Preprocess the samples... " << std::endl;
    std::cout << std::endl;

    //makePUHistos(samples);

    for(auto &i : samples)
    {
        if(i.second->sampleType != "data") continue;

        // Output some info about the current file
        std::cout << "  /// Using sample " << i.second->name << std::endl;
        std::cout << std::endl;
        std::cout << "    sample name:       " << i.second->name << std::endl;
        std::cout << "    sample file:       " << i.second->filename << std::endl;
        std::cout << "    pileup file:       " << i.second->pileupfile << std::endl;
        std::cout << "    nOriginal:         " << i.second->nOriginal << std::endl;
        std::cout << "    N:                 " << i.second->N << std::endl;
        std::cout << "    nOriginalWeighted: " << i.second->nOriginalWeighted << std::endl;
        std::cout << std::endl;

        i.second->setBranchAddresses(1);
        samplevec.push_back(i.second);
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    //std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 
    
    std::cout << "@@@ nCPUs Available: " << getNumCPUs() << std::endl;
    std::cout << "@@@ nCPUs used     : " << nthreads << std::endl;
    std::cout << "@@@ nSamples used  : " << samplevec.size() << std::endl;

    std::cout << std::endl;
    std::cout << "======== Plot Configs ========" << std::endl;
    std::cout << "xvar         : " << xname << std::endl;
    std::cout << "xmin         : " << xmin << std::endl;
    std::cout << "xmax         : " << xmax << std::endl;
    std::cout << "xbins        : " << xbins << std::endl;
    std::cout << std::endl;
    std::cout << "massmin      : " << massmin << std::endl;
    std::cout << "massmax      : " << massmax << std::endl;
    std::cout << "massbins     : " << massbins << std::endl;
    std::cout << std::endl;
  
    auto makePlotsForSample = [varNumber, xname, fitsig, massbins, massmin, massmax, xbins, xmin, xmax, luminosity, reductionFactor](Sample* s)
    {
      // Output some info about the current file
      std::cout << Form("  /// Processing %s \n", s->name.Data());

      // cuts to apply to the events in the sample
      Run2MuonSelectionCuts     run2MuonSelection;
      Run2EventSelectionCuts80X run2EventSelectionData(true);

      // will probably want one for each type of mass
      ZCalibration* zcal_pf   = new ZCalibration(xname, s->name+"_mass_PF", fitsig, massmin, massmax, massbins, xmin, xmax, xbins);
      ZCalibration* zcal_roch = new ZCalibration(xname, s->name+"_mass_Roch", fitsig, massmin, massmax, massbins, xmin, xmax, xbins);
      ZCalibration* zcal_kamu = new ZCalibration(xname, s->name+"_mass_KaMu", fitsig, massmin, massmax, massbins, xmin, xmax, xbins);

      ///////////////////////////////////////////////////////////////////
      // HISTOGRAMS TO FILL ---------------------------------------------
      ///////////////////////////////////////////////////////////////////

      for(unsigned int i=0; i<s->N/reductionFactor; i++)
      {
        ///////////////////////////////////////////////////////////////////
        // GET INFORMATION ------------------------------------------------
        ///////////////////////////////////////////////////////////////////
        
        // only load essential information for the first set of cuts 
        s->branches.recoDimuCands->GetEntry(i);
        s->branches.recoMuons->GetEntry(i);

        // need a dimuon candidate to fill the zmass
        if(s->vars.recoDimuCands->size() < 1) continue;
        bool found_good_dimuon = false;

        // find the first good dimuon candidate and fill info
        for(auto& dimu: (*s->vars.recoDimuCands))
        {
          s->vars.dimuCand = &dimu;

          ///////////////////////////////////////////////////////////////////
          // CUTS  ----------------------------------------------------------
          ///////////////////////////////////////////////////////////////////

          if(!run2EventSelectionData.evaluate(s->vars) && s->sampleType.Contains("data"))
          { 
              continue; 
          }
          if(!s->vars.recoMuons->at(s->vars.dimuCand->iMu1).isTightID || !s->vars.recoMuons->at(s->vars.dimuCand->iMu2).isTightID)
          { 
              continue; 
          }
          if(!run2MuonSelection.evaluate(s->vars)) 
          {
              continue; 
          }

          found_good_dimuon = true;

          MuonInfo& mu1 = s->vars.recoMuons->at(s->vars.dimuCand->iMu1);
          MuonInfo& mu2 = s->vars.recoMuons->at(s->vars.dimuCand->iMu2);

          float phi_plus = (mu1.charge==1)?mu1.phi:mu2.phi;
          float phi_minus = (mu1.charge==-1)?mu1.phi:mu2.phi;
          
          float eta_plus = (mu1.charge==1)?mu1.eta:mu2.eta;
          float eta_minus = (mu1.charge==-1)?mu1.eta:mu2.eta;
           
          if(varNumber == 0) 
          {
              zcal_pf->fill(phi_plus, s->vars.dimuCand->mass_PF);
              zcal_roch->fill(phi_plus, s->vars.dimuCand->mass_Roch);
              zcal_kamu->fill(phi_plus, s->vars.dimuCand->mass_KaMu);
          }

          if(varNumber == 1) 
          {
              zcal_pf->fill(phi_minus, s->vars.dimuCand->mass_PF);
              zcal_roch->fill(phi_minus, s->vars.dimuCand->mass_Roch);
              zcal_kamu->fill(phi_minus, s->vars.dimuCand->mass_KaMu);
          }

          if(varNumber == 2) 
          {
              zcal_pf->fill(eta_plus, s->vars.dimuCand->mass_PF);
              zcal_roch->fill(eta_plus, s->vars.dimuCand->mass_Roch);
              zcal_kamu->fill(eta_plus, s->vars.dimuCand->mass_KaMu);
          }

          if(varNumber == 3) 
          {
              zcal_pf->fill(eta_minus, s->vars.dimuCand->mass_PF);
              zcal_roch->fill(eta_minus, s->vars.dimuCand->mass_Roch);
              zcal_kamu->fill(eta_minus, s->vars.dimuCand->mass_KaMu);
          }

          if(found_good_dimuon) break;

        } // end dimucand loop
      } // end events loop
      
      std::vector<ZCalibration*>* returnVector = new std::vector<ZCalibration*>();
      returnVector->push_back(zcal_pf);
      returnVector->push_back(zcal_roch);
      returnVector->push_back(zcal_kamu);
 
      std::cout << Form("  /// Checking histos within thread %s \n", s->name.Data());
      for(auto& zcal: (*returnVector))
      {
        for(unsigned int i=0; i<zcal->histos.size(); i++)
        {
            std::cout << zcal->histos[i]->GetName() << ": " << zcal->histos[i]->Integral() << std::endl;
        }
      }

      std::cout << Form("  /// Done processing %s \n", s->name.Data());
      return returnVector;

    };

   ///////////////////////////////////////////////////////////////////
   // SAMPLE PARALLELIZATION------ ----------------------------------
   ///////////////////////////////////////////////////////////////////

    TList* netlist = new TList();
    ThreadPool pool(nthreads);
    std::vector< std::future< std::vector<ZCalibration*>* > > results;

    TStopwatch timerWatch;
    timerWatch.Start();

    for(auto &s : samplevec)
        results.push_back(pool.enqueue(makePlotsForSample, s));

    // Create the stack and ratio plot    
    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile("rootfiles/zcalibration_"+xname+"_data_8_0_X.root", "RECREATE");
    savefile->cd();

    for(auto && result: results)
    {
        for(auto zcal: (*result.get()))
        {
            std::cout << "Checking histos before plotting... " << std::endl;
            for(unsigned int i=0; i<zcal->histos.size(); i++)
            {
                std::cout << zcal->histos[i]->GetName() << ": " << zcal->histos[i]->Integral() << std::endl;
            }
            std::cout << "Fitting histos, e.g. " << std::endl;
            zcal->plot();
            std::cout << "Done fitting histos, e.g. " << std::endl;

            // Add all of the plots to the tlist
            std::cout << "Adding histos, e.g. " << zcal->histos[0]->GetName() << std::endl;
            for(unsigned int i=0; i<zcal->histos.size(); i++)
            {
                std::cout << zcal->histos[i]->GetName() << ": " << zcal->histos[i]->Integral() << std::endl;
                zcal->histos[i]->Write();
            }
            
            std::cout << "Adding plots ... " << std::endl;
            std::cout << "mean" << std::endl;
            std::cout << "mean: " << zcal->mean_vs_x->GetName() << std::endl;
            std::cout << "resolution" << std::endl;
            std::cout << "res:  " << zcal->resolution_vs_x->GetName() << std::endl;
            zcal->mean_vs_x->Write();
            zcal->resolution_vs_x->Write();
            std::cout << "Done adding plots ... " << std::endl;
        }
    }

    savefile->Write();
    savefile->Close();

    timerWatch.Stop();
    std::cout << "### DONE " << timerWatch.RealTime() << " seconds" << std::endl;

    return 0;
}
