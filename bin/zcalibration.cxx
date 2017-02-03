// zcalibration.cxx, plot zmass and resolution in data vs phi+, phi-, eta+, eta-

#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "EventSelection.h"
#include "MuonSelection.h"
#include "CategorySelection.h"
#include "SampleDatabase.cxx"

#include "TLorentzVector.h"
#include "TSystem.h"

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
    int varNumber = 0;
    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) ss >> varNumber;
    }   

    float luminosity = 36814;
    float reductionFactor = 1;
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
    getSamples(luminosity, samples);

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
    auto makeHistoForSample = [varNumber, xname, fitsig, massbins, massmin, massmax, xbins, xmin, xmax, luminosity, reductionFactor](Sample* s)
    {
      // Output some info about the current file
      std::cout << std::endl;
      std::cout << Form("  /// Processing %s \n", s->name.Data());

      // cuts to apply to the events in the sample
      Run2MuonSelectionCuts     run2MuonSelection;
      Run2EventSelectionCuts80X run2EventSelectionData(true);

      // will probably want one for each type of mass
      ZCalibration* zcal_pf = new ZCalibration(xname, "mass_PF", fitsig, massmin, massmax, massbins, xmin, xmax, xbins);
      ZCalibration* zcal_roch = new ZCalibration(xname, "mass_Roch", fitsig, massmin, massmax, massbins, xmin, xmax, xbins);
      ZCalibration* zcal_kamu = new ZCalibration(xname, "mass_KaMu", fitsig, massmin, massmax, massbins, xmin, xmax, xbins);

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

      // need to return all the zcals
      return zcal_pf;
    };

 
    // Access plots in zcal pointer

    // Create the stack and ratio plot    
    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile("rootfiles/"+xname+"_data_8_0_X.root", "RECREATE");
    TDirectory* graphs = savefile->mkdir("graphs");
    TDirectory* hists = savefile->mkdir("hists");

    //// save the different histos in the appropriate directories in the tfile
    //hists->cd();
    //for(unsigned int i=0; i<zcal->histos.size(); i++)
    //    zcal->histos[i]->Write();

    //graphs->cd();
    //zcal->plot();
    //zcal.mean_vs_x->Write();
    //zcal.resolution_vs_x->Write();

    savefile->Close();

    return 0;
}
