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

void combineDataRunInfo(std::map<TString, ZCalibration*>& zcalmap)
{
    // get an item from the map, so that we have the zcal histo and plot settings
    ZCalibration* zfirst = 0;
    for(auto& zc: zcalmap)
        zfirst = zc.second;

    bool variableBinning = zfirst->variableBinning;

    std::cout << Form("  /// Running combineDataRunInfo \n");

    TString massname = "Net_Data_mass_";
    TString massname_pf = massname+"PF";
    TString massname_roch = massname+"Roch";
    TString massname_kamu = massname+"KaMu";

    ZCalibration* zcal_net_pf = 0;
    ZCalibration* zcal_net_roch = 0;
    ZCalibration* zcal_net_kamu = 0;

    if(!variableBinning)
    {
        zcal_net_pf = new ZCalibration(zfirst->xname, massname_pf, zfirst->fitsig, 
                                         zfirst->massmin, zfirst->massmax, zfirst->massbins, 
                                         zfirst->xmin, zfirst->xmax, zfirst->xbins);

        zcal_net_roch = new ZCalibration(zfirst->xname, massname_roch, zfirst->fitsig, 
                                         zfirst->massmin, zfirst->massmax, zfirst->massbins, 
                                         zfirst->xmin, zfirst->xmax, zfirst->xbins);

        zcal_net_kamu = new ZCalibration(zfirst->xname, massname_kamu, zfirst->fitsig, 
                                         zfirst->massmin, zfirst->massmax, zfirst->massbins, 
                                         zfirst->xmin, zfirst->xmax, zfirst->xbins);
    }
    else
    {
        zcal_net_pf = new ZCalibration(zfirst->xname, massname_pf, zfirst->fitsig, 
                                         zfirst->massmin, zfirst->massmax, zfirst->massbins, 
                                         zfirst->binning);

        zcal_net_roch = new ZCalibration(zfirst->xname, massname_roch, zfirst->fitsig, 
                                         zfirst->massmin, zfirst->massmax, zfirst->massbins, 
                                         zfirst->binning);

        zcal_net_kamu = new ZCalibration(zfirst->xname, massname_kamu, zfirst->fitsig, 
                                         zfirst->massmin, zfirst->massmax, zfirst->massbins, 
                                         zfirst->binning);
    }

    for(unsigned int i=0; i<zfirst->histos.size(); i++)
    {
        TList* add_hist_list_pf = new TList();
        TList* add_hist_list_roch = new TList();
        TList* add_hist_list_kamu = new TList();

        for(auto& zc: zcalmap)
        {
            if(!zc.second->massname.Contains("Run")) continue;

            if(zc.second->massname.Contains("PF")) add_hist_list_pf->Add(zc.second->histos[i]);
            else if(zc.second->massname.Contains("Roch")) add_hist_list_roch->Add(zc.second->histos[i]);
            else if(zc.second->massname.Contains("KaMu")) add_hist_list_kamu->Add(zc.second->histos[i]);
        }
 
        if(add_hist_list_pf->GetSize() != 0)
        {
            TH1D* net_hist_i_pf = DiMuPlottingSystem::addHists(add_hist_list_pf, zcal_net_pf->histos[i]->GetName(), zcal_net_pf->histos[i]->GetTitle()); 
            TH1D* temp_pf = zcal_net_pf->histos[i];
            zcal_net_pf->histos[i] = net_hist_i_pf;
            delete temp_pf;
        }

        if(add_hist_list_roch->GetSize() != 0)
        {
            TH1D* net_hist_i_roch = DiMuPlottingSystem::addHists(add_hist_list_roch, zcal_net_roch->histos[i]->GetName(), zcal_net_roch->histos[i]->GetTitle()); 
            TH1D* temp_roch = zcal_net_roch->histos[i];
            zcal_net_roch->histos[i] = net_hist_i_roch;
            delete temp_roch;
        }

        if(add_hist_list_kamu->GetSize() != 0)
        {
            TH1D* net_hist_i_kamu = DiMuPlottingSystem::addHists(add_hist_list_kamu, zcal_net_kamu->histos[i]->GetName(), zcal_net_kamu->histos[i]->GetTitle()); 
            TH1D* temp_kamu = zcal_net_kamu->histos[i];
            zcal_net_kamu->histos[i] = net_hist_i_kamu;
            delete temp_kamu;
        }
    }

    if(zcal_net_pf->histos.size() > 0 && zcal_net_pf->histos[0]->Integral() != 0) zcalmap[massname_pf] = zcal_net_pf;
    if(zcal_net_roch->histos.size() > 0 && zcal_net_roch->histos[0]->Integral() != 0) zcalmap[massname_roch] = zcal_net_roch;
    if(zcal_net_kamu->histos.size() > 0 && zcal_net_kamu->histos[0]->Integral() != 0) zcalmap[massname_kamu] = zcal_net_kamu;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

void initPlotSettings(int varNumber, float& fitsig, float& massmin, float& massmax, int& massbins, 
                      float& xmin, float& xmax, int& xbins, TString& xname, std::vector<Float_t>& binning)
{
    massmin = 82.2;
    massmax = 100.2;
    massbins = 62;
    xbins = 12;
    fitsig = 1;

    if(varNumber == 0)
    {
        xname = "phi_plus";
        xmin = -3.15;
        xmax = 3.15;
    }

    if(varNumber == 1)
    {
        xname = "phi_minus";
        xmin = -3.15;
        xmax = 3.15;
    }

    if(varNumber == 2)
    {
        xname = "eta_plus";
        xmin = -2.41;
        xmax = 2.41;
    }

    if(varNumber == 3)
    {
        xname = "eta_minus";
        xmin = -2.41;
        xmax = 2.41;
    }
    if(varNumber == 4)
    {
        xname = "pt_plus";
        xmin = 10;
        xmax = 210;
        xbins = 50;
        binning = {10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 510};
    }

    if(varNumber == 5)
    {
        xname = "pt_minus";
        xmin = 10;
        xmax = 510;
        xbins = 50;
        binning = {10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 510};
    }
    if(varNumber == 6)
    {
        xname = "dimu_pt";
        xmin = 0;
        xmax = 500;
        xbins=50;
        binning = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 500};
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

    float nthreads = 8;
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
    std::vector<Float_t> binning;

    initPlotSettings(varNumber, fitsig, massmin, massmax, massbins, xmin, xmax, xbins, xname, binning);

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
        if(i.second->sampleType != "data" && i.second->name != "ZJets_AMC") continue;
        //if(i.second->name != "RunE") continue;

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

    for(auto& i: binning)
        std::cout << i << ",";
    std::cout << std::endl << std::endl;


  
    auto makePlotsForSample = [&binning, varNumber, xname, fitsig, massbins, massmin, massmax, xbins, xmin, xmax, luminosity, reductionFactor](Sample* s)
    {
      // Output some info about the current file
      std::cout << Form("  /// Processing %s \n", s->name.Data());

      bool isData = s->sampleType == "data";

      // cuts to apply to the events in the sample
      Run2MuonSelectionCuts     run2MuonSelection;
      Run2EventSelectionCuts80X run2EventSelectionData(true);
      Run2EventSelectionCuts80X run2EventSelectionMC;

      // will probably want one for each type of mass
      ZCalibration* zcal_pf   = 0; 
      ZCalibration* zcal_roch =  0; 

      if(binning.size()==0) 
      {
          zcal_pf = new ZCalibration(xname, s->name+"_mass_PF", fitsig, massmin, massmax, massbins, xmin, xmax, xbins);
          zcal_roch = new ZCalibration(xname, s->name+"_mass_Roch", fitsig, massmin, massmax, massbins, xmin, xmax, xbins);
          //zcal_kamu = new ZCalibration(xname, s->name+"_mass_KaMu", fitsig, massmin, massmax, massbins, xmin, xmax, xbins);
      }
      else
      {
          zcal_pf = new ZCalibration(xname, s->name+"_mass_PF", fitsig, massmin, massmax, massbins, binning);
          zcal_roch = new ZCalibration(xname, s->name+"_mass_Roch", fitsig, massmin, massmax, massbins, binning);
          //zcal_kamu = new ZCalibration(xname, s->name+"_mass_KaMu", fitsig, massmin, massmax, massbins, binning);
      }

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
        if(s->vars.recoMuons->size() < 2) continue;
        bool found_good_dimuon = false;

        // find the first good dimuon candidate and fill info
        for(auto& dimu: (*s->vars.recoDimuCands))
        {
          s->vars.dimuCand = &dimu;

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
          if(!s->vars.recoMuons->at(dimu.iMu1).isTightID || !s->vars.recoMuons->at(dimu.iMu2).isTightID)
          { 
              continue; 
          }
          if(!run2MuonSelection.evaluate(s->vars)) 
          {
              continue; 
          }

          found_good_dimuon = true;

          MuonInfo& mu1 = s->vars.recoMuons->at(dimu.iMu1);
          MuonInfo& mu2 = s->vars.recoMuons->at(dimu.iMu2);

          // mc branches
          if(s->sampleType != "data")
          {
              s->branches.gen_wgt->GetEntry(i);
              s->branches.nPU->GetEntry(i);
              s->branches.pu_wgt->GetEntry(i);
              s->branches.eff_wgt->GetEntry(i);
          }

          if(varNumber == 0) 
          {
              float phi_plus = (mu1.charge==1)?mu1.phi:mu2.phi;
          
              zcal_pf->fill(phi_plus, dimu.mass_PF, s->getWeight());
              zcal_roch->fill(phi_plus, dimu.mass_Roch, s->getWeight());
              //zcal_kamu->fill(phi_plus, dimu.mass_KaMu, s->getWeight());
          }

          if(varNumber == 1) 
          {
              float phi_minus = (mu1.charge==-1)?mu1.phi:mu2.phi;

              zcal_pf->fill(phi_minus, dimu.mass_PF, s->getWeight());
              zcal_roch->fill(phi_minus, dimu.mass_Roch, s->getWeight());
              //zcal_kamu->fill(phi_minus, dimu.mass_KaMu, s->getWeight());
          }

          if(varNumber == 2) 
          {
              float eta_plus = (mu1.charge==1)?mu1.eta:mu2.eta;

              zcal_pf->fill(eta_plus, dimu.mass_PF, s->getWeight());
              zcal_roch->fill(eta_plus, dimu.mass_Roch, s->getWeight());
              //zcal_kamu->fill(eta_plus, dimu.mass_KaMu, s->getWeight());
          }

          if(varNumber == 3) 
          {
              float eta_minus = (mu1.charge==-1)?mu1.eta:mu2.eta;

              zcal_pf->fill(eta_minus, dimu.mass_PF, s->getWeight());
              zcal_roch->fill(eta_minus, dimu.mass_Roch, s->getWeight());
              //zcal_kamu->fill(eta_minus, dimu.mass_KaMu, s->getWeight());
          }

          if(varNumber == 4) 
          {
              float pt_plus_pf = (mu1.charge==1)?mu1.pt_PF:mu2.pt_PF;
              float pt_plus_roch = (mu1.charge==1)?mu1.pt_Roch:mu2.pt_Roch;
              //float pt_plus_kamu = (mu1.charge==1)?mu1.pt_KaMu:mu2.pt_KaMu;

              zcal_pf->fill(pt_plus_pf, dimu.mass_PF, s->getWeight());
              zcal_roch->fill(pt_plus_roch, dimu.mass_Roch, s->getWeight());
              //zcal_kamu->fill(pt_plus_kamu, dimu.mass_KaMu, s->getWeight());
          }

          if(varNumber == 5) 
          {
              float pt_minus_pf = (mu1.charge==-1)?mu1.pt_PF:mu2.pt_PF;
              float pt_minus_roch = (mu1.charge==-1)?mu1.pt_Roch:mu2.pt_Roch;
              //float pt_minus_kamu = (mu1.charge==-1)?mu1.pt_KaMu:mu2.pt_KaMu;

              zcal_pf->fill(pt_minus_pf, dimu.mass_PF, s->getWeight());
              zcal_roch->fill(pt_minus_roch, dimu.mass_Roch, s->getWeight());
              //zcal_kamu->fill(pt_minus_kamu, dimu.mass_KaMu, s->getWeight());
          }

          if(varNumber == 6) 
          {
              float dimu_pt_pf = dimu.pt_PF;
              float dimu_pt_roch = dimu.pt_Roch;
              //float dimu_pt_kamu = dimu.pt_KaMu;

              zcal_pf->fill(dimu_pt_pf, dimu.mass_PF, s->getWeight());
              zcal_roch->fill(dimu_pt_roch, dimu.mass_Roch, s->getWeight());
              //zcal_kamu->fill(dimu_pt_kamu, dimu.mass_KaMu, s->getWeight());
          }

          if(found_good_dimuon) break;

        } // end dimucand loop
      } // end events loop
      
      std::vector<ZCalibration*>* returnVector = new std::vector<ZCalibration*>();
      returnVector->push_back(zcal_pf);
      returnVector->push_back(zcal_roch);
      //returnVector->push_back(zcal_kamu);

      std::cout << Form("  /// Done processing %s \n", s->name.Data());
      return returnVector;

    };

   ///////////////////////////////////////////////////////////////////
   // SAMPLE PARALLELIZATION------ ----------------------------------
   ///////////////////////////////////////////////////////////////////

    ThreadPool pool(nthreads);
    std::vector< std::future< std::vector<ZCalibration*>* > > results;

    TStopwatch timerWatch;
    timerWatch.Start();

    for(auto &s : samplevec)
        results.push_back(pool.enqueue(makePlotsForSample, s));

    // put all of the zcalibration pointers into this map
    std::map<TString, ZCalibration*> zcalMap;  

    // Gather results for each sample
    for(auto && result: results)
    {
        // loop over each sample's PF, Roch, KaMu zcal
        for(auto zcal: (*result.get()))
        {
           zcalMap[zcal->massname] = zcal;
        }
    }

    // add up the individual data runs into net data and add to the map
    combineDataRunInfo(zcalMap);

    // Now get the plots from the map and organize them into lists.
    // Then, we can save them into the different root file directories easily.
    TList* overlays_list = new TList();
    TList* histos_list = new TList();
    TList* plots_list = new TList();
    std::map<TString, TList*> overlayMap;

    for(auto& item: zcalMap)
    {
        TString zcal_key = item.first;
        ZCalibration* zcal = item.second;
        std::cout << zcal_key << ": " << zcal->histos[0]->GetName() << ": " << zcal->histos[1]->Integral() << std::endl;
        zcal->plot();

        TGraph* mean_graph = zcal->mean_vs_x;
        TGraph* res_graph = zcal->resolution_vs_x;
        
        for(auto& h: zcal->histos)
            histos_list->Add(h);

        plots_list->Add(mean_graph);
        plots_list->Add(res_graph);

        // strip extraneous naming information to initialize the overlay map
        TString overlay_name = zcal->massname;
        overlay_name = overlay_name.ReplaceAll("_mass_PF", "");
        overlay_name = overlay_name.ReplaceAll("_mass_KaMu", "");
        overlay_name = overlay_name.ReplaceAll("_mass_Roch", "");

        // Name by mass type only
        //if(zcal->massname.Contains("PF")) 
        //{        
        //    mean_graph->SetTitle("PF");
        //    res_graph->SetTitle("PF");
        //}
        //else if(zcal->massname.Contains("Roch")) 
        //{        
        //    mean_graph->SetTitle("Roch");
        //    res_graph->SetTitle("Roch");
        //}
        //else if(zcal->massname.Contains("KaMu")) 
        //{        
        //    mean_graph->SetTitle("KaMu");
        //    res_graph->SetTitle("KaMu");
        //}

        // name by sample_MASSTYPE 
        TString title = zcal_key;
        title = title.ReplaceAll("_mass", "");

        mean_graph->SetTitle(title);
        res_graph->SetTitle(title);
       
        // overlay pf, roch, kamu mean for sample
        if(overlayMap.count(overlay_name+"_Mean") == 0) 
        {
            overlayMap[overlay_name+"_Mean"] = new TList();
            overlayMap[overlay_name+"_Mean"]->Add(mean_graph);
        }
        else 
        {
            overlayMap[overlay_name+"_Mean"]->Add(mean_graph);
        }

        // overlay pf, roch, kamu res for the sample
        if(overlayMap.count(overlay_name+"_Resolution") == 0) 
        {
            overlayMap[overlay_name+"_Resolution"] = new TList();
            overlayMap[overlay_name+"_Resolution"]->Add(res_graph);
        }
        else 
        {
            overlayMap[overlay_name+"_Resolution"]->Add(res_graph);
        }

        // overlay data and mc for pf
        if(zcal_key.Contains("PF") && (zcal_key.Contains("Net_Data") || zcal_key.Contains("ZJets")))
        {
            if(overlayMap.count("Data_MC_Mean_PF") == 0) 
            {
                overlayMap["Data_MC_Mean_PF"] = new TList();
                overlayMap["Data_MC_Mean_PF"]->Add((TGraph*)mean_graph->Clone());
            }
            else 
            {
                overlayMap["Data_MC_Mean_PF"]->Add((TGraph*)mean_graph->Clone());
            }

            if(overlayMap.count("Data_MC_Resolution_PF") == 0) 
            {
                overlayMap["Data_MC_Resolution_PF"] = new TList();
                overlayMap["Data_MC_Resolution_PF"]->Add((TGraph*)res_graph->Clone());
            }
            else 
            {
                overlayMap["Data_MC_Resolution_PF"]->Add((TGraph*)res_graph->Clone());
            }
        }
        // overlay data and mc for roch
        if(zcal_key.Contains("Roch") && (zcal_key.Contains("Net_Data") || zcal_key.Contains("ZJets")))
        {
            if(overlayMap.count("Data_MC_Mean_Roch") == 0) 
            {
                overlayMap["Data_MC_Mean_Roch"] = new TList();
                overlayMap["Data_MC_Mean_Roch"]->Add((TGraph*)mean_graph->Clone());
            }
            else 
            {
                overlayMap["Data_MC_Mean_Roch"]->Add((TGraph*)mean_graph->Clone());
            }

            if(overlayMap.count("Data_MC_Resolution_Roch") == 0) 
            {
                overlayMap["Data_MC_Resolution_Roch"] = new TList();
                overlayMap["Data_MC_Resolution_Roch"]->Add((TGraph*)res_graph->Clone());
            }
            else 
            {
                overlayMap["Data_MC_Resolution_Roch"]->Add((TGraph*)res_graph->Clone());
            }
        }
        // overlay data and mc for kamu
        if(zcal_key.Contains("KaMu") && (zcal_key.Contains("Net_Data") || zcal_key.Contains("ZJets")))
        {
            if(overlayMap.count("Data_MC_Mean_KaMu") == 0) 
            {
                overlayMap["Data_MC_Mean_KaMu"] = new TList();
                overlayMap["Data_MC_Mean_KaMu"]->Add((TGraph*)mean_graph->Clone());
            }
            else 
            {
                overlayMap["Data_MC_Mean_KaMu"]->Add((TGraph*)mean_graph->Clone());
            }

            if(overlayMap.count("Data_MC_Resolution_KaMu") == 0) 
            {
                overlayMap["Data_MC_Resolution_KaMu"] = new TList();
                overlayMap["Data_MC_Resolution_KaMu"]->Add((TGraph*)res_graph->Clone());
            }
            else 
            {
                overlayMap["Data_MC_Resolution_KaMu"]->Add((TGraph*)res_graph->Clone());
            }
        }


    }

    // overlay the mean and resolution vs x plots
    for(auto& item: overlayMap)
    {   
        TString yname;

        if(item.first.Contains("Mean"))
            yname = "Fit Mean (GeV)";
        else
            yname = "Fit Resolution (GeV)";

        TCanvas* c = DiMuPlottingSystem::overlay(item.second, item.first, item.first, xname, yname, false);
        overlays_list->Add(c);
    }   

    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile("rootfiles/zcalibration_"+xname+"_data_8_0_X.root", "RECREATE");
    savefile->cd();
    TDirectory* overlaydir = savefile->mkdir("overlays");
    TDirectory* plotdir = savefile->mkdir("plots");
    TDirectory* histodir = savefile->mkdir("histos");

    // Create the stack and ratio plot    
    overlaydir->cd();
    overlays_list->Write();

    plotdir->cd();
    plots_list->Write();
 
    histodir->cd();
    histos_list->Write();

    savefile->Write();
    savefile->Close();

    timerWatch.Stop();
    std::cout << "### DONE " << timerWatch.RealTime() << " seconds" << std::endl;

    return 0;
}
