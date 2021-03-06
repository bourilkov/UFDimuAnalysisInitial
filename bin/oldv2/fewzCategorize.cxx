// Make DY_MC or Data plots to compare DY_MC/Data to FEWZ

#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "CutSet.h"
#include "Cut.h"
#include "SelectionCuts.h"
#include "CategorySelection.h"
#include "JetSelectionTools.h"

#include "EventTools.h"
#include "PUTools.h"
#include "ParticleTools.h"
#include "SignificanceMetrics.cxx"

#include "TLorentzVector.h"

#include <sstream>
#include <map>
#include <vector>
#include <utility>

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    int input = 0;
    bool useRecoMuCuts = 0;
    bool useRecoToPlotMu = 0;
    bool useRecoJetCuts = 1;
    bool useRecoToPlotJets = 1;
    bool cutNoGens = 1;
    bool plotData = 1;

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) ss >> input;
        if(i==2) ss >> useRecoMuCuts;
        if(i==3) ss >> useRecoToPlotMu;
        if(i==4) ss >> useRecoJetCuts;
        if(i==5) ss >> useRecoToPlotJets;
        if(i==6) ss >> cutNoGens;
        if(i==7) ss >> plotData;
    }   

    float reductionFactor = 1;     // run over less data (e.g. 20 = 20x less). useful for debugging, runs faster.
    float luminosity = 33598;      // pb-1
    float triggerSF = 0.913;       // no HLT trigger info available for the samples so we scale for the trigger efficiency instead
    float signalSF = 100;          // not using this at the moment

    // map of samples
    std::map<std::string, Sample*> samples;

    // Second container so that we can have a copy sorted by cross section.
    std::vector<Sample*> samplevec;

    // Use this to plot some things if we wish
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

  
    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // ================================================================
    // Data -----------------------------------------------------------
    // ================================================================

    if(plotData)
    {
        TString datafilename = 
        TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data/25ns/golden/CMSSW_8_0_X/stage_1_singleMuon_Run2016BCDEFGH_ALL_etm.root");

        Sample* datasample = new Sample(datafilename, "Data", "data");
        datasample->lumi = luminosity;
        datasample->xsec = 9999;
        datasample->pileupfile = "pu_reweight_trees/8_0_X/PU_2016BCDEFGH_xsec69p2mb_CMSSW_8_0_X.root";
        samples["Data"] = datasample;
    }

    // ================================================================
    // DYJetsToLL -----------------------------------------------------
    // ================================================================

    TString dyfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/bg/dy/CMSSW_8_0_X/stage_1_dy_jetsToLL_ALL_etm.root");
    samples["DYJetsToLL"] = new Sample(dyfilename, "DYJetsToLL", "background");
    samples["DYJetsToLL"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_DYJetsToLL.root"; //nPU
    samples["DYJetsToLL"]->xsec = 6025.2; // pb

    ///////////////////////////////////////////////////////////////////
    // PREPROCESSING---------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // Loop through all of the samples to do some pre-processing
    std::cout << std::endl;
    std::cout << "======== Preprocess the samples... " << std::endl;
    std::cout << std::endl;
    
    for(auto &i : samples)
    {
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

// no pileup reweighting for this study. commenting this out.
/*
        if(!i.second->sampleType.Contains("data"))
        {
            // Pileup reweighting
            std::cout << "    +++ PU Reweighting " << i.second->name << "..."  << std::endl;
            std::cout << std::endl;

            i.second->lumiWeights = new reweight::LumiReWeighting(i.second->pileupfile.Data(), samples["Data"]->pileupfile.Data(), "pileup", "pileup");
            std::cout << "        " << i.first << "->lumiWeights: " << i.second->lumiWeights << std::endl;
            std::cout << std::endl;
        }
*/
        samplevec.push_back(i.second);
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 
    

    ///////////////////////////////////////////////////////////////////
    // Cut and Categorize ---------------------------------------------
    ///////////////////////////////////////////////////////////////////
    
    TString varname;
    int nbins, wbins;
    float nmin, wmin;
    float nmax, wmax;

    // dimu_mass
    if(input == 0)
    {
        nbins = 50;
        nmin = 110;
        nmax = 160;

        wbins = 100;
        wmin = 110;
        wmax = 310;

        varname = "dimu_mass";
    }

    // dimu_pt 
    if(input == 1)
    {
        nbins = 200;
        nmin = 0;
        nmax = 100;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "dimu_pt";
    }

    // recoMu_pt
    if(input == 2)
    {
        nbins = 200;
        nmin = 0;
        nmax = 150;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "mu_pt";
    }
 
    // recoMu_eta
    if(input == 3)
    {
        nbins = 100;
        nmin = -2.5;
        nmax = 2.5;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "mu_eta";
    }

    // NPV
    if(input == 4)
    {
        nbins = 50;
        nmin = 0;
        nmax = 50;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "NPV";
    }

    // jet_pt
    if(input == 5)
    {   
        nbins = 200;
        nmin = 0;
        nmax = 200;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "jet_pt";
    }   

    // jet_eta 
    if(input == 6)
    {   
        nbins = 100;
        nmin = -5; 
        nmax = 5;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "jet_eta";
    }   

    // nValJets
    if(input == 7)
    {   
        nbins = 11;
        nmin = 0; 
        nmax = 11;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "nValJets";
    }   

    // m_jj
    if(input == 8)
    {   
        nbins = 200;
        nmin = 0; 
        nmax = 2000;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "m_jj";
    }   

    // dEta_jj
    if(input == 9)
    {   
        nbins = 100;
        nmin = -10; 
        nmax = 10;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "dEta_jj";
    }   

    // dR_jj
    if(input == 10)
    {   
        nbins = 100;
        nmin = 0; 
        nmax = 5;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "dR_jj";
    }   

    // Figure out the save file name here, since data may change the flags for the settings later.
    TString savefilename = TString("rootfiles/")+Form("%d%d%d%d%d%d_", (int)useRecoMuCuts, (int)useRecoToPlotMu, (int)useRecoJetCuts,
                                 (int)useRecoToPlotJets, (int)cutNoGens, (int)plotData)+ "validate_"+varname+"_DY-FEWZ_MC_categories_"+
                                 Form("%d",(int)luminosity)+".root";


    std::cout << std::endl;
    std::cout << "======== Plot Configs ========" << std::endl;
    std::cout << "useRecoMuCuts    : " << useRecoMuCuts << std::endl;
    std::cout << "useRecoToPlotMu  : " << useRecoToPlotMu << std::endl;
    std::cout << "useRecoJetCuts   : " << useRecoJetCuts << std::endl;
    std::cout << "useRecoToPlotJets: " << useRecoToPlotJets << std::endl;
    std::cout << "cutNoGens        : " << cutNoGens << std::endl;
    std::cout << std::endl;
    std::cout << "var          : " << varname << std::endl;
    std::cout << "nmin         : " << nmin << std::endl;
    std::cout << "nmax         : " << nmax << std::endl;
    std::cout << "nbins        : " << nbins << std::endl;
    std::cout << std::endl;
    std::cout << "wmin         : " << wmin << std::endl;
    std::cout << "wmax         : " << wmax << std::endl;
    std::cout << "wbins        : " << wbins << std::endl;
    std::cout << std::endl;


    // Objects to help with the cuts and selections
    JetSelectionTools jetSelectionTools;
    FEWZCompareCuts FEWZselection(useRecoMuCuts);
    CategorySelectionFEWZ categorySelection(useRecoMuCuts, useRecoJetCuts);

    for(auto &s : samplevec)
    {
      // Output some info about the current file
      std::cout << std::endl;
      std::cout << "  /// Looping over " << s->name << std::endl;

      ///////////////////////////////////////////////////////////////////
      // HISTOGRAMS TO FILL ---------------------------------------------
      ///////////////////////////////////////////////////////////////////

      // Keep track of which histogram (sample and variable) to fill in the category
      TString hkey;

      // We usually fill the same plots for each category, just with different events
      // so we loop through each category and fill the appropriate events for each category's histo
      for(auto &c : categorySelection.categoryMap)
      {
          //number of bins for the histogram, initialize to narrow category binning
          //overwrite later if it's a wide category
          int hbins = nbins;
          double hmin = nmin;
          double hmax = nmax;

          // c.second is the category object, c.first is the category name
          TString hname = c.first+"_"+s->name;
          hkey = s->name;

          // The VBF categories have low stats so we use fewer bins
          // Same goes for 01jet categories with dijet variables
          if(c.first.Contains("Wide")) 
          {
              hbins = wbins; 
              hmin = wmin;
              hmax = wmax;
          }

          // Set up the histogram for the category and variable to plot
          c.second.histoMap[hkey] = new TH1D(hname, hname, hbins, hmin, hmax);
          c.second.histoMap[hkey]->GetXaxis()->SetTitle(varname);
          c.second.histoList->Add(c.second.histoMap[hkey]); // need them ordered by xsec for the stack and ratio plot
      }

      if(s->sampleType.Contains("data"))
      {
          // must use reco information for data
          // data is last so we can overwrite these and it won't affect the drell yan sample
          useRecoMuCuts = 1;
          useRecoToPlotMu = 1;
          useRecoJetCuts = 1;
          useRecoToPlotJets = 1;
          cutNoGens = 0;

          categorySelection.useRecoMu   = useRecoMuCuts;
          categorySelection.useRecoJets = useRecoJetCuts;
          FEWZselection.useReco = useRecoMuCuts; 
      }

      for(unsigned int i=0; i<s->N/reductionFactor; i++)
      {

        ///////////////////////////////////////////////////////////////////
        // GET INFORMATION ------------------------------------------------
        ///////////////////////////////////////////////////////////////////

        s->getEntry(i); 

        s->vars.validJets = std::vector<TLorentzVector>();
        s->vars.validGenJets = std::vector<TLorentzVector>();

        jetSelectionTools.getValidJetsdR(s->vars, s->vars.validJets);

        // Can't get the MC info until we know the sample is MC
        _TrackInfo gen_mu0;
        _TrackInfo gen_mu1;
        TLorentzVector gen_dimu;

        std::pair<int,int> e(s->vars.eventInfo.run, s->vars.eventInfo.event); // create a pair that identifies the event uniquely

        // can get the gen info since we have a MC sample
        if(!s->sampleType.Contains("data"))
        {
            jetSelectionTools.getValidGenValJets(s->vars, s->vars.validGenJets);
            // get first postFSR DY gen muon
            gen_mu0 = ParticleTools::getGenMuDY(0, 1, s->vars);
            // get second postFSR DY gen muon
            gen_mu1 = ParticleTools::getGenMuDY(1, 1, s->vars);
            // get dimuon candidate from DY gen muons
            gen_dimu = ParticleTools::getMotherPtEtaPhiM(gen_mu0.pt, gen_mu0.eta, gen_mu0.phi, MASS_MUON, gen_mu1.pt, gen_mu1.eta, gen_mu1.phi, MASS_MUON);

            if(gen_mu0.pt < 0 || gen_mu1.pt < 0 && cutNoGens)
            {   
                continue;
                // was printing this for debugging, don't want to print it right now, but leaving it here in case I do later
                std::cout << "run, lumi, event" << std::endl;
                std::cout << s->vars.eventInfo.run << ", " << s->vars.eventInfo.lumi << ", " << s->vars.eventInfo.event << std::endl;
                std::cout << "gen muons gathered" << std::endl;
                std::cout << gen_mu0.pt << "," << gen_mu0.eta << "," << gen_mu0.phi << std::endl;
                std::cout << gen_mu1.pt << "," << gen_mu1.eta << "," << gen_mu1.phi << std::endl;
                std::cout << "gen muons from Z" << std::endl;
                std::cout << s->vars.genM1ZpostFSR.pt << "," << s->vars.genM1ZpostFSR.eta << "," << s->vars.genM1ZpostFSR.phi << std::endl;
                std::cout << s->vars.genM2ZpostFSR.pt << "," << s->vars.genM2ZpostFSR.eta << "," << s->vars.genM2ZpostFSR.phi << std::endl;
                std::cout << "gen muons from gamma*" << std::endl;
                std::cout << s->vars.genM1GpostFSR.pt << "," << s->vars.genM1GpostFSR.eta << "," << s->vars.genM1GpostFSR.phi << std::endl;
                std::cout << s->vars.genM2GpostFSR.pt << "," << s->vars.genM2GpostFSR.eta << "," << s->vars.genM2GpostFSR.phi << std::endl;
                std::cout << "Dimuon gen mother" << std::endl;
                std::cout << s->vars.genZpostFSR.pt << "," << s->vars.genZpostFSR.eta << "," << s->vars.genZpostFSR.phi << std::endl;
                std::cout << s->vars.genGpostFSR.pt << "," << s->vars.genGpostFSR.eta << "," << s->vars.genGpostFSR.phi << std::endl;
                std::cout << std::endl;
            }
        }

        ///////////////////////////////////////////////////////////////////
        // CUTS  ----------------------------------------------------------
        ///////////////////////////////////////////////////////////////////
        
        // If cutNoGens is true, we don't want to plot events where there are no gen muons from a Z or gamma*
        
        // don't apply iso cuts when using reco, already turned off for gen
        //FEWZselection.cutset.cuts[7].on = false;
        //FEWZselection.cutset.cuts[8].on = false;

        //if(!s->vars.muons.isTightMuon[0] || !s->vars.muons.isTightMuon[1] && useRecoMuCuts)
        if(!s->vars.muons.isTightMuon[0] || !s->vars.muons.isTightMuon[1])
        {
            continue; 
        }
        if(!FEWZselection.evaluate(s->vars)) 
        {
            continue; 
        }

        // Figure out which category the event belongs to
        categorySelection.evaluate(s->vars);

        // Look at each category
        for(auto &c : categorySelection.categoryMap)
        {
            // dimuCand.recoCandMass
            if(varname.EqualTo("dimu_mass")) 
            {
                float varvalue = useRecoToPlotMu?s->vars.dimuCand.recoCandMassPF:gen_dimu.M();
                // blind the signal region for data but not for MC
                if(!(s->sampleType.Contains("data") && varvalue >= 110 && varvalue < 140))
                    // if the event is in the current category then fill the category's histogram for the given sample and variable
                    if(c.second.inCategory) c.second.histoMap[hkey]->Fill(varvalue, s->getWeight());
                    //std::cout << "    " << c.first << ": " << varvalue;
            }

            if(varname.EqualTo("dimu_pt"))
            {
                // if the event is in the current category then fill the category's histogram for the given sample and variable
                if(c.second.inCategory) c.second.histoMap[hkey]->Fill(useRecoToPlotMu?s->vars.dimuCand.recoCandPtPF:gen_dimu.Pt(), s->getWeight());
            }

            if(varname.EqualTo("mu_pt"))
            {
                if(c.second.inCategory)
                {
                    c.second.histoMap[hkey]->Fill(useRecoToPlotMu?s->vars.muons.pt[0]:gen_mu0.pt, s->getWeight());
                    c.second.histoMap[hkey]->Fill(useRecoToPlotMu?s->vars.muons.pt[1]:gen_mu1.pt, s->getWeight());
                }
            }

            // recoMu_Eta
            if(varname.EqualTo("mu_eta"))
            {
                if(c.second.inCategory) c.second.histoMap[hkey]->Fill(useRecoToPlotMu?s->vars.muons.eta[0]:gen_mu0.eta, s->getWeight());
                if(c.second.inCategory) c.second.histoMap[hkey]->Fill(useRecoToPlotMu?s->vars.muons.eta[1]:gen_mu1.eta, s->getWeight());
            }

            // NPV
            if(varname.EqualTo("NPV"))
            {
                 if(c.second.inCategory) c.second.histoMap[hkey]->Fill(useRecoToPlotMu?s->vars.vertices.nVertices:s->vars.nPU, s->getWeight());
            }

            // jet_pt all valid
            if(varname.EqualTo("jet_pt"))
            {
                if(c.second.inCategory && useRecoToPlotJets) 
                {
                    for(unsigned int j=0; j<s->vars.validJets.size(); j++)
                        c.second.histoMap[hkey]->Fill(s->vars.validJets[j].Pt(), s->getWeight());
                }
                if(c.second.inCategory && !useRecoToPlotJets) 
                {
                    for(unsigned int j=0; j<s->vars.validGenJets.size(); j++)
                        c.second.histoMap[hkey]->Fill(s->vars.validGenJets[j].Pt(), s->getWeight());
                }
            }

            // jet_eta all valid
            if(varname.EqualTo("jet_eta"))
            {
                if(c.second.inCategory && useRecoToPlotJets) 
                {
                    for(unsigned int j=0; j<s->vars.validJets.size(); j++)
                        c.second.histoMap[hkey]->Fill(s->vars.validJets[j].Eta(), s->getWeight());
                }
                if(c.second.inCategory && !useRecoToPlotJets) 
                {
                    for(unsigned int j=0; j<s->vars.validGenJets.size(); j++)
                        c.second.histoMap[hkey]->Fill(s->vars.validGenJets[j].Eta(), s->getWeight());
                }
            }

            // nValJets
            if(varname.EqualTo("nValJets"))
            {
                 if(c.second.inCategory && useRecoToPlotJets) c.second.histoMap[hkey]->Fill(s->vars.validJets.size(), s->getWeight());
                 if(c.second.inCategory && !useRecoToPlotJets) c.second.histoMap[hkey]->Fill(s->vars.validGenJets.size(), s->getWeight());
            }

            // m_jj leading two
            if(varname.EqualTo("m_jj"))
            {
                 if(s->vars.validJets.size() >= 2 && useRecoToPlotJets)
                 {
                     TLorentzVector dijet = s->vars.validJets[0] + s->vars.validJets[1];
                     if(c.second.inCategory) c.second.histoMap[hkey]->Fill(dijet.M(), s->getWeight());
                 }
                 if(s->vars.validGenJets.size() >= 2 && !useRecoToPlotJets)
                 {
                     TLorentzVector dijet = s->vars.validGenJets[0] + s->vars.validGenJets[1];
                     if(c.second.inCategory) c.second.histoMap[hkey]->Fill(dijet.M(), s->getWeight());
                 }
            }

            // dEta_jj leading two
            if(varname.EqualTo("dEta_jj"))
            {
                 if(s->vars.validJets.size() >= 2 && useRecoToPlotJets)
                 {
                     float dEta = s->vars.validJets[0].Eta() - s->vars.validJets[1].Eta();
                     if(c.second.inCategory) c.second.histoMap[hkey]->Fill(dEta, s->getWeight());
                 }
                 if(s->vars.validGenJets.size() >= 2 && !useRecoToPlotJets)
                 {
                     float dEta = s->vars.validGenJets[0].Eta() - s->vars.validGenJets[1].Eta();
                     if(c.second.inCategory) c.second.histoMap[hkey]->Fill(dEta, s->getWeight());
                 }
            }

            // dR_jj
            if(varname.EqualTo("dR_jj"))
            {
                 if(s->vars.validJets.size() >= 2 && useRecoToPlotJets)
                 {
                     for(unsigned int j=0; j<s->vars.validJets.size(); j++)
                     {                
                       for(unsigned int k=j+1; k<s->vars.validJets.size(); k++)
                       {                
                           float dR = JetSelectionTools::dR(s->vars.validJets[j].Eta(), s->vars.validJets[j].Phi(), 
                                                            s->vars.validJets[k].Eta(),  s->vars.validJets[k].Phi());
                           if(c.second.inCategory) c.second.histoMap[hkey]->Fill(dR, s->getWeight());
                       }
                     }
                 }
                 if(s->vars.validGenJets.size() >= 2 && !useRecoToPlotJets)
                 {
                     for(unsigned int j=0; j<s->vars.validGenJets.size(); j++)
                     {                
                       for(unsigned int k=j+1; k<s->vars.validGenJets.size(); k++)
                       {                
                           float dR = JetSelectionTools::dR(s->vars.validGenJets[j].Eta(), s->vars.validGenJets[j].Phi(), 
                                                            s->vars.validGenJets[k].Eta(),  s->vars.validGenJets[k].Phi());
                           if(c.second.inCategory) c.second.histoMap[hkey]->Fill(dR, s->getWeight());
                       }
                     }
                 }
            }

        } // end category loop

        if(false)
          // ouput pt, mass info etc for the event
          EventTools::outputEvent(s->vars, categorySelection);

        // Reset the flags in preparation for the next event
        categorySelection.reset();

      } // end event loop

      // Scale according to luminosity and sample xsec now that the histograms are done being filled for that sample
      for(auto &c : categorySelection.categoryMap)
      {
          c.second.histoMap[hkey]->Scale(s->getScaleFactor(luminosity));
      }

    } // end sample loop

    TList* histolist = new TList();     // list to save all of the histos

    // For each category we want to compare the dimu_mass for DY/Data with the dimu mass for FEWZ
    for(auto &c : categorySelection.categoryMap)
    {
        histolist->Add(c.second.histoList);
    }
   
    //std::cout << std::endl;

    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile(savefilename, "RECREATE");

    TDirectory* histos = savefile->mkdir("histos");

    // save the different histos and stacks in the appropriate directories in the tfile
    histos->cd();
    histolist->Write();
    savefile->Close();
    return 0;
}
