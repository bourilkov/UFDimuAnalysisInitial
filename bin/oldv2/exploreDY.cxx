// Used to debug some stuff in DY_MC since there were discrepancies between
// DY_MC and fewz in the 1jet category during initial comparisons.
// Check out gen muons from Z/gamma*, plot variables and output info to terminal
// for all of the objects when plotting jet_pt or jet_eta.
// Wanted to debug the jet discrepancies hence the printouts for those.
// Turns out that that fewz emulates antikt jets so we need to compare fewz jets to
// DY_MC reco jets.


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
    bool useReco = 0;

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) ss >> input;
        if(i==2) ss >> useReco;
    }   

    // Not sure that we need a map if we have a vector
    // Should use this as the main database and choose from it to make the vector
    std::map<std::string, Sample*> samples;

    // Second container so that we can have a copy sorted by cross section.
    std::vector<Sample*> samplevec;

    // Use this to plot some things if we wish
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

  
    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    float luminosity = 3990;       // pb-1
    //float luminosity = 12900;      // pb-1 (ICHEP)
    //float luminosity = 15900;      // pb-1   (2016-08-03)
    float triggerSF = 0.913;       // no HLT trigger info available for the samples so we scale for the trigger efficiency instead
    float signalSF = 100;

    // ================================================================
    // Data -----------------------------------------------------------
    // ================================================================
/* 

    TString datafilename = 
    TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data/25ns/golden/CMSSW_8_0_X/old/stage_1_singleMuon_Run2016B_ALL.root");
    //TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data/25ns/golden/CMSSW_8_0_X/stage_1_singleMuon_Run2016BCD_ICHEP_ALL.root");
    //TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data/25ns/golden/CMSSW_8_0_X/stage_1_singleMuon_Run2016BCDE_ALL.root");

    Sample* datasample = new Sample(datafilename, "Data", "data");
    datasample->lumi = luminosity;
    datasample->xsec = 9999;
    datasample->pileupfile = "pu_reweight_trees/8_0_X/PU_2016B_xsec69mb_CMSSW_8_0_X.root";
    //datasample->pileupfile = "pu_reweight_trees/8_0_X/PU_2016BCD_ICHEP_xsec69p2mb_CMSSW_8_0_X.root";
    //datasample->pileupfile = "pu_reweight_trees/8_0_X/PU_2016BCDE_xsec69p2mb_CMSSW_8_0_X.root";
    samples["Data"] = datasample;
*/
    // ================================================================
    // DYJetsToLL -----------------------------------------------------
    // ================================================================

    TString dyfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/bg/dy/CMSSW_8_0_X/stage_1_dy_jetsToLL_bonusMuons_ALL.root");
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
    
    // Objects to help with the cuts and selections
    JetSelectionTools jetSelectionTools;
    FEWZCompareCuts FEWZselection(useReco);
    CategorySelectionFEWZ categorySelection(useReco);

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

    // mu_pt
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
 
    // mu_eta
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

    // N__jets
    if(input == 7)
    {   
        nbins = 11;
        nmin = 0; 
        nmax = 11;

        wbins = nbins;
        wmin = nmin;
        wmax = nmax;

        varname = "N_jets";
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

    std::cout << std::endl;
    std::cout << "======== Plot Configs ========" << std::endl;
    std::cout << "var          : " << varname << std::endl;
    std::cout << "nmin         : " << nmin << std::endl;
    std::cout << "nmax         : " << nmax << std::endl;
    std::cout << "nbins        : " << nbins << std::endl;
    std::cout << std::endl;
    std::cout << "wmin         : " << wmin << std::endl;
    std::cout << "wmax         : " << wmax << std::endl;
    std::cout << "wbins        : " << wbins << std::endl;
    std::cout << std::endl;

    TList* histolist = new TList();

    // Not sure how to deal with the scaling correctly when using a subset of events
    float reductionFactor = 1;

    // only sample is DY, so it's not really a loop
    for(auto &s : samplevec)
    {
      // Output some info about the current file
      std::cout << std::endl;
      std::cout << "  /// Looping over " << s->name << std::endl;

      ///////////////////////////////////////////////////////////////////
      // HISTOGRAMS TO FILL ---------------------------------------------
      ///////////////////////////////////////////////////////////////////

      // Set up the histogram for the category and variable to plot
      TString hname = TString(useReco?"":"gen_")+varname+"_dy";
      TH1F* hist = new TH1F(hname, hname, wbins, wmin, wmax);
      histolist->Add(hist);

      for(unsigned int i=0; i<s->N/reductionFactor; i++)
      {

        ///////////////////////////////////////////////////////////////////
        // GET INFORMATION ------------------------------------------------
        ///////////////////////////////////////////////////////////////////

        s->getEntry(i); 

        s->vars.validJets = std::vector<TLorentzVector>();
        s->vars.validGenValJets = std::vector<TLorentzVector>();

        jetSelectionTools.getValidJetsdR(s->vars, s->vars.validJets);
        jetSelectionTools.getValidGenValJets(s->vars, s->vars.validGenValJets);

        std::pair<int,int> e(s->vars.eventInfo.run, s->vars.eventInfo.event); // create a pair that identifies the event uniquely

        // get first postFSR DY gen muon
        _TrackInfo gen_mu0 = ParticleTools::getGenMuDY(0, 1, s->vars);
        // get second postFSR DY gen muon
        _TrackInfo gen_mu1 = ParticleTools::getGenMuDY(1, 1, s->vars);
        // get dimuon candidate from DY gen muons
        TLorentzVector gen_dimu = ParticleTools::getMotherPtEtaPhiM(gen_mu0.pt, gen_mu0.eta, gen_mu0.phi, MASS_MUON, gen_mu1.pt, gen_mu1.eta, gen_mu1.phi, MASS_MUON);

        if(gen_mu0.pt < 0 || gen_mu1.pt < 0)
        {
            continue;
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

        ///////////////////////////////////////////////////////////////////
        // CUTS  ----------------------------------------------------------
        ///////////////////////////////////////////////////////////////////
        
        // don't apply iso cuts when using reco, already turned off for gen
        FEWZselection.cutset.cuts[7].on = false;
        FEWZselection.cutset.cuts[8].on = false;

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

        // recoCandMass
        if(varname.EqualTo("dimu_mass")) 
        {
            float varvalue = useReco?s->vars.recoCandMassPF:gen_dimu.M();
            // blind the signal region for data but not for MC
            if(!(s->sampleType.Contains("data") && varvalue >= 110 && varvalue < 140))
                // if the event is in the current category then fill the category's histogram for the given sample and variable
                hist->Fill(varvalue, s->getWeight());
        }

        if(varname.EqualTo("dimu_pt"))
        {
            // if the event is in the current category then fill the category's histogram for the given sample and variable
            hist->Fill(useReco?s->vars.recoCandPtPF:gen_dimu.Pt(), s->getWeight());
        }

        if(varname.EqualTo("mu_pt"))
        {
           hist->Fill(useReco?s->vars.muons.pt[0]:gen_mu0.pt, s->getWeight());
           hist->Fill(useReco?s->vars.muons.pt[1]:gen_mu1.pt, s->getWeight());
        }

        // recoMu_Eta
        if(varname.EqualTo("mu_eta"))
        {
            hist->Fill(useReco?s->vars.muons.eta[0]:gen_mu0.eta, s->getWeight());
            hist->Fill(useReco?s->vars.muons.eta[1]:gen_mu1.eta, s->getWeight());
        }

        // NPV
        if(varname.EqualTo("NPV"))
        {
             hist->Fill(useReco?s->vars.vertices.nVertices:s->vars.nPU, s->getWeight());
        }

        // jet_pt
        if(varname.EqualTo("jet_pt"))
        {
            if(useReco) 
            {
                if(i<2000) EventTools::outputEvent(s->vars, categorySelection);
                for(unsigned int j=0; j<s->vars.jets.nValJets && j<N_JET_INFO; j++)
                    hist->Fill(s->vars.jets.pt[j], s->getWeight());
            }
            if(!useReco) 
            {
                if(i<2000) EventTools::outputEvent(s->vars, categorySelection);
                for(unsigned int j=0; j<s->vars.genValJets.nValJets && j<N_JET_INFO; j++)
                    hist->Fill(s->vars.genValJets.pt[j], s->getWeight());
            }
        }

        // jet_eta
        if(varname.EqualTo("jet_eta"))
        {
            if(useReco) 
            {
                if(i<2000) EventTools::outputEvent(s->vars, categorySelection);
                for(unsigned int j=0; j<s->vars.jets.nValJets && j<N_JET_INFO; j++)
                    hist->Fill(s->vars.jets.eta[j], s->getWeight());
            }
            if(!useReco) 
            {
                if(i<2000) EventTools::outputEvent(s->vars, categorySelection);
                for(unsigned int j=0; j<s->vars.genValJets.nValJets && j<N_JET_INFO; j++)
                    hist->Fill(s->vars.genValJets.eta[j], s->getWeight());
            }
        }

        // N_jets
        if(varname.EqualTo("N_jets"))
        {
             if(useReco)  hist->Fill(s->vars.jets.nValJets, s->getWeight());
             if(!useReco) hist->Fill(s->vars.genValJets.nValJets, s->getWeight());
        }

        // m_jj
        if(varname.EqualTo("m_jj"))
        {
             if(s->vars.validJets.size() >= 2 && useReco)
             {
                 TLorentzVector dijet = s->vars.validJets[0] + s->vars.validJets[1];
                 hist->Fill(dijet.M(), s->getWeight());
             }
             if(s->vars.validGenValJets.size() >= 2 && !useReco)
             {
                 TLorentzVector dijet = s->vars.validGenValJets[0] + s->vars.validGenValJets[1];
                 hist->Fill(dijet.M(), s->getWeight());
             }
        }

        // dEta_jj
        if(varname.EqualTo("dEta_jj"))
        {
             if(s->vars.validJets.size() >= 2 && useReco)
             {
                 float dEta = s->vars.validJets[0].Eta() - s->vars.validJets[1].Eta();
                 hist->Fill(dEta, s->getWeight());
             }
             if(s->vars.validGenValJets.size() >= 2 && !useReco)
             {
                 float dEta = s->vars.validGenValJets[0].Eta() - s->vars.validGenValJets[1].Eta();
                 hist->Fill(dEta, s->getWeight());
             }
        }

        if(false)
          // ouput pt, mass info etc for the event
          EventTools::outputEvent(s->vars, categorySelection);

        // Reset the flags in preparation for the next event
        categorySelection.reset();

      } // end event loop

      // Scale according to luminosity and sample xsec now that the histograms are done being filled for that sample
      hist->Scale(s->getScaleFactor(luminosity));

    } // end sample loop
   
    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile(TString("rootfiles/")+TString(useReco?"":"gen_")+"explore_dy_"+varname
                                +"_DY-FEWZ_with_selection_"+Form("%d",(int)luminosity)+".root", "RECREATE");
    //TDirectory* stacks = savefile->mkdir("stacks");
    TDirectory* histos = savefile->mkdir("histos");

    // save the different histos and stacks in the appropriate directories in the tfile
    //stacks->cd();
    //varstacklist->Write();
    histos->cd();
    histolist->Write();
    savefile->Close();
    return 0;
}
