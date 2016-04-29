#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "CutSet.h"
#include "Cut.h"
#include "SelectionCuts.h"
#include "CategorySelection.h"
#include "JetSelectionTools.h"

#include "EventTools.cxx"
#include "PUTools.cxx"
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
    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) ss >> input;
    }   

    // Not sure that we need a map if we have a vector
    // Should use this as the main database and choose from it to make the vector
    std::map<std::string, Sample*> samples;

    // Second container so that we can have a copy sorted by cross section.
    std::vector<Sample*> samplevec;

    // Use this to plot some things if we wish
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

    //CutSet cuts();
    //Cut cut();  // Abstract class, can't initialize
  
    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    float luminosity = 2169;
    float signalSF = 100;

    // ================================================================
    // Data -----------------------------------------------------------
    // ================================================================
 

    TString datafilename = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data_from_json/25ns/golden/stage_1_singleMuon_RunCD_GOLDEN_ALL.root");
    //TString datafilename = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data_from_json/25ns/golden/stage_1_singleMuon_RunDBoth_MINIAOD_GOLDEN_ALL.root");
    Sample* datasample = new Sample(datafilename, "Data", "data");
    datasample->lumi = luminosity;
    datasample->xsec = 9999;
    //datasample->pileupfile = "/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data_from_json/25ns/golden/pileup/old/PUCalib_Golden_71mb.root";
    datasample->pileupfile = "/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data_from_json/25ns/golden/pileup/pileup_data_71mb.root";
    samples["Data"] = datasample;

    // ================================================================
    // DYJetsToLL -----------------------------------------------------
    // ================================================================

    TString dyfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/monte_carlo/bg/dy/stage_1_dy_jetsToLL_ALL.root");
    samples["DYJetsToLL"] = new Sample(dyfilename, "DYJetsToLL", "background");
    samples["DYJetsToLL"]->pileupfile = "./pu_reweight_trees/PUCalib_gw_DYJetsToLL.root"; //nPU
    samples["DYJetsToLL"]->xsec = 6025.2; // pb

    // ================================================================
    // TTJets ---------------------------------------------------------
    // ================================================================

    TString ttbarfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/monte_carlo/bg/ttbar/stage_1_ttJets_ALL.root");
    samples["TTJets"] = new Sample(ttbarfilename, "TTJets", "background");
    samples["TTJets"]->pileupfile = "./pu_reweight_trees/PUCalib_gw_TTJets.root"; //nPU
    samples["TTJets"]->xsec = 831.76; // pb

    // ================================================================
    // VBF ---------------------------------------------------------
    // ================================================================

    TString vbffilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/monte_carlo/signal/stage_1_vbf_HToMuMu_ALL.root");
    samples["VBF"] = new Sample(vbffilename, "VBF", "signal");
    samples["VBF"]->pileupfile = "./pu_reweight_trees/PUCalib_gw_VBF.root"; //nPU
    samples["VBF"]->xsec = 3.727*0.00022; // pb

    // ================================================================
    // GGF ---------------------------------------------------------
    // ================================================================

    TString ggfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/monte_carlo/signal/stage_1_gg_HToMuMu_ALL.root");
    samples["GGF"] = new Sample(ggfilename, "GGF", "signal");
    samples["GGF"]->pileupfile = "./pu_reweight_trees/PUCalib_gw_GGF.root"; //nPU
    samples["GGF"]->xsec = 43.62*0.00022; // pb

    ///////////////////////////////////////////////////////////////////
    // PREPROCESSING---------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // Loop through all of the samples to do some pre-processing
    std::cout << std::endl;
    std::cout << "======== Preprocess the samples... " << std::endl;
    std::cout << std::endl;

    //makePUHistos(samples);
    
    for(auto const &i : samples)
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

        if(!i.second->sampleType.Contains("data"))
        {
            // Pileup reweighting
            std::cout << "    +++ PU Reweighting " << i.second->name << "..."  << std::endl;
            std::cout << std::endl;

            i.second->lumiWeights = new reweight::LumiReWeighting(i.second->pileupfile.Data(), samples["Data"]->pileupfile.Data(), "pileup", "pileup");
            std::cout << "        " << i.first << "->lumiWeights: " << i.second->lumiWeights << std::endl;
            std::cout << std::endl;
        }
        samplevec.push_back(i.second);
    }

    // Sort the samples by xsec. Useful when making the histogram stack.
    std::sort(samplevec.begin(), samplevec.end(), [](Sample* a, Sample* b){ return a->xsec < b->xsec; }); 
    

    ///////////////////////////////////////////////////////////////////
    // Cut and Categorize ---------------------------------------------
    ///////////////////////////////////////////////////////////////////
    
    // Objects to help with the cuts and selections
    TightMuonIdCuts tightMuonId;
    JetSelectionTools jetSelectionTools;
    CategorySelection categorySelection;
    Run1EventSelectionCuts run1EventSelection;
    Run1MuonSelectionCuts run1MuonSelection;

    // vars that tell us the settings for the histograms
    bool muonselection = false;
    bool eventselection = false;
    int bins, min, max;
    TString varname;

    // event selection
    // 0 = charge1 != charge2
    // 1 = cTrigMuPtMin
    // 2 = cDimuMassMin
    if(input >=1 && input <=2)
    {
        run1EventSelection.cutset.cuts[input].on = false;
        eventselection = true;
        bins = run1EventSelection.cutset.cuts[input].bins;
        min = run1EventSelection.cutset.cuts[input].min;
        max = run1EventSelection.cutset.cuts[input].max;
        varname = run1EventSelection.cutset.cuts[input].name;
        std::cout << std::endl;
        std::cout << "======== N-1 Settings ========" << std::endl;
        std::cout << "Selection Type: Event Selection" << std::endl;
        std::cout << "cut #         : " << input << std::endl;
        std::cout << "cut name      : " << varname << std::endl;
        std::cout << std::endl;
    }
    // muon selection
    // 0 = cMinPt
    // 1 = cMaxEta
    // 2 = cMaxRelIso
    // 3 = cMinPt
    // 4 = cMaxEta
    // 5 = cMaxRelIso
    // first three for reco1, second three for reco2, though both use the same value for the cut
    else if(input >= 3 && input <= 8)
    {
        // adjust the input by three so that input=3 is the 0th muon selection cut
        input-=3;
        run1MuonSelection.cutset.cuts[input].on = false;
        muonselection = true;
        bins = run1MuonSelection.cutset.cuts[input].bins;
        min = run1MuonSelection.cutset.cuts[input].min;
        max = run1MuonSelection.cutset.cuts[input].max;
        varname = run1MuonSelection.cutset.cuts[input].name;
        std::cout << std::endl;
        std::cout << "======== N-1 Settings ========" << std::endl;
        std::cout << "Selection Type: Muon Selection" << std::endl;
        std::cout << "cut #         : " << input << std::endl;
        std::cout << "cut name      : " << varname << std::endl;
        std::cout << std::endl;
    }
    else
    {
        input=-1;
        bins = 1;
        min = 0;
        max = 1;
        varname = "nm0";
        std::cout << "======== N-1 Settings ========" << std::endl;
        std::cout << "Selection Type: N-0" << std::endl;
        std::cout << "cut #         : " << input << std::endl;
        std::cout << "cut name      : " << varname << std::endl;
        std::cout << std::endl;
    }
    // remove dots and replace with underscores 
    varname.ReplaceAll(".","_");

    TList* varlist = new TList();   // the list of histograms used to make the stack for var
    TList* nlist = new TList();     // the list of histograms used to make the stack for the counts
    double nsignal = 0;             // the number of signal events
    double nbackground = 0;         // the number of background events
    double ndata = 0;               // the number of background events

    // Not sure how to deal with the scaling correctly when using a subset of events
    float reductionFactor = 100;

    // make histograms to count everything
    TH1F* nhistosignal = new TH1F("counts_"+TString("signal"), "counts_"+TString("signal"), 1, 0, 1);
    TH1F* nhistobg = new TH1F("counts_"+TString("background"), "counts_"+TString("background"), 1, 0, 1);
    TH1F* nhistodata = new TH1F("counts_"+TString("data"), "counts_"+TString("data"), 1, 0, 1);

    // make histograms for the variable of interest
    TH1F* varhistosignal = new TH1F(varname+"_"+TString("signal"), varname+"_"+TString("signal"), bins, min, max);
    TH1F* varhistobg = new TH1F(varname+"_"+TString("background"), varname+"_"+TString("background"), bins, min, max);
    TH1F* varhistodata = new TH1F(varname+"_"+TString("data"), varname+"_"+TString("data"), bins, min, max);

    for(auto const &s : samplevec)
    {
      // Output some info about the current file
      std::cout << std::endl;
      std::cout << "  /// Looping over " << s->name << std::endl;

      ///////////////////////////////////////////////////////////////////
      // HISTOGRAMS TO FILL ---------------------------------------------
      ///////////////////////////////////////////////////////////////////

      TH1F* varhisto = new TH1F(varname+"_"+s->name, varname+"_"+s->name, bins, min, max);
      TH1F* nhisto = new TH1F("counts_"+s->name, "counts_"+s->name, 1, 0, 1);

      for(unsigned int i=0; i<s->N/reductionFactor; i++)
      {

        ///////////////////////////////////////////////////////////////////
        // GET INFORMATION ------------------------------------------------
        ///////////////////////////////////////////////////////////////////

        s->getEntry(i); 
        s->vars.validJets = std::vector<TLorentzVector>();
        jetSelectionTools.getValidJetsdR(s->vars, s->vars.validJets);
        std::pair<int,int> e(s->vars.eventInfo.run, s->vars.eventInfo.event); // create a pair that identifies the event uniquely

        ///////////////////////////////////////////////////////////////////
        // CUTS  ----------------------------------------------------------
        ///////////////////////////////////////////////////////////////////

        if(!tightMuonId.evaluate(s->vars))
        { 
            continue; 
        }
        if(!run1EventSelection.evaluate(s->vars))
        { 
            continue; 
        }
        if(!run1MuonSelection.evaluate(s->vars)) 
        {
            continue; 
        }
   
        float varvalue = 0.5; 
        if(input >= 0)
        {
            if(eventselection) varvalue = run1EventSelection.cutset.cuts[input].value;
            if(muonselection) varvalue = run1MuonSelection.cutset.cuts[input].value;
        }

        varhisto->Fill(varvalue, s->getWeight());   
        nhisto->Fill(0.5, s->getWeight());   
        //varhisto->Fill(s->vars.reco2.var, s->getWeight());   
        
        // The event passes the cuts, see which category it falls into
        //categorySelection.evaluate(s->vars);
        //if(categorySelection.isVBFTight) 
        //{
        //}
        //else if(categorySelection.isGGFTight)
        //{
        //}
        //else if(categorySelection.isVBFLoose)
        //{
        //}
        //else if (categorySelection.isTight01)
        //{
        //}
        //else if (categorySelection.isLoose01)
        //{
        //}

        //if(eventInVector(e, eventsToCheck) || i < numToPrint) // Adrian gave a list of events to look at for synch purposes
        if(false)
          // ouput pt, mass info etc
          outputEvent(s->vars, categorySelection);

        // Reset the flags in preparation for the next event
        categorySelection.reset();
      }
      varhisto->Scale(s->getScaleFactor(luminosity));
      nhisto->Scale(s->getScaleFactor(luminosity));
      if(s->sampleType.Contains("signal"))
      {
          // scale the signal so that it's easier to see on the plots
          varhisto->Scale(signalSF);
          varhistosignal->Add(varhisto); 
          nhistosignal->Add(nhisto); 
      }
      if(s->sampleType.Contains("background"))
      {
          varhistobg->Add(varhisto); 
          nhistobg->Add(nhisto); 
      }
      if(s->sampleType.Contains("data")) 
      {
          varhistodata->Add(varhisto); 
          nhistodata->Add(nhisto); 
      }
      varlist->Add(varhisto);
      nlist->Add(nhisto);
    }

    // Create the stack and ratio plot    
    TCanvas* varstackcanvas;
    varstackcanvas = dps->stackedHistogramsAndRatio(varlist, varname+"_stack", varname, "Num Entries");
    varstackcanvas->SetName("c_"+varname);
    std::cout << std::endl;

    TCanvas* nstackcanvas;
    nstackcanvas = dps->stackedHistogramsAndRatio(nlist, "nstack", "counted", "Num Entries");
    nstackcanvas->SetName("c_counts");
    std::cout << std::endl;

    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile("rootfiles/nm1_"+varname+".root", "RECREATE");

    // save the different histos in the appropriate directories in the tfile
    TDirectory* stacks = savefile->mkdir("stacks");
    stacks->cd();
    varstackcanvas->Write();
    nstackcanvas->Write();

    TDirectory* datahistos = savefile->mkdir("data");
    datahistos->cd();
    varhistodata->Write();
    nhistodata->Write();

    TDirectory* signalhistos = savefile->mkdir("signal");
    signalhistos->cd();
    varhistosignal->Write();
    nhistosignal->Write();

    TDirectory* bghistos = savefile->mkdir("bg");
    bghistos->cd();
    varhistobg->Write();
    nhistobg->Write();

    savefile->Close();


    // ////////////////////////////////////////////////////////////////////////////
    // ========= Total Counts =====================================================
    // ////////////////////////////////////////////////////////////////////////////
    
    nsignal = nhistosignal->Integral();
    nbackground = nhistobg->Integral();
    ndata = nhistodata->Integral();

    double nsignalv = varhistosignal->Integral()/signalSF;
    double nbackgroundv = varhistobg->Integral();
    double ndatav = varhistodata->Integral();

    std::cout << std::endl;
    std::cout << "=========== Total Counts ============" << std::endl;
    std::cout << "Signal:     " << nsignal << std::endl;
    std::cout << "Background: " << nbackground << std::endl;
    std::cout << "Total:      " << nsignal + nbackground << std::endl;
    std::cout << std::endl;
    std::cout << "vSignal:     " << nsignalv << std::endl;
    std::cout << "vBackground: " << nbackgroundv << std::endl;
    std::cout << "vTotal:      " << nsignalv + nbackgroundv << std::endl;
    std::cout << std::endl;

    // ////////////////////////////////////////////////////////////////////////////
    // ========= Significance =====================================================
    // ////////////////////////////////////////////////////////////////////////////
    
    AsimovSignificance asimov;
    PoissonSignificance poisson;

    // unc = 0.06 for 12596 background events, scales like 1/sqrt(N)
    double asimovZ0 = asimov.significance(nsignal, nbackground, 0);
    double asimovZ1 = asimov.significance(nsignal, nbackground, 0.051);
    double rootbZ0 = poisson.significance(nsignal, nbackground, 0);
    double rootbZ1 = poisson.significance(nsignal, nbackground, 0.051);

    std::cout << std::endl;
    std::cout << "=========== Significance ============" << std::endl;
    std::cout << "Asimov w/o bgfit:    " << asimovZ0 << std::endl;
    std::cout << "Asimov w/ bgfit:     " << asimovZ1 << std::endl;
    std::cout << "s/sqrt(b):           " << rootbZ0 << std::endl;
    std::cout << "s/sqrt(b + varb):    " << rootbZ1 << std::endl;
    std::cout << std::endl;


    // ////////////////////////////////////////////////////////////////////////////
    // ========= Total Counts =====================================================
    // ////////////////////////////////////////////////////////////////////////////
    
    //std::cout << std::endl;
    //std::cout << "=========== Category Counts ============" << std::endl;
    //std::cout << "VBFTight: " << VBFTightVec.size() << std::endl;
    //std::cout << "GGFTight: " << GGFTightVec.size() << std::endl;
    //std::cout << "VBFLoose: " << VBFLooseVec.size() << std::endl;
    //std::cout << "10Tight: " << tight01Vec.size() <<  std::endl;
    //std::cout << "10Loose: " << loose01Vec.size() <<  std::endl;
    //std::cout << std::endl;

    //int dimutotal = VBFTightVec.size() + GGFTightVec.size() + VBFLooseVec.size() + tight01Vec.size() + loose01Vec.size();
    //std::cout << "total: " << dimutotal << std::endl;
    //std::cout << std::endl;


    return 0;
}
