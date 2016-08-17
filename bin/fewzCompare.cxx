// Missing HLT trigger info in CMSSW_8_0_X MC so we have to compare Data and MC in a different manner.
// We apply triggers to data but not to MC. Then scale MC for trigger efficiency.


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

  
    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    float luminosity = 3990;       // pb-1
    //float luminosity = 12900;      // pb-1 (ICHEP)
    //float luminosity = 15900;      // pb-1   (2016-08-03)
    float triggerSF = 0.913;       // no HLT trigger info available for the samples so we scale for the trigger efficiency instead
    float signalSF = 100;
/*
    // ================================================================
    // Data -----------------------------------------------------------
    // ================================================================
 

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

    TString dyfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/bg/dy/CMSSW_8_0_X/stage_1_dy_jetsToLL_ALL.root");
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
    FEWZCompareCuts FEWZselection;
    CategorySelectionFEWZ categorySelectionFEWZ;

    TString varname;
    int nbins, wbins;
    float nmin, wmin;
    float nmax, wmax;

    // recoCandMass
    if(input == 0)
    {
        nbins = 50;
        wbins = 100;

        nmin = 110;
        nmax = 160;

        wmin = 110;
        wmax = 310;

        varname = "dimuMass";
    }

    // recoPt
    if(input == 1)
    {
        nbins = 200;
        wbins = 200;

        nmin = 0;
        nmax = 200;

        wmin = 0;
        wmax = 200;
        varname = "recoMu_Pt";
    }
 
    // recoEta
    if(input == 2)
    {
        nbins = 50;
        wbins = 50;

        nmin = -3;
        nmax = 3;

        wmin = -3;
        wmax = 3;
        varname = "recoMu_Eta";
    }

    // jetPt
    if(input == 3)
    {
        nbins = 200;
        wbins = 200;

        nmin = 0;
        nmax = 200;

        wmin = 0;
        wmax = 200;
        varname = "jet_Pt";
    }

    // jet_Eta 
    if(input == 4)
    {
        nbins = 100;
        wbins = 50;

        nmin = -5;
        nmax = 5;

        wmin = -5;
        wmax = 5;
        varname = "jet_Eta";
    }

    std::cout << std::endl;
    std::cout << "======== Plot Configs ========" << std::endl;
    std::cout << "var         : " << varname << std::endl;
    std::cout << "nmin         : " << nmin << std::endl;
    std::cout << "nmax         : " << nmax << std::endl;
    std::cout << "nbins        : " << nbins << std::endl;
    std::cout << std::endl;
    std::cout << "wmin         : " << wmin << std::endl;
    std::cout << "wmax         : " << wmax << std::endl;
    std::cout << "wbins        : " << wbins << std::endl;
    std::cout << std::endl;

    // The lists of histograms used to make the stacks
    TList* varlistalln = new TList();   // all events passing cuts, narrow mass window 
    TList* varlistallw = new TList();   // all events passing cuts, wide mass window
    TList* varlistccn = new TList();    // central central narrow
    TList* varlistccw = new TList();    // central central wide
    TList* varlistcncn = new TList();   // central not central narrow
    TList* varlistcncw = new TList();   // central not central wide
    TList* varlist1jetn = new TList();  // 1 jet narrow
    TList* varlist1jetw = new TList();  // 1 jet wide

    // Not sure how to deal with the scaling correctly when using a subset of events
    float reductionFactor = 1;

    for(auto const &s : samplevec)
    {
      // Output some info about the current file
      std::cout << std::endl;
      std::cout << "  /// Looping over " << s->name << std::endl;

      ///////////////////////////////////////////////////////////////////
      // HISTOGRAMS TO FILL ---------------------------------------------
      ///////////////////////////////////////////////////////////////////

      TH1F* varhistoalln  = new TH1F(varname+"_all_n_"+s->name, varname+"_all_n_"+s->name, nbins, nmin, nmax);
      TH1F* varhistoallw  = new TH1F(varname+"_all_w_"+s->name, varname+"_all_w_"+s->name, wbins, wmin, wmax);

      // Different categories for the analysis
      TH1F* varhistoccn = new TH1F(varname+"_ccn_"+s->name, varname+"_ccn_"+s->name, nbins, nmin, nmax);
      TH1F* varhistoccw = new TH1F(varname+"_ccw_"+s->name, varname+"_ccw_"+s->name, wbins, wmin, wmax);
      TH1F* varhistocncn = new TH1F(varname+"_cncn_"+s->name, varname+"_cncn_"+s->name, nbins, nmin, nmax);
      TH1F* varhistocncw  = new TH1F(varname+"_cncw_"+s->name, varname+"_cncw_"+s->name, wbins, wmin, wmax);
      TH1F* varhisto1jetn  = new TH1F(varname+"_1jetn_"+s->name, varname+"_1jetn_"+s->name, nbins, nmin, nmax);
      TH1F* varhisto1jetw  = new TH1F(varname+"_1jetw_"+s->name, varname+"_1jetw_"+s->name, wbins, wmin, wmax);

      for(unsigned int i=0; i<s->N/reductionFactor; i++)
      {

        ///////////////////////////////////////////////////////////////////
        // GET INFORMATION ------------------------------------------------
        ///////////////////////////////////////////////////////////////////

        s->getEntry(i); 
        s->vars.validJets = std::vector<TLorentzVector>();
        jetSelectionTools.getValidJets(s->vars, s->vars.validJets);
        std::pair<int,int> e(s->vars.eventInfo.run, s->vars.eventInfo.event); // create a pair that identifies the event uniquely

        ///////////////////////////////////////////////////////////////////
        // CUTS  ----------------------------------------------------------
        ///////////////////////////////////////////////////////////////////

        if(!FEWZselection.evaluate(s->vars)) 
        {
            continue; 
        }

        // Figure out which category the event belongs to
        categorySelectionFEWZ.evaluate(s->vars);

        // recoCandMass
        if(varname.Contains("dimuMass")) 
        {
            float varvalue = -9999;
            varvalue = s->vars.recoCandMass;
            if(!(s->sampleType.Contains("data") && varvalue >= 110 && varvalue < 140))
            {
                if(categorySelectionFEWZ.isWide) 
                    varhistoallw->Fill(varvalue, s->getWeight());

                if(categorySelectionFEWZ.isNarrow) 
                    varhistoalln->Fill(varvalue, s->getWeight());

                if(categorySelectionFEWZ.isCentralCentralWide) 
                    varhistoccw->Fill(varvalue, s->getWeight());

                if(categorySelectionFEWZ.isCentralCentralNarrow)  
                    varhistoccn->Fill(varvalue, s->getWeight());

                if(categorySelectionFEWZ.isCentralNotCentralWide)  
                    varhistocncw->Fill(varvalue, s->getWeight());

                if(categorySelectionFEWZ.isCentralNotCentralNarrow)  
                    varhistocncn->Fill(varvalue, s->getWeight());

                if(categorySelectionFEWZ.isOneJetInclusiveWide)  
                    varhisto1jetw->Fill(varvalue, s->getWeight());

                if(categorySelectionFEWZ.isOneJetInclusiveNarrow)  
                    varhisto1jetn->Fill(varvalue, s->getWeight());
            }
        }

        // recoMu_Pt
        if(varname.Contains("recoMu_Pt"))
        {
            float varvalue = -9999;

            // reco1
            varvalue = s->vars.reco1.pt;
            if(categorySelectionFEWZ.isWide) 
                varhistoallw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isNarrow) 
                varhistoalln->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralCentralWide) 
                varhistoccw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralCentralNarrow)  
                varhistoccn->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralNotCentralWide)  
                varhistocncw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralNotCentralNarrow)  
                varhistocncn->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isOneJetInclusiveWide)  
                varhisto1jetw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isOneJetInclusiveNarrow)  
                varhisto1jetn->Fill(varvalue, s->getWeight());

            // reco2  
            varvalue = s->vars.reco2.pt;
            if(categorySelectionFEWZ.isWide) 
                varhistoallw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isNarrow) 
                varhistoalln->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralCentralWide) 
                varhistoccw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralCentralNarrow)  
                varhistoccn->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralNotCentralWide)  
                varhistocncw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralNotCentralNarrow)  
                varhistocncn->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isOneJetInclusiveWide)  
                varhisto1jetw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isOneJetInclusiveNarrow)  
                varhisto1jetn->Fill(varvalue, s->getWeight());
        }

        // recoMu_Eta
        if(varname.Contains("recoMu_Eta"))
        {
            float varvalue = -9999;

            // reco1
            varvalue = s->vars.reco1.eta;
            if(categorySelectionFEWZ.isWide) 
                varhistoallw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isNarrow) 
                varhistoalln->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralCentralWide) 
                varhistoccw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralCentralNarrow)  
                varhistoccn->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralNotCentralWide)  
                varhistocncw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralNotCentralNarrow)  
                varhistocncn->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isOneJetInclusiveWide)  
                varhisto1jetw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isOneJetInclusiveNarrow)  
                varhisto1jetn->Fill(varvalue, s->getWeight());

            // reco2
            varvalue = s->vars.reco2.eta;
            if(categorySelectionFEWZ.isWide) 
                varhistoallw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isNarrow) 
                varhistoalln->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralCentralWide) 
                varhistoccw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralCentralNarrow)  
                varhistoccn->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralNotCentralWide)  
                varhistocncw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isCentralNotCentralNarrow)  
                varhistocncn->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isOneJetInclusiveWide)  
                varhisto1jetw->Fill(varvalue, s->getWeight());

            if(categorySelectionFEWZ.isOneJetInclusiveNarrow)  
                varhisto1jetn->Fill(varvalue, s->getWeight());
        }

        // jet_Pt
        if(varname.Contains("jet_Pt"))
        {
            float varvalue = -9999;

            if(categorySelectionFEWZ.isWide) 
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Pt();
                    varhistoallw->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelectionFEWZ.isNarrow) 
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Pt();
                    varhistoalln->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelectionFEWZ.isCentralCentralWide) 
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Pt();
                    varhistoccw->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelectionFEWZ.isCentralCentralNarrow)  
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Pt();
                    varhistoccn->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelectionFEWZ.isCentralNotCentralWide)  
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Pt();
                    varhistocncw->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelectionFEWZ.isCentralNotCentralNarrow)  
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Pt();
                    varhistocncn->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelectionFEWZ.isOneJetInclusiveWide)  
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Pt();
                    varhisto1jetw->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelectionFEWZ.isOneJetInclusiveNarrow)  
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Pt();
                    varhisto1jetn->Fill(varvalue, s->getWeight());
                }
            }
        }

        // jet_Eta
        if(varname.Contains("jet_Eta"))
        {
            float varvalue = -9999;

            if(categorySelectionFEWZ.isWide) 
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Eta();
                    varhistoallw->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelectionFEWZ.isNarrow) 
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Eta();
                    varhistoalln->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelectionFEWZ.isCentralCentralWide) 
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Eta();
                    varhistoccw->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelectionFEWZ.isCentralCentralNarrow)  
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Eta();
                    varhistoccn->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelectionFEWZ.isCentralNotCentralWide)  
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Eta();
                    varhistocncw->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelectionFEWZ.isCentralNotCentralNarrow)  
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Eta();
                    varhistocncn->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelectionFEWZ.isOneJetInclusiveWide)  
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Eta();
                    varhisto1jetw->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelectionFEWZ.isOneJetInclusiveNarrow)  
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Eta();
                    varhisto1jetn->Fill(varvalue, s->getWeight());
                }
            }
        }

        // Reset the flags in preparation for the next event
        categorySelectionFEWZ.reset();
      }

      // Scale according to luminosity
      varhistoallw->Scale(s->getScaleFactor(luminosity));
      varhistoalln->Scale(s->getScaleFactor(luminosity));
      varhistoccw->Scale(s->getScaleFactor(luminosity));
      varhistoccn->Scale(s->getScaleFactor(luminosity));
      varhistocncw->Scale(s->getScaleFactor(luminosity));
      varhistocncn->Scale(s->getScaleFactor(luminosity));
      varhisto1jetw->Scale(s->getScaleFactor(luminosity));
      varhisto1jetn->Scale(s->getScaleFactor(luminosity));

      // No trigger info in 80X MC samples, scale for trigger efficiency
      if(!s->sampleType.Contains("data"))
      {
          //varhistoallw->Scale(s->getScaleFactor(triggerSF));
          //varhistoalln->Scale(s->getScaleFactor(triggerSF));
          //varhistoccw->Scale(s->getScaleFactor(triggerSF);
          //varhistoccn->Scale(s->getScaleFactor(triggerSF);
          //varhistocncw->Scale(s->getScaleFactor(triggerSF));
          //varhistocncn->Scale(s->getScaleFactor(triggerSF));
          //varhisto1jetw->Scale(s->getScaleFactor(triggerSF));
          //varhisto1jetn->Scale(s->getScaleFactor(triggerSF));
      }

      // Add to the appropriate list
      varlistallw->Add(varhistoallw);
      varlistalln->Add(varhistoalln);
      varlistccw->Add(varhistoccw);
      varlistccn->Add(varhistoccn);
      varlist1jetw->Add(varhisto1jetw);
      varlist1jetn->Add(varhisto1jetn);
    }

    // ////////////////////////////////////////////////////////////////////////////
    // ========= Scale, Stack, Save ===============================================
    // ////////////////////////////////////////////////////////////////////////////

    //TIter next(varlist);
    //TObject* object = 0;
    //while( (object = next()) )
    //{
    //  TH1F* varhisto = (TH1F*) object;
    //  if(TString(varhisto->GetName()).Contains("signal"))
    //  {
    //      // scale the signal so that it's easier to see on the plots
    //      // only do this right before saving or it would skew the significance results
    //      varhisto->Scale(signalSF);
    //  }
    //}

    // Create the stack and ratio plot    
    //TCanvas* varstackcanvasall = dps->stackedHistogramsAndRatio(varlistall, "c_all_"+varname, varname+"_all_stack", varname, "Num Entries");
    //TCanvas* varstackcanvasVBFt = dps->stackedHistogramsAndRatio(varlistVBFt, "c_vbft_"+varname, varname+"_vbft_stack", varname, "Num Entries");
    //TCanvas* varstackcanvasVBFl = dps->stackedHistogramsAndRatio(varlistVBFl, "c_vbfl_"+varname, varname+"_vbfl_stack", varname, "Num Entries");
    //TCanvas* varstackcanvasGGFt = dps->stackedHistogramsAndRatio(varlistGGFt, "c_ggft_"+varname, varname+"_ggft_stack", varname, "Num Entries");
    //TCanvas* varstackcanvas01t = dps->stackedHistogramsAndRatio(varlist01t, "c_01t_"+varname, varname+"_01t_stack", varname, "Num Entries");
    //TCanvas* varstackcanvas01l = dps->stackedHistogramsAndRatio(varlist01l, "c_01l_"+varname, varname+"_01l_stack", varname, "Num Entries");
    //std::cout << std::endl;

    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile("rootfiles/dy_hists_for_fewz_compare_"+varname+".root", "RECREATE");
    //TDirectory* stacks = savefile->mkdir("stacks");
    TDirectory* histos = savefile->mkdir("histos");

    // save the different histos in the appropriate directories in the tfile
    //stacks->cd();
    //varstackcanvasall->Write();
    //varstackcanvasVBFt->Write();
    //varstackcanvasGGFt->Write();
    //varstackcanvasVBFl->Write();
    //varstackcanvas01t->Write();
    //varstackcanvas01l->Write();

    histos->cd();
    varlistallw->Write();
    varlistalln->Write();
    varlistccw->Write();
    varlistccn->Write();
    varlist1jetw->Write();
    varlist1jetn->Write();

    savefile->Close();

    return 0;
}
