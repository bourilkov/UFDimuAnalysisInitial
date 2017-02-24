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

    // ================================================================
    // DYJetsToLL -----------------------------------------------------
    // ================================================================

    TString dyfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/bg/dy/CMSSW_8_0_X/stage_1_dy_jetsToLL_ALL.root");
    samples["DYJetsToLL"] = new Sample(dyfilename, "DYJetsToLL", "background");
    samples["DYJetsToLL"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_DYJetsToLL.root"; //nPU
    samples["DYJetsToLL"]->xsec = 6025.2; // pb

    // ================================================================
    // TTJets ---------------------------------------------------------
    // ================================================================

    TString ttbarfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/bg/ttbar/CMSSW_8_0_X/stage_1_ttJets2_ALL.root");
    samples["TTJets"] = new Sample(ttbarfilename, "TTJets", "background");
    samples["TTJets"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_TTJets.root"; //nPU
    samples["TTJets"]->xsec = 831.76; // pb

    // ================================================================
    // VBF ---------------------------------------------------------
    // ================================================================

    TString vbffilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/signal/CMSSW_8_0_X/stage_1_vbf_HToMuMu_ALL.root");
    samples["VBF"] = new Sample(vbffilename, "VBF", "signal");
    samples["VBF"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_VBF.root"; //nPU
    samples["VBF"]->xsec = 3.727*0.00022; // pb

    // ================================================================
    // GGF ---------------------------------------------------------
    // ================================================================

    TString ggfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/signal/CMSSW_8_0_X/stage_1_gg_HToMuMu2_ALL.root");
    samples["GGF"] = new Sample(ggfilename, "GGF", "signal");
    samples["GGF"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_GGF.root"; //nPU
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
    JetSelectionTools jetSelectionTools;
    CategorySelection categorySelection;
    Run1MuonSelectionCuts run1MuonSelection;
    Run1EventSelectionCuts80X run1EventSelectionData(true);
    Run1EventSelectionCuts80X run1EventSelectionMC;

    TString varname;
    int bins;
    float min;
    float max;

    // recoCandMass
    if(input == 0)
    {
        bins = 150;
        min = 50;
        max = 200;
        varname = "dimuMass";
    }

    // dimuPt 
    if(input == 1)
    {
        bins = 200;
        min = 0;
        max = 100;
        varname = "dimuPt";
    }

    // recoPt
    if(input == 2)
    {
        bins = 200;
        min = 0;
        max = 150;
        varname = "recoMu_Pt";
    }
 
    // recoEta
    if(input == 3)
    {
        bins = 100;
        min = -2.5;
        max = 2.5;
        varname = "recoMu_Eta";
    }

    // NPV
    if(input == 4)
    {
        bins = 50;
        min = 0;
        max = 50;
        varname = "NPV";
    }

    // jet_Pt
    if(input == 5)
    {   
        bins = 200;
        min = 0;
        max = 200;
        varname = "jet_Pt";
    }   

    // jet_Eta 
    if(input == 6)
    {   
        bins = 100;
        min = -5; 
        max = 5;
        varname = "jet_Eta";
    }   

    // N_valid_jets
    if(input == 7)
    {   
        bins = 11;
        min = 0; 
        max = 11;
        varname = "N_valid_jets";
    }   

    // Dijet_mass
    if(input == 8)
    {   
        bins = 200;
        min = 0; 
        max = 2000;
        varname = "Dijet_mass";
    }   

    // dEta_jj
    if(input == 9)
    {   
        bins = 100;
        min = -10; 
        max = 10;
        varname = "dEta_jj";
    }   



    std::cout << std::endl;
    std::cout << "======== Plot Configs ========" << std::endl;
    std::cout << "var         : " << varname << std::endl;
    std::cout << "min         : " << min << std::endl;
    std::cout << "max         : " << max << std::endl;
    std::cout << "bins        : " << bins << std::endl;
    std::cout << std::endl;

    // The lists of histograms used to make the stacks
    TList* varlistall = new TList();   // all events passing cuts 
    TList* varlistVBFt = new TList();  // VBF Tight category
    TList* varlistVBFl = new TList();  // VBF Loose
    TList* varlistGGFt = new TList();  // GGF Tight
    TList* varlist01t = new TList();   // 01 Jet Tight
    TList* varlist01l = new TList();   // 01 Jet Loose

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

      TH1F* varhistoall  = new TH1F(varname+"_all_"+s->name, varname+"_all_"+s->name, bins, min, max);

      int lowstatsbins = bins/5;
      if(varname.Contains("N")) lowstatsbins = bins;

      // Different categories for the analysis
      TH1F* varhistoVBFt = new TH1F(varname+"_vbft_"+s->name, varname+"_vbft_"+s->name, lowstatsbins, min, max);
      TH1F* varhistoVBFl = new TH1F(varname+"_vbfl_"+s->name, varname+"_vbfl_"+s->name, lowstatsbins, min, max);
      TH1F* varhistoGGFt = new TH1F(varname+"_ggft_"+s->name, varname+"_ggft_"+s->name, lowstatsbins, min, max);
      TH1F* varhisto01t  = new TH1F(varname+"_01t_"+s->name, varname+"_01t_"+s->name, bins, min, max);
      TH1F* varhisto01l  = new TH1F(varname+"_01l_"+s->name, varname+"_01l_"+s->name, bins, min, max);

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

        if(!s->vars.reco1.isTightMuon || !s->vars.reco2.isTightMuon)
        { 
            continue; 
        }
        if(!run1EventSelectionData.evaluate(s->vars) && s->sampleType.Contains("data"))
        { 
            continue; 
        }
        if(!run1EventSelectionMC.evaluate(s->vars) && !s->sampleType.Contains("data"))
        { 
            continue; 
        }
        if(!run1MuonSelection.evaluate(s->vars)) 
        {
            continue; 
        }

        // Figure out which category the event belongs to
        categorySelection.evaluate(s->vars);

        // recoCandMass
        if(varname.Contains("dimuMass")) 
        {
            float varvalue = -9999;
            varvalue = s->vars.recoCandMassPF;
            if(!(s->sampleType.Contains("data") && varvalue >= 110 && varvalue < 140))
            {
                varhistoall->Fill(varvalue, s->getWeight());   
                if(categorySelection.isVBFTight) varhistoVBFt->Fill(varvalue, s->getWeight());
                if(categorySelection.isGGFTight) varhistoGGFt->Fill(varvalue, s->getWeight());
                if(categorySelection.isVBFLoose) varhistoVBFl->Fill(varvalue, s->getWeight());
                if(categorySelection.isTight01)  varhisto01t->Fill(varvalue, s->getWeight());
                if(categorySelection.isLoose01)  varhisto01l->Fill(varvalue, s->getWeight());
            }
        }

        // recoCandPt
        if(varname.Contains("dimuPt"))
        {
             varhistoall->Fill(s->vars.recoCandPtPF, s->getWeight());
             if(categorySelection.isVBFTight) varhistoVBFt->Fill(s->vars.recoCandPtPF, s->getWeight());
             if(categorySelection.isGGFTight) varhistoGGFt->Fill(s->vars.recoCandPtPF, s->getWeight());
             if(categorySelection.isVBFLoose) varhistoVBFl->Fill(s->vars.recoCandPtPF, s->getWeight());
             if(categorySelection.isTight01)  varhisto01t->Fill(s->vars.recoCandPtPF, s->getWeight());
             if(categorySelection.isLoose01)  varhisto01l->Fill(s->vars.recoCandPtPF, s->getWeight());
        }

        // recoMu_Pt
        if(varname.Contains("recoMu_Pt"))
        {
             varhistoall->Fill(s->vars.reco1.pt, s->getWeight());
             varhistoall->Fill(s->vars.reco2.pt, s->getWeight());

             if(categorySelection.isVBFTight) varhistoVBFt->Fill(s->vars.reco1.pt, s->getWeight());
             if(categorySelection.isGGFTight) varhistoGGFt->Fill(s->vars.reco1.pt, s->getWeight());
             if(categorySelection.isVBFLoose) varhistoVBFl->Fill(s->vars.reco1.pt, s->getWeight());
             if(categorySelection.isTight01)  varhisto01t->Fill(s->vars.reco1.pt, s->getWeight());
             if(categorySelection.isLoose01)  varhisto01l->Fill(s->vars.reco1.pt, s->getWeight());

             if(categorySelection.isVBFTight) varhistoVBFt->Fill(s->vars.reco2.pt, s->getWeight());
             if(categorySelection.isGGFTight) varhistoGGFt->Fill(s->vars.reco2.pt, s->getWeight());
             if(categorySelection.isVBFLoose) varhistoVBFl->Fill(s->vars.reco2.pt, s->getWeight());
             if(categorySelection.isTight01)  varhisto01t->Fill(s->vars.reco2.pt, s->getWeight());
             if(categorySelection.isLoose01)  varhisto01l->Fill(s->vars.reco2.pt, s->getWeight());
        }

        // recoMu_Eta
        if(varname.Contains("recoMu_Eta"))
        {
             varhistoall->Fill(s->vars.reco1.eta, s->getWeight());
             varhistoall->Fill(s->vars.reco2.eta, s->getWeight());

             if(categorySelection.isVBFTight) varhistoVBFt->Fill(s->vars.reco1.eta, s->getWeight());
             if(categorySelection.isGGFTight) varhistoGGFt->Fill(s->vars.reco1.eta, s->getWeight());
             if(categorySelection.isVBFLoose) varhistoVBFl->Fill(s->vars.reco1.eta, s->getWeight());
             if(categorySelection.isTight01)  varhisto01t->Fill(s->vars.reco1.eta, s->getWeight());
             if(categorySelection.isLoose01)  varhisto01l->Fill(s->vars.reco1.eta, s->getWeight());

             if(categorySelection.isVBFTight) varhistoVBFt->Fill(s->vars.reco2.eta, s->getWeight());
             if(categorySelection.isGGFTight) varhistoGGFt->Fill(s->vars.reco2.eta, s->getWeight());
             if(categorySelection.isVBFLoose) varhistoVBFl->Fill(s->vars.reco2.eta, s->getWeight());
             if(categorySelection.isTight01)  varhisto01t->Fill(s->vars.reco2.eta, s->getWeight());
             if(categorySelection.isLoose01)  varhisto01l->Fill(s->vars.reco2.eta, s->getWeight());
        }

        // NPV
        if(varname.Contains("NPV"))
        {
             varhistoall->Fill(s->vars.vertices.nVertices, s->getWeight());

             if(categorySelection.isVBFTight) varhistoVBFt->Fill(s->vars.vertices.nVertices, s->getWeight());
             if(categorySelection.isGGFTight) varhistoGGFt->Fill(s->vars.vertices.nVertices, s->getWeight());
             if(categorySelection.isVBFLoose) varhistoVBFl->Fill(s->vars.vertices.nVertices, s->getWeight());
             if(categorySelection.isTight01)  varhisto01t->Fill(s->vars.vertices.nVertices, s->getWeight());
             if(categorySelection.isLoose01)  varhisto01l->Fill(s->vars.vertices.nVertices, s->getWeight());
        }

        // jet_Pt
        if(varname.Contains("jet_Pt"))
        {
            float varvalue = -9999;

            for(unsigned int i=0; i<s->vars.validJets.size(); i++)
            {
                varvalue = s->vars.validJets[i].Pt();
                varhistoall->Fill(varvalue, s->getWeight());
            }

            if(categorySelection.isVBFTight)
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Pt();
                    varhistoVBFt->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelection.isGGFTight)
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Pt();
                    varhistoGGFt->Fill(varvalue, s->getWeight());
                }
            }
            if(categorySelection.isVBFLoose)
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Pt();
                    varhistoVBFl->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelection.isTight01)
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Pt();
                    varhisto01t->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelection.isLoose01)
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Pt();
                    varhisto01l->Fill(varvalue, s->getWeight());
                }
            }
        }

        // jet_Eta
        if(varname.Contains("jet_Eta"))
        {
            float varvalue = -9999;

            for(unsigned int i=0; i<s->vars.validJets.size(); i++)
            {
                varvalue = s->vars.validJets[i].Eta();
                varhistoall->Fill(varvalue, s->getWeight());
            }

            if(categorySelection.isVBFTight)
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Eta();
                    varhistoVBFt->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelection.isGGFTight)
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Eta();
                    varhistoGGFt->Fill(varvalue, s->getWeight());
                }
            }
            if(categorySelection.isVBFLoose)
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Eta();
                    varhistoVBFl->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelection.isTight01)
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Eta();
                    varhisto01t->Fill(varvalue, s->getWeight());
                }
            }

            if(categorySelection.isLoose01)
            {
                for(unsigned int i=0; i<s->vars.validJets.size(); i++)
                {
                    varvalue = s->vars.validJets[i].Eta();
                    varhisto01l->Fill(varvalue, s->getWeight());
                }
            }
        }

        // N_valid_jets
        if(varname.Contains("N_valid_jets"))
        {
             varhistoall->Fill(s->vars.validJets.size(), s->getWeight());

             if(categorySelection.isVBFTight) varhistoVBFt->Fill(s->vars.validJets.size(), s->getWeight());
             if(categorySelection.isGGFTight) varhistoGGFt->Fill(s->vars.validJets.size(), s->getWeight());
             if(categorySelection.isVBFLoose) varhistoVBFl->Fill(s->vars.validJets.size(), s->getWeight());
             if(categorySelection.isTight01)  varhisto01t->Fill(s->vars.validJets.size(), s->getWeight());
             if(categorySelection.isLoose01)  varhisto01l->Fill(s->vars.validJets.size(), s->getWeight());
        }

        // Dijet_mass
        if(varname.Contains("Dijet_mass"))
        {
             if(s->vars.validJets.size() >= 2)
             {
                 TLorentzVector dijet = s->vars.validJets[0] + s->vars.validJets[1];
                 varhistoall->Fill(dijet.M(), s->getWeight());

                 if(categorySelection.isVBFTight) varhistoVBFt->Fill(dijet.M(), s->getWeight());
                 if(categorySelection.isGGFTight) varhistoGGFt->Fill(dijet.M(), s->getWeight());
                 if(categorySelection.isVBFLoose) varhistoVBFl->Fill(dijet.M(), s->getWeight());
                 if(categorySelection.isTight01)  varhisto01t->Fill(dijet.M(), s->getWeight());
                 if(categorySelection.isLoose01)  varhisto01l->Fill(dijet.M(), s->getWeight());
             }
        }

        // Dijet_mass
        if(varname.Contains("dEta_jj"))
        {
             if(s->vars.validJets.size() >= 2)
             {
                 float dEta = s->vars.validJets[0].Eta() - s->vars.validJets[1].Eta();
                 varhistoall->Fill(dEta, s->getWeight());

                 if(categorySelection.isVBFTight) varhistoVBFt->Fill(dEta, s->getWeight());
                 if(categorySelection.isGGFTight) varhistoGGFt->Fill(dEta, s->getWeight());
                 if(categorySelection.isVBFLoose) varhistoVBFl->Fill(dEta, s->getWeight());
                 if(categorySelection.isTight01)  varhisto01t->Fill(dEta, s->getWeight());
                 if(categorySelection.isLoose01)  varhisto01l->Fill(dEta, s->getWeight());
             }
        }



        if(false)
          // ouput pt, mass info etc
          outputEvent(s->vars, categorySelection);

        // Reset the flags in preparation for the next event
        categorySelection.reset();
      }

      // Scale according to luminosity
      varhistoall->Scale(s->getScaleFactor(luminosity));
      varhistoVBFt->Scale(s->getScaleFactor(luminosity));
      varhistoGGFt->Scale(s->getScaleFactor(luminosity));
      varhistoVBFl->Scale(s->getScaleFactor(luminosity));
      varhisto01t->Scale(s->getScaleFactor(luminosity));
      varhisto01l->Scale(s->getScaleFactor(luminosity));

      // No trigger info in 80X MC samples, scale for trigger efficiency
      if(!s->sampleType.Contains("data"))
      {
          varhistoall->Scale(triggerSF);
          varhistoVBFt->Scale(triggerSF);
          varhistoGGFt->Scale(triggerSF);
          varhistoVBFl->Scale(triggerSF);
          varhisto01t->Scale(triggerSF);
          varhisto01l->Scale(triggerSF);
      }

      // Add to the appropriate list
      varlistall->Add(varhistoall);
      varlistVBFt->Add(varhistoVBFt);
      varlistGGFt->Add(varhistoGGFt);
      varlistVBFl->Add(varhistoVBFl);
      varlist01t->Add(varhisto01t);
      varlist01l->Add(varhisto01l);
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
    TCanvas* varstackcanvasall = dps->stackedHistogramsAndRatio(varlistall, "c_all_"+varname, varname+"_all_stack", varname, "Num Entries");
    TCanvas* varstackcanvasVBFt = dps->stackedHistogramsAndRatio(varlistVBFt, "c_vbft_"+varname, varname+"_vbft_stack", varname, "Num Entries");
    TCanvas* varstackcanvasVBFl = dps->stackedHistogramsAndRatio(varlistVBFl, "c_vbfl_"+varname, varname+"_vbfl_stack", varname, "Num Entries");
    TCanvas* varstackcanvasGGFt = dps->stackedHistogramsAndRatio(varlistGGFt, "c_ggft_"+varname, varname+"_ggft_stack", varname, "Num Entries");
    TCanvas* varstackcanvas01t = dps->stackedHistogramsAndRatio(varlist01t, "c_01t_"+varname, varname+"_01t_stack", varname, "Num Entries");
    TCanvas* varstackcanvas01l = dps->stackedHistogramsAndRatio(varlist01l, "c_01l_"+varname, varname+"_01l_stack", varname, "Num Entries");
    std::cout << std::endl;

    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile("rootfiles/validate_"+varname+"_x69_8_0_X_MC_categories_"+Form("%d",(int)luminosity)+".root", "RECREATE");
    TDirectory* stacks = savefile->mkdir("stacks");
    TDirectory* histos = savefile->mkdir("histos");

    // save the different histos in the appropriate directories in the tfile
    stacks->cd();
    varstackcanvasall->Write();
    varstackcanvasVBFt->Write();
    varstackcanvasGGFt->Write();
    varstackcanvasVBFl->Write();
    varstackcanvas01t->Write();
    varstackcanvas01l->Write();

    varstackcanvasall ->SaveAs("imgs/"+TString(varstackcanvasall->GetName())+".png");
    varstackcanvasVBFt->SaveAs("imgs/"+TString(varstackcanvasVBFt->GetName())+".png");
    varstackcanvasGGFt->SaveAs("imgs/"+TString(varstackcanvasGGFt->GetName())+".png");
    varstackcanvasVBFl->SaveAs("imgs/"+TString(varstackcanvasVBFl->GetName())+".png");
    varstackcanvas01t ->SaveAs("imgs/"+TString(varstackcanvas01t->GetName())+".png");
    varstackcanvas01l ->SaveAs("imgs/"+TString(varstackcanvas01l->GetName())+".png");

    histos->cd();
    varlistall->Write();
    varlistVBFt->Write();
    varlistVBFl->Write();
    varlistGGFt->Write();
    varlist01t->Write();
    varlist01l->Write();

    savefile->Close();

    return 0;
}
