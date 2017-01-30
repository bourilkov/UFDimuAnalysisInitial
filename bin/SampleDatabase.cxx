#include "Sample.h"
#include "DiMuPlottingSystem.h"

#include <sstream>
#include <map>
#include <vector>
#include <utility>

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

std::vector<Sample*>& getSamples(float luminosity, std::map<TString, Sample*>& samples) 
{
    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // ================================================================
    // Data -----------------------------------------------------------
    // ================================================================

    //TString datafilename = 
    //TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data/25ns/golden/CMSSW_8_0_X/stage_1_singleMuon_Run2016BCDEFGH_ALL_etm.root");

    //Sample* datasample = new Sample(datafilename, "Data", "data");
    //datasample->lumi = luminosity;
    //datasample->xsec = 9999;
    //datasample->pileupfile = "pu_reweight_trees/8_0_X/PU_2016BCDEFGH_xsec69p2mb_CMSSW_8_0_X.root";
    //samples["Data"] = datasample;

    /////// INDIVIDUAL RUNS /////////////////////////////////////////////////////
    std::vector<TString> runs = {"B", "C", "D", "E", "F", "G", "H"};
    for(auto X: runs)
    {
      TString datafilenameX = 
      TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data/25ns/golden/CMSSW_8_0_X/individual_runs/stage1_singleMuon_2016"+X+"_ALL_etm.root");

      Sample* datasample = new Sample(datafilenameX, "Run"+X, "data");
      datasample->lumi = luminosity;
      datasample->xsec = 9999;
      datasample->pileupfile = "pu_reweight_trees/8_0_X/PU_2016BCDEFGH_xsec69p2mb_CMSSW_8_0_X.root";
      samples["Run"+X] = datasample;
    }

    // ================================================================
    // DYJetsToLL -----------------------------------------------------
    // ================================================================

    TString dyfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/bg/dy/CMSSW_8_0_X/stage_1_dy_jetsToLL_ALL_etm.root");
    samples["DYJetsToLL"] = new Sample(dyfilename, "DYJetsToLL", "background");
    samples["DYJetsToLL"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_DYJetsToLL.root"; //nPU
    samples["DYJetsToLL"]->xsec = 6025.2; // pb

    // ================================================================
    // TTJets ---------------------------------------------------------
    // ================================================================

    TString ttbarfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/bg/ttbar/CMSSW_8_0_X/stage_1_ttJets_ALL_etm.root");
    samples["TTJets"] = new Sample(ttbarfilename, "TTJets", "background");
    samples["TTJets"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_TTJets.root"; //nPU
    samples["TTJets"]->xsec = 831.76; // pb

    // ================================================================
    // VBF ---------------------------------------------------------
    // ================================================================

    TString vbffilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/signal/CMSSW_8_0_X/stage_1_vbf_HToMuMu_ALL_etm.root");
    samples["VBF"] = new Sample(vbffilename, "VBF", "signal");
    samples["VBF"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_VBF.root"; //nPU
    samples["VBF"]->xsec = 3.727*0.00022; // pb

    // ================================================================
    // GGF ---------------------------------------------------------
    // ================================================================

    TString ggfilename   = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/mc/signal/CMSSW_8_0_X/stage_1_gg_HToMuMu_ALL_etm.root");
    samples["GGF"] = new Sample(ggfilename, "GGF", "signal");
    samples["GGF"]->pileupfile = "./pu_reweight_trees/8_0_X/PUCalib_GGF.root"; //nPU
    samples["GGF"]->xsec = 43.62*0.00022; // pb

    TString ggAWBfilename   = TString("HToMuMu_NTuple_1.root");
    samples["GGF_AWB"] = new Sample(ggAWBfilename, "GGF_AWB", "signal");
    samples["GGF_AWB"]->pileupfile = "";      //nPU
    samples["GGF_AWB"]->xsec = 43.62*0.00022; // pb

    ///////////////////////////////////////////////////////////////////
    // Output Sample Info ---------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // Loop through all of the samples to output the information
    std::cout << std::endl;
    std::cout << "++++++++ Luminosity = " << luminosity << std::endl;
    std::cout << "======== Samples Available... " << std::endl;
    std::cout << std::endl;

    //makePUHistos(samples);
    
    for(auto &i : samples)
    {
        // Output info about the current file
        std::cout << std::endl;
        std::cout << "    sample name:       " << i.second->name << std::endl;
        std::cout << "    sample file:       " << i.second->filename << std::endl;
        std::cout << "    pileup file:       " << i.second->pileupfile << std::endl;
        std::cout << "    nOriginalWeighted: " << i.second->nOriginalWeighted << std::endl;
        std::cout << "    nOriginal:         " << i.second->nOriginal << std::endl;
        std::cout << "    N:                 " << i.second->N << std::endl;
        std::cout << std::endl;
    }
}
