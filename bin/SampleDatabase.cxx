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
    
    TString sampledir = "/cms/data/store/user/t2/users/acarnes/h2mumu/awb_samples/simplified/"; 

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

    float dyxsec_m50 = 6025.2; //5765;           // pb // old value = 6025.2
    float dy_m100to200_factor = 1.235; 

    TString dyfile   = TString(sampledir+"dy/DYJetsToLL_M-50_amcatnlo.root");
    samples["DYJetsToLL"] = new Sample(dyfile, "DYJetsToLL", "background");
    samples["DYJetsToLL"]->pileupfile = ""; //nPU
    samples["DYJetsToLL"]->xsec = dyxsec_m50; 

    dyfile   = TString(sampledir+"dy/DYJetsToLL_M-50_madgraph.root");
    samples["DYJetsToLL_madgraph"] = new Sample(dyfile, "DYJetsToLL_madgraph", "background");
    samples["DYJetsToLL_madgraph"]->pileupfile = ""; //nPU
    samples["DYJetsToLL_madgraph"]->xsec = dyxsec_m50;

    dyfile   = TString(sampledir+"dy/DYJetsToLL_M-100to200_amcatnlo.root");
    samples["DYJetsToLL_M100to200"] = new Sample(dyfile, "DYJetsToLL_M100to200", "background");
    samples["DYJetsToLL_M100to200"]->pileupfile = ""; //nPU
    samples["DYJetsToLL_M100to200"]->xsec = dyxsec_m50; // 7117

    dyfile = TString(sampledir+"/DYJetsToLL_M-100to200_amcatnlo_SpringPU.root");
    samples["DYJetsToLL_M100to200_SpringPU"] = new Sample(dyfile, "DYJetsToLL_M100to200_SpringPU", "background");
    samples["DYJetsToLL_M100to200_SpringPU"]->pileupfile = ""; //nPU
    samples["DYJetsToLL_M100to200_SpringPU"]->xsec = dyxsec_m50; // 7117

    // ================================================================
    // TTJets ---------------------------------------------------------
    // ================================================================

    TString ttbarfilename   = TString(sampledir+"ttjets/TTJetsToLL_amcatnlo.root");
    samples["TTJetsToLL"] = new Sample(ttbarfilename, "TTJetsToLL", "background");
    samples["TTJetsToLL"]->pileupfile = ""; //nPU
    samples["TTJetsToLL"]->xsec = 85.656; // pb

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

    TString ggfilename   = TString("HToMuMu_NTuple_1.root");
    samples["GGF"] = new Sample(ggfilename, "GGF_AWB", "signal");
    samples["GGF"]->pileupfile = "";      //nPU
    samples["GGF"]->xsec = 43.62*0.00022; // pb

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
