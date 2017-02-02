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
    std::cout << std::endl;
    std::cout << "++++++++ Luminosity = " << luminosity << std::endl;
    std::cout << "======== Samples Available... " << std::endl;
    std::cout << std::endl;

    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////
    
    TString sampledir = "/cms/data/store/user/t2/users/acarnes/h2mumu/awb_samples/simplified/"; 

    // ================================================================
    // Data -----------------------------------------------------------
    // ================================================================

    /////// INDIVIDUAL RUNS /////////////////////////////////////////////////////
    // Should split F, G, H in half
    //std::vector<TString> runs = {"B", "C", "D", "E", "F1", "F2", "G1", "G2", "H1", "H2"};
    std::vector<TString> runs = {"B", "C", "D", "E", "F", "G", "H"};
    for(auto X: runs)
    {
      TString datafilenameX = TString(sampledir+"data/SingleMuon_SingleMu_2016"+X+".root");
      std::cout << datafilenameX << std::endl;

      Sample* datasample = new Sample(datafilenameX, "Run"+X, "data");
      datasample->lumi = luminosity;
      datasample->xsec = 9999;
      datasample->pileupfile = "";
      samples["Run"+X] = datasample;
    }

    // ================================================================
    // VBF ---------------------------------------------------------
    // ================================================================

    TString vbffilename   = TString(sampledir+"signal/VBF_HToMuMu_M125_13TeV_powheg_pythia8_H2Mu_VBF.root");
    std::cout << vbffilename << std::endl;
    samples["VBF"] = new Sample(vbffilename, "VBF", "signal");
    samples["VBF"]->pileupfile = ""; //nPU
    samples["VBF"]->xsec = 0.0008208; // pb

    // ================================================================
    // GGF ---------------------------------------------------------
    // ================================================================

    TString ggfilename   = TString(sampledir+"signal/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8_H2Mu_gg.root");
    std::cout << ggfilename << std::endl;
    samples["GGF"] = new Sample(ggfilename, "GGF", "signal");
    samples["GGF"]->pileupfile = "";      //nPU
    samples["GGF"]->xsec = 0.009618; // pb

    // ================================================================
    // VH ---------------------------------------------------------
    // ================================================================

    TString zhfilename   = TString(sampledir+"signal/ZH_HToMuMu_M125_13TeV_powheg_pythia8_H2Mu_ZH.root");
    std::cout << zhfilename << std::endl;
    samples["ZH"] = new Sample(zhfilename, "ZH", "signal");
    samples["ZH"]->pileupfile = "";      //nPU
    samples["ZH"]->xsec = 0.0002136;     // pb

    TString wphfilename   = TString(sampledir+"signal/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8_H2Mu_WH_pos.root");
    std::cout << wphfilename << std::endl;
    samples["WplusH"] = new Sample(wphfilename, "WplusH", "signal");
    samples["WplusH"]->pileupfile = "";      //nPU
    samples["WplusH"]->xsec = 0.0001858; // 0.851*0.0002176;     // pb

    TString wmhfilename   = TString(sampledir+"signal/WMinusH_HToMuMu_M125_13TeV_powheg_pythia8_H2Mu_WH_neg.root");
    std::cout << wmhfilename << std::endl;
    samples["WminusH"] = new Sample(zhfilename, "WminusH", "signal");
    samples["WminusH"]->pileupfile = "";      //nPU
    samples["WminusH"]->xsec = 0.0001164; //0.5331*0.0002176;     // pb

    // ================================================================
    // DYJetsToLL -----------------------------------------------------
    // ================================================================

    float dyxsec_m50 = 5765; // pb // old value = 6025.2
    float dy_m100to200_factor = 1.235; 

    TString dyfile   = TString(sampledir+"dy/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ZJets_AMC.root");
    std::cout << dyfile << std::endl;
    samples["DYJetsToLL"] = new Sample(dyfile, "DYJetsToLL", "background");
    samples["DYJetsToLL"]->pileupfile = ""; //nPU
    samples["DYJetsToLL"]->xsec = dyxsec_m50; 

    //dyfile   = TString(sampledir+"dy/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ZJets_MG.root");
    //samples["DYJetsToLL_madgraph"] = new Sample(dyfile, "DYJetsToLL_madgraph", "background");
    //samples["DYJetsToLL_madgraph"]->pileupfile = ""; //nPU
    //samples["DYJetsToLL_madgraph"]->xsec = dyxsec_m50;

    //dyfile   = TString(sampledir+"dy/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ZJets_hiM.root");
    //samples["DYJetsToLL_M100to200"] = new Sample(dyfile, "DYJetsToLL_M100to200", "background");
    //samples["DYJetsToLL_M100to200"]->pileupfile = ""; //nPU
    //samples["DYJetsToLL_M100to200"]->xsec = dyxsec_m50; // 7117

    //dyfile = TString(sampledir+"dy/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ZJets_hiM_SpringPU.root");
    //samples["DYJetsToLL_M100to200_SpringPU"] = new Sample(dyfile, "DYJetsToLL_M100to200_SpringPU", "background");
    //samples["DYJetsToLL_M100to200_SpringPU"]->pileupfile = ""; //nPU
    //samples["DYJetsToLL_M100to200_SpringPU"]->xsec = dyxsec_m50; // 7117

    // ================================================================
    // TTJets ---------------------------------------------------------
    // ================================================================

    TString ttbarfilename   = TString(sampledir+"ttjets/TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8_tt_ll_AMC.root");
    std::cout << ttbarfilename << std::endl;
    samples["TTJetsToLL"] = new Sample(ttbarfilename, "TTJetsToLL", "background");
    samples["TTJetsToLL"]->pileupfile = ""; //nPU
    samples["TTJetsToLL"]->xsec = 85.656; // pb

    //ttbarfilename   = TString(sampledir+"ttjets/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_tt_ll_MG_2.root");
    //samples["TTJetsToLL_madgraph"] = new Sample(ttbarfilename, "TTJetsToLL_madgraph", "background");
    //samples["TTJetsToLL_madgraph"]->pileupfile = ""; //nPU
    //samples["TTJetsToLL_madgraph"]->xsec = 85.656; // pb

    // ================================================================
    // SingleTop ------------------------------------------------------
    // ================================================================

    //TString tbarwfilename   = TString(sampledir+"singletop/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1_tW_neg.root");
    //std::cout << tbarwfilename << std::endl;
    //samples["tbarW"] = new Sample(tbarwfilename, "tbarW", "background");
    //samples["tbarW"]->pileupfile = ""; //nPU
    //samples["tbarW"]->xsec = 35.85; // pb

    //TString twfilename   = TString(sampledir+"singletop/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1_tW_pos.root");
    //std::cout << twfilename << std::endl;
    //samples["tW"] = new Sample(twfilename, "tW", "background");
    //samples["tW"]->pileupfile = ""; //nPU
    //samples["tW"]->xsec = 35.85; // pb

    //TString twzfilename   = TString(sampledir+"singletop/ST_tWll_5f_LO_13TeV-MadGraph-pythia8_tZW.root");
    //samples["tWZ"] = new Sample(twzfilename, "tWZ", "background");
    //samples["tWZ"]->pileupfile = ""; //nPU
    //samples["tWZ"]->xsec = -999; // pb

    //TString tzqfilename   = TString(sampledir+"singletop/tZq_ll_4f_13TeV-amcatnlo-pythia8_tZq.root");
    //std::cout << tzqfilename << std::endl;
    //samples["tZq"] = new Sample(tzqfilename, "tZq", "background");
    //samples["tZq"]->pileupfile = ""; //nPU
    //samples["tZq"]->xsec = 0.0758; // pb

    // ================================================================
    // TTV ---------------------------------------------------------
    // ================================================================
    
    //TString ttwfilename   = TString(sampledir+"ttv/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_ttW.root");
    //std::cout << ttwfilename << std::endl;
    //samples["ttW"] = new Sample(ttwfilename, "ttW", "background");
    //samples["ttW"]->pileupfile = ""; //nPU
    //samples["ttW"]->xsec = 0.2043; // pb

    //TString ttzfilename   = TString(sampledir+"ttv/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_ttZ.root");
    //std::cout << ttzfilename << std::endl;
    //samples["ttZ"] = new Sample(ttzfilename, "ttZ", "background");
    //samples["ttZ"]->pileupfile = ""; //nPU
    //samples["ttZ"]->xsec = 0.2529; // pb

    // ================================================================
    // Diboson ------------------------------------------------------
    // ================================================================

    //TString wwfilename   = TString(sampledir+"diboson/WWTo2L2Nu_13TeV-powheg_WW.root");
    //std::cout << wwfilename << std::endl;
    //samples["WWTo2L2Nu"] = new Sample(wwfilename, "WWTo2L2Nu", "background");
    //samples["WWTo2L2Nu"]->pileupfile = ""; //nPU
    //samples["WWTo2L2Nu"]->xsec = 12.46; // pb

    //TString wzfilename   = TString(sampledir+"diboson/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_WZ_2l.root");
    //samples["WZTo2L2Q"] = new Sample(wzfilename, "WZTo2L2Q", "background");
    //samples["WZTo2L2Q"]->pileupfile = ""; //nPU
    //samples["WZTo2L2Q"]->xsec = ; // pb

    //TString wz2filename   = TString(sampledir+"diboson/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_WZ_3l_AMC.root");
    //std::cout << wz2filename << std::endl;
    //samples["WZTo3LNu"] = new Sample(wz2filename, "WZTo3LNu", "background");
    //samples["WZTo3LNu"]->pileupfile = ""; //nPU
    //samples["WZTo3LNu"]->xsec = 2.113; // pb

    //TString zz2l2nufilename   = TString(sampledir+"diboson/ZZTo2L2Nu_13TeV_powheg_pythia8_ZZ_2l_2v.root");
    //std::cout << zz2l2nufilename << std::endl;
    //samples["ZZTo2L2Nu"] = new Sample(zz2l2nufilename, "ZZTo2L2Nu", "background");
    //samples["ZZTo2L2Nu"]->pileupfile = ""; //nPU
    //samples["ZZTo2L2Nu"]->xsec = 0.564; // pb

    //TString zz2l2qfilename   = TString(sampledir+"diboson/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_ZZ_2l_2q.root");
    //std::cout << zz2l2qfilename << std::endl;
    //samples["ZZTo2L2Q"] = new Sample(zz2l2qfilename, "ZZTo2L2Q", "background");
    //samples["ZZTo2L2Q"]->pileupfile = ""; //nPU
    //samples["ZZTo2L2Q"]->xsec = 3.22; // pb

    //TString zz4lfilename   = TString(sampledir+"diboson/ZZTo4L_13TeV-amcatnloFXFX-pythia8_ZZ_4l_AMC.root");
    //std::cout << zz4lfilename << std::endl;
    //samples["ZZTo4L"] = new Sample(zz4lfilename, "ZZTo4L", "background");
    //samples["ZZTo4L"]->pileupfile = ""; //nPU
    //samples["ZZTo4L"]->xsec = 1.212; // pb

    // ================================================================
    // Triboson ------------------------------------------------------
    // ================================================================

    //TString wwwfilename   = TString(sampledir+"triboson/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8_WWW.root");
    //samples["WWW"] = new Sample(wwwfilename, "WWW", "background");
    //samples["WWW"]->pileupfile = ""; //nPU
    //samples["WWW"]->xsec = -999; // pb

    //TString wwzfilename   = TString(sampledir+"triboson/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_WWZ.root");
    //samples["WWZ"] = new Sample(wwzfilename, "WWZ", "background");
    //samples["WWZ"]->pileupfile = ""; //nPU
    //samples["WWZ"]->xsec = -999; // pb

    //TString wzzfilename   = TString(sampledir+"triboson/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_WZZ.root");
    //samples["WZZ"] = new Sample(wzzfilename, "WZZ", "background");
    //samples["WZZ"]->pileupfile = ""; //nPU
    //samples["WZZ"]->xsec = -999; // pb

    //TString zzzfilename   = TString(sampledir+"triboson/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_ZZZ.root");
    //samples["ZZZ"] = new Sample(zzzfilename, "ZZZ", "background");
    //samples["ZZZ"]->pileupfile = ""; //nPU
    //samples["ZZZ"]->xsec = -999; // pb
}
