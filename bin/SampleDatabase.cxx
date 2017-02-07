#include "Sample.h"
#include "DiMuPlottingSystem.h"

#include <sstream>
#include <map>
#include <vector>
#include <utility>

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

std::vector<Sample*>& GetSamples(std::map<TString, Sample*>& samples, TString location, TString select="ALL") {
  
  std::cout << "\n======== Getting samples: " << select << "\n" << std::endl;
  
  ///////////////////////////////////////////////////////////////////
  // SAMPLES---------------------------------------------------------
  ///////////////////////////////////////////////////////////////////
  
  TString in_dir;
  if (location == "UF")
    in_dir = "/cms/data/store/user/t2/users/acarnes/h2mumu/awb_samples/simplified/"; 
  else if (location == "CERN")
    in_dir = "root://eoscms.cern.ch//store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Jan28";
  else
    std::cout << "\n\nInput location is " << location << ", not UF or CERN.  NOT AN OPTION!!!\n\n" << std::endl;
  
  // ================================================================
  // Data -----------------------------------------------------------
  // ================================================================
  
  /////// INDIVIDUAL ERAS /////////////////////////////////////////////////////
  std::vector< std::tuple< TString, float, std::vector<int>, std::vector<TString> > > eras;
  // Era tuple has name, luminosity, numbers of files, and locations
  // Very rough lumi splitting between eras; needs to be updated - AWB 01.02.17
  eras.push_back( std::make_tuple("B", 5800, std::vector<int>{133},     std::vector<TString>{"B/170128_235241"}) );
  eras.push_back( std::make_tuple("C", 2600, std::vector<int>{ 44},     std::vector<TString>{"C/170128_235256"}) );
  eras.push_back( std::make_tuple("D", 4300, std::vector<int>{ 73},     std::vector<TString>{"D/170128_235314"}) );
  eras.push_back( std::make_tuple("E", 4100, std::vector<int>{ 62},     std::vector<TString>{"E/170128_235329"}) );
  eras.push_back( std::make_tuple("F", 3200, std::vector<int>{ 91, 46}, std::vector<TString>{"F_1/170130_063615", "F_2/170128_235359"}) );
  eras.push_back( std::make_tuple("G", 7800, std::vector<int>{107},     std::vector<TString>{"G/170128_235417"}) );
  //eras.push_back( std::make_tuple("H", 9014, std::vector<int>{116,  4}, std::vector<TString>{"H_1/170128_235433", "H_2/170128_235449"}) );
  
  for (auto era: eras) {
    if (select != "DATA" && select!="ALL" && select != "Run"+std::get<0>(era))
      continue;
    std::cout << "\nAdding files for Run" << std::get<0>(era) << " ...." << std::endl;
    
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"data/SingleMuon_SingleMu_2016"+std::get<0>(era)+".root") );
    } else {
      for (int i = 0; i < std::get<2>(era).size(); i++) {
    	for (int j = 1; j <= std::get<2>(era).at(i); j++) {
    	  in_file.Form( "%s/SingleMuon/SingleMu_2016%s/0000/tuple_%d.root", in_dir.Data(), std::get<3>(era).at(i).Data(), j );
    	  in_files.push_back(in_file);
    	}
      }
    } 
    
    Sample* data_sample = new Sample(in_files, "Run"+std::get<0>(era), "data");
    data_sample->lumi = std::get<1>(era);
    data_sample->xsec = 9999;
    samples["Run"+std::get<0>(era)] = data_sample;
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }
  
  // ================================================================
  // H2Mu_gg ---------------------------------------------------------
  // ================================================================
  
  if (select == "ALL" || select == "MC" || select == "SIGNAL" || select == "H2Mu_gg") {
    std::cout << "\nAdding files for H2Mu_gg ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"signal/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8_H2Mu_gg.root") );
    } else {
      for (int i = 1; i <= 1; i++) {
	in_file.Form( "%s/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_gg/170128_235508/0000/tuple_%d.root", in_dir.Data(), i);
	in_files.push_back(in_file);
      }
    }
    samples["H2Mu_gg"] = new Sample(in_files, "H2Mu_gg", "signal");
    samples["H2Mu_gg"]->xsec = 0.009618; // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  // ================================================================
  // H2Mu_VBF ---------------------------------------------------------
  // ================================================================
  
  if (select == "ALL" || select == "MC" || select == "SIGNAL" || select == "H2Mu_VBF") {
    std::cout << "\nAdding files for H2Mu_VBF ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"signal/VBF_HToMuMu_M125_13TeV_powheg_pythia8_H2Mu_VBF.root") );
    } else {
      for (int i = 1; i <= 1; i++) {
	in_file.Form( "%s/VBF_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_VBF/170128_235530/0000/tuple_%d.root", in_dir.Data(), i);
	in_files.push_back(in_file);
      }
    }
    samples["H2Mu_VBF"] = new Sample(in_files, "H2Mu_VBF", "signal");
    samples["H2Mu_VBF"]->xsec = 0.0008208; // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  } 
 
  // ================================================================
  // H2Mu_VH ---------------------------------------------------------
  // ================================================================
  
  if (select == "ALL" || select == "MC" || select == "SIGNAL" || select == "H2Mu_VH" || select == "H2Mu_ZH") {
    std::cout << "\nAdding files for H2Mu_ZH ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"signal/ZH_HToMuMu_M125_13TeV_powheg_pythia8_H2Mu_ZH.root") );
    } else {
      for (int i = 1; i <= 1; i++) {
	in_file.Form( "%s/ZH_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_ZH/170128_235617/0000/tuple_%d.root", in_dir.Data(), i);
	in_files.push_back(in_file);
      }
    }
    samples["H2Mu_ZH"] = new Sample(in_files, "H2Mu_ZH", "signal");
    samples["H2Mu_ZH"]->xsec = 0.0002136;     // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  } 
 
  if (select == "ALL" || select == "MC" || select == "SIGNAL" || select == "H2Mu_VH" || select == "H2Mu_WH"  || select == "H2Mu_WH_pos") {
    std::cout << "\nAdding files for H2Mu_WH_pos ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"signal/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8_H2Mu_WH_pos.root") );
    } else {
      for (int i = 1; i <= 1; i++) {
	in_file.Form( "%s/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_WH_pos/170128_235546/0000/tuple_%d.root", in_dir.Data(), i);
	in_files.push_back(in_file);
      }
    }
    samples["H2Mu_WH_pos"] = new Sample(in_files, "H2Mu_WH_pos", "signal");
    samples["H2Mu_WH_pos"]->xsec = 0.0001858; // 0.851*0.0002176;     // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  } 
 
  if (select == "ALL" || select == "MC" || select == "SIGNAL" || select == "H2Mu_VH" || select == "H2Mu_WH" || select == "H2Mu_WH_neg") {
    std::cout << "\nAdding files for H2Mu_WH_neg ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"signal/WMinusH_HToMuMu_M125_13TeV_powheg_pythia8_H2Mu_WH_neg.root") );
    } else {
      for (int i = 1; i <= 1; i++) {
	in_file.Form( "%s/WMinusH_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_WH_neg/170128_235601/0000/tuple_%d.root", in_dir.Data(), i);
	in_files.push_back(in_file);
      }
    }
    samples["H2Mu_WH_neg"] = new Sample(in_files, "H2Mu_WH_neg", "signal");
    samples["H2Mu_WH_neg"]->xsec = 0.0001164; //0.5331*0.0002176;     // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  } 

 
  // ================================================================
  // DYJetsToLL -----------------------------------------------------
  // ================================================================
  
  float DY_xsec_m50 = 5765; // pb // old value = 6025.2
  float DY_m100to200_factor = 1.235; 
  
  if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "ZJets" || select == "ZJets_AMC") {
    std::cout << "\nAdding files for ZJets_AMC ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ZJets_AMC.root") );
    } else {
      for (int i = 1; i <= 15; i++) {
	in_file.Form( "%s/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ZJets_AMC/170128_235652/0000/tuple_%d.root", in_dir.Data(), i);
	in_files.push_back(in_file);
      }
    }
    samples["ZJets_AMC"] = new Sample(in_files, "ZJets_AMC", "background");
    samples["ZJets_AMC"]->xsec = DY_xsec_m50; 
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  if (select == "MC" || select == "BACKGROUND" || select == "ZJets" || select == "ZJets_MG") {
    std::cout << "\nAdding files for ZJets_MG ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ZJets_MG.root") );
    } else {
      for (int i = 1; i <= 48; i++) {
	in_file.Form( "%s/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ZJets_MG/170130_063203/0000/tuple_%d.root", in_dir.Data(), i);
	in_files.push_back(in_file);
      }
    }
    samples["ZJets_MG"] = new Sample(in_files, "ZJets_MG", "background");
    samples["ZJets_MG"]->xsec = DY_xsec_m50;
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  if (select == "MC" || select == "BACKGROUND" || select == "ZJets" || select == "ZJets_hiM") {
    std::cout << "\nAdding files for ZJets_hiM ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ZJets_hiM.root") );
    } else {
      for (int i = 1; i <= 2; i++) {
	in_file.Form( "%s/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ZJets_hiM/170128_235708/0000/tuple_%d.root", in_dir.Data(), i);
	in_files.push_back(in_file);
      }
    }
    samples["ZJets_hiM"] = new Sample(in_files, "ZJets_hiM", "background");
    samples["ZJets_hiM"]->xsec = DY_xsec_m50 * DY_m100to200_factor; // 7117
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  if (select == "MC" || select == "BACKGROUND" || select == "ZJets" || select == "ZJets_hiM_SpringPU") {
    std::cout << "\nAdding files for ZJets_hiM_SpringPU ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"dy/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ZJets_hiM_SpringPU.root") );
    } else {
      for (int i = 1; i <= 2; i++) {
	in_file.Form( "%s/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ZJets_hiM_SpringPU/170130_062620/0000/tuple_%d.root", in_dir.Data(), i);
	in_files.push_back(in_file);
      }
    }
    samples["ZJets_hiM_SpringPU"] = new Sample(in_files, "ZJets_hiM_SpringPU", "background");
    samples["ZJets_hiM_SpringPU"]->xsec = DY_xsec_m50 * DY_m100to200_factor; // 7117
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }


  // ================================================================
  // TTJets ---------------------------------------------------------
  // ================================================================
  
  if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "ttbar" || select == "tt_ll_AMC") {
    std::cout << "\nAdding files for tt_ll_AMC ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"ttjets/TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8_tt_ll_AMC.root") );
    } else {
      for (int i = 1; i <= 10; i++) {
	in_file.Form( "%s/TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/tt_ll_AMC/170128_235902/0000/tuple_%d.root", in_dir.Data(), i);
	in_files.push_back(in_file);
      }
    }
    samples["tt_ll_AMC"] = new Sample(in_files, "tt_ll_AMC", "background");
    samples["tt_ll_AMC"]->xsec = 85.656; // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  if (select == "MC" || select == "BACKGROUND" || select == "ttbar" || select == "tt_ll_MG") {
    std::cout << "\nAdding files for tt_ll_MG ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"ttjets/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_tt_ll_MG_2.root") );
    } else {
      for (int i = 1; i <= 5; i++) {
	in_file.Form( "%s/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/tt_ll_MG_1/170128_235829/0000/tuple_%d.root", in_dir.Data(), i);
	in_files.push_back(in_file);
      }
      for (int i = 1; i <= 18; i++) {
	in_file.Form( "%s/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/tt_ll_MG_2/170128_235847/0000/tuple_%d.root", in_dir.Data(), i);
	in_files.push_back(in_file);
      }
    }
    samples["tt_ll_MG"] = new Sample(in_files, "tt_ll_MG", "background");
    samples["tt_ll_MG"]->xsec = 85.656; // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }


  // ================================================================
  // SingleTop ------------------------------------------------------
  // ================================================================
  
  if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "singleTop" || select == "tW" || select == "tW_pos") {
    std::cout << "\nAdding files for tW_pos ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"singletop/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1_tW_pos.root") );
    } else {
      // for (int i = 1; i <= 1; i++) {
      // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
      // 	in_files.push_back(in_file);
      // }
    }
    samples["tW_pos"] = new Sample(in_files, "tW_pos", "background");
    samples["tW_pos"]->xsec = 35.85; // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "singleTop" || select == "tW" || select == "tW_neg") {
    std::cout << "\nAdding files for tW_neg ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"singletop/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1_tW_neg.root") );
    } else {
      // for (int i = 1; i <= 1; i++) {
      // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
      // 	in_files.push_back(in_file);
      // }
    }
    samples["tW_neg"] = new Sample(in_files, "tW_neg", "background");
    samples["tW_neg"]->xsec = 35.85; // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "singleTop" || select == "tZq") {
    std::cout << "\nAdding files for tZq ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"singletop/tZq_ll_4f_13TeV-amcatnlo-pythia8_tZq.root") );
    } else {
      // for (int i = 1; i <= 1; i++) {
      // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
      // 	in_files.push_back(in_file);
    }
    samples["tZq"] = new Sample(in_files, "tZq", "background");
    samples["tZq"]->xsec = 0.0758; // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  //if (select == "MC" || select == "BACKGROUND" || select == "singleTop" || select == "tZW") {
  //  std::cout << "\nAdding files for tZW ..." << std::endl;
  //  std::vector<TString> in_files;
  //  TString in_file;
  //  if (location == "UF") {
  //    in_files.push_back( TString(in_dir+"singletop/ST_tWll_5f_LO_13TeV-MadGraph-pythia8_tZW.root") );
  //  } else {
  //    // for (int i = 1; i <= 1; i++) {
  //    // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
  //    // 	in_files.push_back(in_file);
  //    // }
  //  }
  //  samples["tZW"] = new Sample(in_files, "tZW", "background");
  //  samples["tZW"]->xsec = -999; // pb
  //  std::cout << ".... " << in_files.size() << " files added." << std::endl;
  //}

    // ================================================================
    // TTV ---------------------------------------------------------
    // ================================================================
    
  if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "ttX" || select == "ttV" || select == "ttW") {
    std::cout << "\nAdding files for ttW ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"ttv/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_ttW.root") );
    } else {
      // for (int i = 1; i <= 1; i++) {
      // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
      // 	in_files.push_back(in_file);
      // }
    }
    samples["ttW"] = new Sample(in_files, "ttW", "background");
    samples["ttW"]->xsec = 0.2043; // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "ttX" || select == "ttV" || select == "ttZ") {
    std::cout << "\nAdding files for ttZ ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"ttv/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_ttZ.root") );
    } else {
      // for (int i = 1; i <= 1; i++) {
      // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
      // 	in_files.push_back(in_file);
      // }
    }
    samples["ttZ"] = new Sample(in_files, "ttZ", "background");
    samples["ttZ"]->xsec = 0.2529; // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  // ================================================================
  // Diboson ------------------------------------------------------
  // ================================================================
  
  if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "VV" || select == "WW") {
    std::cout << "\nAdding files for WW ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"diboson/WWTo2L2Nu_13TeV-powheg_WW.root") );
    } else {
      // for (int i = 1; i <= 1; i++) {
      // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
      // 	in_files.push_back(in_file);
      // }
    }
    samples["WW"] = new Sample(in_files, "WW", "background");
    samples["WW"]->xsec = 12.46; // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

   //if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "VV" || select == "WZ" || select == "WZ_2l") {
   //  std::cout << "\nAdding files for WZ_2l ..." << std::endl;
   //  std::vector<TString> in_files;
   //  TString in_file;
   //  if (location == "UF") {
   //    in_files.push_back( TString(in_dir+"diboson/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_WZ_2l.root") );
   //  } else {
   //    // for (int i = 1; i <= 1; i++) {
   //    // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
   //    // 	in_files.push_back(in_file);
   //    // }
   //  }
   //  samples["WZ_2l"] = new Sample(in_files, "WZ_2l", "background");
   //  samples["WZ_2l"]->xsec = ; // pb
   //  std::cout << ".... " << in_files.size() << " files added." << std::endl;
   //}

  if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "VV" || select == "WZ" || select == "WZ_3l") {
    std::cout << "\nAdding files for WZ_3l ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"diboson/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_WZ_3l_AMC.root") );
    } else {
      // for (int i = 1; i <= 1; i++) {
      // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
      // 	in_files.push_back(in_file);
      // }
    }
    samples["WZ_3l"] = new Sample(in_files, "WZ_3l", "background");
    samples["WZ_3l"]->xsec = 2.113; // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "VV" || select == "ZZ" || select == "ZZ_2l_2v") {
    std::cout << "\nAdding files for ZZ_2l_2v ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"diboson/ZZTo2L2Nu_13TeV_powheg_pythia8_ZZ_2l_2v.root") );
    } else {
      // for (int i = 1; i <= 1; i++) {
      // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
      // 	in_files.push_back(in_file);
      // }
    }
    samples["ZZ_2l_2v"] = new Sample(in_files, "ZZ_2l_2v", "background");
    samples["ZZ_2l_2v"]->xsec = 0.564; // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "VV" || select == "ZZ" || select == "ZZ_2l_2q") {
    std::cout << "\nAdding files for ZZ_2l_2q ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"diboson/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_ZZ_2l_2q.root") );
    } else {
      // for (int i = 1; i <= 1; i++) {
      // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
      // 	in_files.push_back(in_file);
      // }
    }
    samples["ZZ_2l_2q"] = new Sample(in_files, "ZZ_2l_2q", "background");
    samples["ZZ_2l_2q"]->xsec = 3.22; // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "VV" || select == "ZZ" || select == "ZZ_4l") {
    std::cout << "\nAdding files for ZZ_4l ..." << std::endl;
    std::vector<TString> in_files;
    TString in_file;
    if (location == "UF") {
      in_files.push_back( TString(in_dir+"diboson/ZZTo4L_13TeV-amcatnloFXFX-pythia8_ZZ_4l_AMC.root") );
    } else {
      // for (int i = 1; i <= 1; i++) {
      // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
      // 	in_files.push_back(in_file);
      // }
    }
    samples["ZZTo4L"] = new Sample(in_files, "ZZTo4L", "background");
    samples["ZZTo4L"]->xsec = 1.212; // pb
    std::cout << ".... " << in_files.size() << " files added." << std::endl;
  }

  // ================================================================
  // Triboson ------------------------------------------------------
  // ================================================================
  
  // if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "VVV" || select == "WWW") {
  //   std::cout << "\nAdding files for WWW ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location == "UF") {
  //     in_files.push_back( TString(in_dir+"triboson/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8_WWW.root") );
  //   } else {
  //     // for (int i = 1; i <= 1; i++) {
  //     // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
  //     // 	in_files.push_back(in_file);
  //     // }
  //   }
  //   samples["WWW"] = new Sample(in_files, "WWW", "background");
  //   samples["WWW"]->xsec = -999; // pb
  //   std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

  // if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "VVV" || select == "WWZ") {
  //   std::cout << "\nAdding files for WWZ ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location == "UF") {
  //     in_files.push_back( TString(in_dir+"triboson/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_WWZ.root") );
  //   } else {
  //     // for (int i = 1; i <= 1; i++) {
  //     // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
  //     // 	in_files.push_back(in_file);
  //     // }
  //   }
  //   samples["WWZ"] = new Sample(in_files, "WWZ", "background");
  //   samples["WWZ"]->xsec = -999; // pb
  //   std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

  // if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "VVV" || select == "WZZ") {
  //   std::cout << "\nAdding files for WZZ ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location == "UF") {
  //     in_files.push_back( TString(in_dir+"triboson/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_WZZ.root") );
  //   } else {
  //     // for (int i = 1; i <= 1; i++) {
  //     // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
  //     // 	in_files.push_back(in_file);
  //     // }
  //   }
  //   samples["WZZ"] = new Sample(in_files, "WZZ", "background");
  //   samples["WZZ"]->xsec = -999; // pb
  //   std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

  // if (select == "ALL" || select == "MC" || select == "BACKGROUND" || select == "VVV" || select == "ZZZ") {
  //   std::cout << "\nAdding files for ZZZ ..." << std::endl;
  //   std::vector<TString> in_files;
  //   TString in_file;
  //   if (location == "UF") {
  //     in_files.push_back( TString(in_dir+"triboson/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_ZZZ.root") );
  //   } else {
  //     // for (int i = 1; i <= 1; i++) {
  //     // 	in_file.Form( "%s///0000/tuple_%d.root", in_dir.Data(), i);
  //     // 	in_files.push_back(in_file);
  //     // }
  //   }
  //   samples["ZZZ"] = new Sample(in_files, "ZZZ", "background");
  //   samples["ZZZ"]->xsec = -999; // pb
  //   std::cout << ".... " << in_files.size() << " files added." << std::endl;
  // }

}
