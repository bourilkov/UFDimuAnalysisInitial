
/////////////////////////////////////////////////////////////////////////////
///  Higgs vs. background classification for multiple factories and MVAs  ///
///                      Andrew Brinkerhoff 23.01.17                      ///
///                                                                       ///
///  Adapted from ROOT TMVAClassification.C                               ///
///  Run using "root -l HiggsClassification_v0.C                          /// 
/////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <cassert>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVAGui.h"

// Specific includes for UFDimuAnalysis
#include "Sample.h"
#include "SampleDatabase.cxx"
#include "EventSelection.h"
#include "MuonSelection.h"
#include "CategorySelection.h"

// Extra tools - AWB 13.03.17
#include "TMVA_helper.h"

const int MAX_EVT    = 1000;
const int MAX_SIG    =  500;
const int MAX_BKG    =  500;
const int REPORT_EVT =   10;

const int MAX_TR_SIG = 1000;
const int MAX_TR_BKG = 1000;

const double PI = 3.14159265359;
const double BIT = 0.000001; // Tiny value or offset

using namespace TMVA;

void TMVAClassification_H2Mu ( TString myMethodList = "" ) {

   // TMVA::Tools::Instance();  // Broken - AWB 17.03.17

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDEFoam"]         = 0;
   Use["KNN"]             = 0;
   //
   // Linear Discriminant Analysis
   Use["LD"]		  = 0;
   //
   // Function Discriminant analysis
   Use["FDA_GA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   //
   // Neural Network
   Use["MLP"]             = 0;
   Use["DNN"]             = 0;
   //
   // Support Vector Machine
   Use["SVM"]             = 0;
   //
   // Boosted Decision Trees
   Use["BDT"]                     = 0;

   Use["BDTG_default"]            = 1;

   Use["BDTG_AWB"]                = 0;
   Use["BDTG_AWB_lite"]           = 0;

   Use["BDTG_AWB_50_trees"]       = 0;
   Use["BDTG_AWB_100_trees"]      = 0;
   Use["BDTG_AWB_200_trees"]      = 0;
   Use["BDTG_AWB_400_trees"]      = 0;
   Use["BDTG_AWB_800_trees"]      = 0;

   Use["BDTG_AWB_3_deep"]         = 0;
   Use["BDTG_AWB_4_deep"]         = 0;
   Use["BDTG_AWB_5_deep"]         = 0;
   Use["BDTG_AWB_6_deep"]         = 0;

   Use["BDTG_Carnes"]             = 0;

   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification_H2Mu" << std::endl;

   // Select methods (don't look at this code - not of interest)
   std::vector<TString> mlist;
   if (myMethodList != "") {
     // Should be "" for normal use cases - AWB 17.03.17
     std::cout << "\n\nWHY ARE WE USING myMethodList = " << myMethodList << " !!!  Exiting." << std::endl;
     exit(1);
     for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

     // mlist = gTools().SplitString( myMethodList, ',' ); // Broken - AWB 17.03.17
     for (UInt_t i=0; i<mlist.size(); i++) {
       std::string regMethod(mlist[i]);
       
       if (Use.find(regMethod) == Use.end()) {
   	 std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
   	 for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
   	 std::cout << std::endl;
   	 return;
       }
       Use[regMethod] = 1;
     }
   }
   
   // --------------------------------------------------------------------------------------------------

   // Here the preparation phase begins

   // Create a new root output file
   TString out_dir = "/afs/cern.ch/work/a/abrinke1/public/EMTF/HiggsToMuMu";
   out_dir = ".";
   TString out_file_name;
   out_file_name.Form( "%s/TMVAClassification_H2Mu_18_03_13_test.root", out_dir.Data() );
   TFile* out_file = TFile::Open( out_file_name, "RECREATE" );


   ///////////////////////////////////////////////////////
   ///  Input samples: MC signal, MC background, data  ///
   ///////////////////////////////////////////////////////

   std::map<TString, Sample*> samples;
   std::cout << "\nAbout to get the samples" << std::endl;

   // Load samples into our map. "UF" if you are on the UF servers
   // or "CERN" if you at CERN. "ALL" specifies that we want to load the Data
   // and all of the MC samples. Can loop through and remove the ones you don't want 
   // to use if you desire or just grab the ones you care about from the map.
   GetSamples(samples, "CERN", "H2Mu_gg" );
   GetSamples(samples, "CERN", "ZJets_MG_incl");

   std::cout << std::endl << "\nGot the samples" << std::endl;

   // Grab the gluon gluon fusion Higgs to mu mu sample
   Sample* samp_sig = samples["H2Mu_gg"];
   Sample* samp_bkg = samples["ZJets_MG"];

   std::cout << "Using signal sample " << samp_sig->name << std::endl;
   std::cout << "Using background sample " << samp_bkg->name << std::endl;

   // Get branches, set addresses
   // Tells the TTree that it should load the event information into samp->vars
   samp_sig->setBranchAddresses(2);
   samp_bkg->setBranchAddresses(2);

   // Tuple of sample, sample ID, and MC weight
   std::vector< std::tuple<Sample*, int, float> > in_samps;
   in_samps.push_back( std::make_tuple(samp_sig, -1, 1.0) );
   in_samps.push_back( std::make_tuple(samp_bkg, +1, 1.0) );


   //////////////////////////////////////////////////////////////////
   ///  Factories: Use different sets of variables, weights, etc. ///
   //////////////////////////////////////////////////////////////////
   
   TString fact_set = "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification";
   std::vector<TString> var_names; // Holds names of variables for a given factory and permutation
   std::vector<Double_t> var_vals; // Holds values of variables for a given factory and permutation
   TMVA::Factory* nullF = new TMVA::Factory("NULL", out_file, fact_set); // Placeholder factory
   TMVA::DataLoader* nullL = new TMVA::DataLoader("NULL");                 // Placeholder loader

   // Tuple is defined by the factory and dataloader,  followed by a name, 
   // var name and value vectors, and hex bit masks for input variables.
   // Each hex bit represents four variables, e.g. 0x1 would select only the 1st variable, 
   // 0xf the 1st 4, 0xff the 1st 8, 0xa the 2nd and 4th, 0xf1 the 1st and 5th-8th, etc.
   std::vector< std::tuple<TMVA::Factory*, TMVA::DataLoader*, TString, std::vector<TString>, std::vector<Double_t>, int> > factories;

   factories.push_back( std::make_tuple( nullF, nullL, "f_0xffff",
   					 var_names, var_vals, 0xffff) );



   // Initialize factories and dataloaders
   for (UInt_t iFact = 0; iFact < factories.size(); iFact++) {
     std::get<0>(factories.at(iFact)) = new TMVA::Factory( std::get<2>(factories.at(iFact)), out_file, fact_set );
     std::get<1>(factories.at(iFact)) = new TMVA::DataLoader( std::get<2>(factories.at(iFact)) );
   }

   // Defined in interface/MVA_helper.h
   // TMVA_var(TString name, TString descr, TString unit, TString type, Double_t def_val)
   std::vector<TMVA_var> in_vars;   // All input variables
   std::vector<TMVA_var> spec_vars; // All spectator variables
   std::vector<TMVA_var> all_vars;  // All variables
   
   /////////////////////////////////////////////////////////
   ///  Input variables: used in BDT to estimate the pT  ///
   /////////////////////////////////////////////////////////
   
   in_vars.push_back( TMVA_var( "dimu_pt",   "p_{T}(#mu#mu)",       "GeV", 'F', -88 ) ); // 0x0000 0001
   in_vars.push_back( TMVA_var( "dimu_eta",  "#eta(#mu#mu)",           "", 'F', -88 ) ); // 0x0000 0002  
   in_vars.push_back( TMVA_var( "mu1_eta",   "#eta(#mu1)",             "", 'F', -88 ) ); // 0x0000 0004
   in_vars.push_back( TMVA_var( "mu2_eta",   "#eta(#mu2)",             "", 'F', -88 ) ); // 0x0000 0008


   /////////////////////////////////////////////////////////////////////////////
   ///  Spectator variables: not used in training, but saved in output tree  ///
   /////////////////////////////////////////////////////////////////////////////
   
   spec_vars.push_back( TMVA_var( "dimu_mass",    "mass(#mu#mu)",            "GeV", 'F', -77 ) ); // 0x01

   assert( in_vars.size() > 0 );   // You need at least one input variable
   // Order is important: input variables first, then specator
   all_vars.insert( all_vars.end(), in_vars.begin(), in_vars.end() );
   all_vars.insert( all_vars.end(), spec_vars.begin(), spec_vars.end() );


   // Fill each factory with the correct set of variables
   for (UInt_t iFact = 0; iFact < factories.size(); iFact++) {
     std::cout << "\n*** Factory " << std::get<2>(factories.at(iFact)) << " variables ***" << std::endl;
       
     std::cout << "*** Input ***" << std::endl;
     for (UInt_t i = 0; i < in_vars.size(); i++) {
       if ( 0x1 & (std::get<5>(factories.at(iFact)) >> i) ) { // Hex bit mask for in_vars
   	 TMVA_var v = in_vars.at(i);
   	 std::cout << v.name << std::endl;
   	 std::get<1>(factories.at(iFact))->AddVariable( v.name, v.descr, v.unit, v.type ); // Add var to dataloader 
   	 std::get<3>(factories.at(iFact)).push_back( v.name );    // Add to vector of var names
   	 std::get<4>(factories.at(iFact)).push_back( v.def_val ); // Add to vector of var values
       }
     }

     std::cout << "*** Spectator ***" << std::endl;
     for (UInt_t i = 0; i < spec_vars.size(); i++) {
       TMVA_var v = spec_vars.at(i);
       std::cout << v.name << std::endl;
       std::get<1>(factories.at(iFact))->AddSpectator( v.name, v.descr, v.unit, v.type );
       std::get<3>(factories.at(iFact)).push_back( v.name );
       std::get<4>(factories.at(iFact)).push_back( v.def_val );
     }
   } // End loop: for (UInt_t iFact = 0; iFact < factories.size(); iFact++)


   std::cout << "\n******* About to loop over samples *******" << std::endl;
   UInt_t nEvt = 0;
   UInt_t nEvt_sig = 0;
   UInt_t nEvt_bkg = 0;
   UInt_t nTrain_sig = 0;
   UInt_t nTrain_bkg = 0;
   UInt_t nTest_sig  = 0;
   UInt_t nTest_bkg  = 0;

   for (int iSamp = 0; iSamp < in_samps.size(); iSamp++) {
     Sample* samp   = std::get<0>(in_samps.at(iSamp));
     int   samp_ID  = std::get<1>(in_samps.at(iSamp));
     float samp_wgt = std::get<2>(in_samps.at(iSamp));
     
     std::cout << "Looping over the events in the sample" << std::endl;
     // for (UInt_t iEvt = 0; iEvt < samp->GetEntries(); iEvt++) { // TODO: fix. - AWB 18.03.17
     for (UInt_t iEvt = 0; iEvt < MAX_EVT; iEvt++) {
       if (nEvt > MAX_EVT) break;
       if (samp_ID < 0 && nEvt_sig > MAX_SIG) break;
       if (samp_ID > 0 && nEvt_bkg > MAX_BKG) break;
       if ( (nEvt % REPORT_EVT) == 0 )
	 std::cout << "Looking at event " << nEvt << std::endl;

       // Load info from the ttree into samp->vars
       // samp->branches.object (load info) <-> samp->vars.object (access info)
       samp->branches.recoDimuCands->GetEntry(iEvt);
       samp->branches.recoMuons->GetEntry(iEvt);

       if (samp->vars.recoDimuCands->size() == 0)
	 continue;

       MuPairInfo& dimu = samp->vars.recoDimuCands->at(0); 
       samp->vars.dimuCand = &dimu;
       MuonInfo& mu1 = samp->vars.recoMuons->at(dimu.iMu1);
       MuonInfo& mu2 = samp->vars.recoMuons->at(dimu.iMu2);

       std::cout << "Dimuon mass = " << dimu.mass << ", pT = " << dimu.pt << ", eta = " << dimu.eta << std::endl;
       std::cout << "Muon 1 pT = " << mu1.pt << ", eta = " << mu1.eta << std::endl;
       std::cout << "Muon 2 pT = " << mu2.pt << ", eta = " << mu2.eta << std::endl;

       Double_t dimu_mass = dimu.mass;
       Double_t dimu_pt   = dimu.pt;
       Double_t dimu_eta  = dimu.eta;

       Double_t mu1_pt  = mu1.pt;
       Double_t mu1_eta = mu1.eta;
       Double_t mu2_pt  = mu2.pt;
       Double_t mu2_eta = mu2.eta;

       // if ( mu_pt < PTMIN && isMC ) continue;
       // if ( mu_pt > PTMAX && isMC ) continue;
       // if ( fabs( mu_eta ) < ETAMIN && isMC ) continue;
       // if ( fabs( mu_eta ) > ETAMAX && isMC ) continue;
	   
       /////////////////////////////////////////////////////
       ///  Loop over factories and set variable values  ///
       /////////////////////////////////////////////////////
       for (UInt_t iFact = 0; iFact < factories.size(); iFact++) {
	 
	 // Set vars equal to default vector of variables for this factory
	 var_names = std::get<3>(factories.at(iFact));
	 var_vals = std::get<4>(factories.at(iFact));
	 
	 // Fill all variables
	 for (UInt_t iVar = 0; iVar < var_names.size(); iVar++) {
	   TString vName = var_names.at(iVar);
	   
	   /////////////////////////
	   ///  Input variables  ///
	   /////////////////////////
	   
	   // Muon kinematics
	   if ( vName == "dimu_pt" )
	     var_vals.at(iVar) = dimu_pt;
	   if ( vName == "dimu_eta" )
	     var_vals.at(iVar) = dimu_eta;
	   if ( vName == "mu1_eta" )
	     var_vals.at(iVar) = mu1_eta;
	   if ( vName == "mu2_eta" )
	     var_vals.at(iVar) = mu2_eta;


	   /////////////////////////////
	   ///  Spectator variables  ///
	   /////////////////////////////
	   
	   if ( vName == "dimu_mass" )
	     var_vals.at(iVar) = dimu_mass;
	   
	 } // End loop: for (UInt_t iVar = 0; iVar < var_names.size(); iVar++)
	     
	 // Unweighted signal and background
	 Double_t sig_evt_weight = 1.0;
	 Double_t bkg_evt_weight = 1.0;
	     
	 // Load values into event
	 if ( (iEvt % 2) == 0 && nTrain_sig < MAX_TR_SIG && nTrain_bkg < MAX_TR_BKG ) {
	   if (samp_ID < 0) {
	     std::get<1>(factories.at(iFact))->AddSignalTrainingEvent( var_vals, sig_evt_weight );
	     if (iFact == 0) nTrain_sig += 1;
	   }
	   else if (samp_ID > 0) {
	     std::get<1>(factories.at(iFact))->AddBackgroundTrainingEvent( var_vals, bkg_evt_weight );
	     if (iFact == 0) nTrain_bkg += 1;
	   }
	 }
	 else {
	   if (samp_ID < 0) {
	     std::get<1>(factories.at(iFact))->AddSignalTestEvent( var_vals, sig_evt_weight );
	     if (iFact == 0) nTest_sig += 1;
	   } 
	   else if (samp_ID > 0) {
	     std::get<1>(factories.at(iFact))->AddBackgroundTestEvent( var_vals, bkg_evt_weight );
	     if (iFact == 0) nTest_bkg += 1;
	   }
	 }

       } // End loop: for (UInt_t iFact = 0; iFact < factories.size(); iFact++) 
       
       nEvt += 1;
       if (samp_ID < 0) nEvt_sig += 1;
       if (samp_ID > 0) nEvt_bkg += 1;
     } // End loop: for (UInt_t iEvt = 0; iEvt < samp->GetEntries(); iEvt++)
   } // End loop: for (int iSamp = 0; iSamp < in_samps.size(); iSamp++)


   std::cout << "******* Made it out of the event loop *******" << std::endl;
   
   std::string NTrS;
   std::string NTrB;

   std::ostringstream convertTrS;
   convertTrS << nTrain_sig;
   NTrS = convertTrS.str();
   std::ostringstream convertTrB;
   convertTrB << nTrain_bkg;
   NTrB = convertTrB.str();

   std::string numTrainStr = "nTrain_Signal="+NTrS+":nTrain_Background="+NTrB+":";
   std::cout << "NTrS: " << NTrS << ", NTrB: " << NTrB << std::endl;

   // // global event weights per tree (see below for setting event-wise weights)
   // Double_t regWeight  = 1.0;

   for (UInt_t iFact = 0; iFact < factories.size(); iFact++) {
     
     TMVA::Factory* factX = std::get<0>(factories.at(iFact));
     TMVA::DataLoader* loadX = std::get<1>(factories.at(iFact));
     
     // // You can add an arbitrary number of regression trees
     // loadX->AddRegressionTree( regTree, regWeight );
     
     // // This would set individual event weights (the variables defined in the
     // // expression need to exist in the original TTree)
     // loadX->SetWeightExpression( "var1", "Regression" );
     loadX->SetWeightExpression( 1.0 );
     
     // // Apply additional cuts on the signal and background samples (can be different)
     // TCut mycut = "( abs(muon.eta[0]) > 1.25 && abs(muon.eta[1]) < 2.4 )"; // && track.mode[0] == 15 )"; 
     
     // Tell the dataloader how to use the training and testing events
     loadX->PrepareTrainingAndTestTree( "", "", numTrainStr+"SplitMode=Random:NormMode=NumEvents:!V" );   
     // loadX->PrepareTrainingAndTestTree( mycut, "nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V" );
     
     // If no numbers of events are given, half of the events in the tree are used
     // for training, and the other half for testing:
     //
     //     loadX->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
     
     // Book MVA methods
     //
     // Please lookup the various method configuration options in the corresponding cxx files, eg:
     // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
     // it is possible to preset ranges in the option string in which the cut optimisation should be done:
     // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable
     
     // Linear discriminant
     if (Use["LD"])
       factX->BookMethod( loadX,  TMVA::Types::kLD, "LD",
   			  "!H:!V:VarTransform=None" );
     
     // Neural network (MLP)
     if (Use["MLP"])
       factX->BookMethod( loadX,  TMVA::Types::kMLP, "MLP", (std::string)
   			  "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:"+
   			  "TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:"+
   			  "ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );
     
     if (Use["DNN"])
       {
	 
   	 // TString layoutString ("Layout=TANH|(N+100)*2,LINEAR");
   	 // TString layoutString ("Layout=SOFTSIGN|100,SOFTSIGN|50,SOFTSIGN|20,LINEAR");
   	 // TString layoutString ("Layout=RELU|300,RELU|100,RELU|30,RELU|10,LINEAR");
   	 // TString layoutString ("Layout=SOFTSIGN|50,SOFTSIGN|30,SOFTSIGN|20,SOFTSIGN|10,LINEAR");
   	 // TString layoutString ("Layout=TANH|50,TANH|30,TANH|20,TANH|10,LINEAR");
   	 // TString layoutString ("Layout=SOFTSIGN|50,SOFTSIGN|20,LINEAR");
   	 // TString layoutString ("Layout=TANH|100,TANH|30,LINEAR");

   	 TString layoutString ("Layout=TANH|100,LINEAR");
	 
   	 TString training0 ( (std::string) "LearningRate=1e-5,Momentum=0.5,Repetitions=1,"+
   			     "ConvergenceSteps=500,BatchSize=50,TestRepetitions=7,WeightDecay=0.01,"+
   			     "Regularization=NONE,DropConfig=0.5+0.5+0.5+0.5,DropRepetitions=2");
   	 TString training1 ( (std::string) "LearningRate=1e-5,Momentum=0.9,Repetitions=1,"+
   			     "ConvergenceSteps=170,BatchSize=30,TestRepetitions=7,WeightDecay=0.01,"+
   			     "Regularization=L2,DropConfig=0.1+0.1+0.1,DropRepetitions=1");
   	 TString training2 ( (std::string) "LearningRate=1e-5,Momentum=0.3,Repetitions=1,ConvergenceSteps=150,"+
   			     "BatchSize=40,TestRepetitions=7,WeightDecay=0.01,Regularization=NONE");
   	 TString training3 ( (std::string) "LearningRate=1e-6,Momentum=0.1,Repetitions=1,ConvergenceSteps=500,"+
   			     "BatchSize=100,TestRepetitions=7,WeightDecay=0.0001,Regularization=NONE");
	 
   	 TString trainingStrategyString ("TrainingStrategy=");
   	 trainingStrategyString += training0 + "|" + training1 + "|" + training2 + "|" + training3;
	 
	 
   	 // TString trainingStrategyString ( (std::string) "TrainingStrategy=LearningRate=1e-1,Momentum=0.3,"+
   	 // 				  "Repetitions=3,ConvergenceSteps=20,BatchSize=30,TestRepetitions=7,"+
   	 // 				  "WeightDecay=0.0,L1=false,DropFraction=0.0,DropRepetitions=5");
	 
   	 TString nnOptions ("!H:V:ErrorStrategy=SUMOFSQUARES:VarTransform=G:WeightInitialization=XAVIERUNIFORM");
   	 // TString nnOptions ("!H:V:VarTransform=Normalize:ErrorStrategy=CHECKGRADIENTS");
   	 nnOptions.Append (":"); nnOptions.Append (layoutString);
   	 nnOptions.Append (":"); nnOptions.Append (trainingStrategyString);
	 
   	 factX->BookMethod(loadX, TMVA::Types::kDNN, "DNN", nnOptions ); // NN
       }


     // Support Vector Machine
     if (Use["SVM"])
       factX->BookMethod( loadX,  TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );
     
     // Boosted Decision Trees
     if (Use["BDT"])
       factX->BookMethod( loadX,  TMVA::Types::kBDT, "BDT", (std::string)
   			  "!H:!V:NTrees=100:MinNodeSize=1.0%:BoostType=AdaBoostR2:"+
   			  "nCuts=20:PruneMethod=CostComplexity:PruneStrength=30" );
     
     // Default TMVA settings
     if (Use["BDTG_default"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_default", (std::string)
   			  "!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
   			  "BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );

     // AWB settings - AbsoluteDeviation
     if (Use["BDTG_AWB"]) // Optimized settings
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB", (std::string)
   			  "!H:!V:NTrees=400::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.000001" );

     if (Use["BDTG_AWB_lite"]) // Fast, simple BDT
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_lite", (std::string)
   			  "!H:!V:NTrees=40::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=3:MinNodeSize=0.000001" );

     if (Use["BDTG_AWB_50_trees"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_50_trees", (std::string)
   			  "!H:!V:NTrees=50::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=3:MinNodeSize=0.001" );
     if (Use["BDTG_AWB_100_trees"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_100_trees", (std::string)
   			  "!H:!V:NTrees=100::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=3:MinNodeSize=0.001" );
     if (Use["BDTG_AWB_200_trees"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_200_trees", (std::string)
   			  "!H:!V:NTrees=200::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=3:MinNodeSize=0.001" );
     if (Use["BDTG_AWB_400_trees"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_400_trees", (std::string)
   			  "!H:!V:NTrees=400::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=3:MinNodeSize=0.001" );
     if (Use["BDTG_AWB_800_trees"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_800_trees", (std::string)
   			  "!H:!V:NTrees=800::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=3:MinNodeSize=0.001" );
     
     if (Use["BDTG_AWB_3_deep"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_3_deep", (std::string)
   			  "!H:!V:NTrees=100::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=4:MinNodeSize=0.001" );
     if (Use["BDTG_AWB_4_deep"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_4_deep", (std::string)
   			  "!H:!V:NTrees=100::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=4:MinNodeSize=0.001" );
     if (Use["BDTG_AWB_5_deep"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_5_deep", (std::string)
   			  "!H:!V:NTrees=100::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.001" );
     if (Use["BDTG_AWB_6_deep"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_6_deep", (std::string)
   			  "!H:!V:NTrees=100::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=6:MinNodeSize=0.001" );
     
     // Factory settings from Andrew Carnes ... what do they do? - AWB 04.01.17
     if (Use["BDTG_Carnes"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_Carnes", (std::string)
   			  "!H:!V:NTrees=64::BoostType=Grad:Shrinkage=0.3:nCuts=99999:MaxDepth=4:MinNodeSize=0.001:"+
   			  "NegWeightTreatment=IgnoreNegWeightsInTraining:PruneMethod=NoPruning" );
     
     
     // --------------------------------------------------------------------------------------------------
     
     // Now you can tell the factory to train, test, and evaluate the MVAs
     
     // Train MVAs using the set of training events
     factX->TrainAllMethods();
     
     // Evaluate all MVAs using the set of test events
     factX->TestAllMethods();
     
     // Evaluate and compare performance of all configured MVAs
     factX->EvaluateAllMethods();

     // // Instead of "EvaluateAllMethods()", just write out the training and testing trees
     // // Skip unnecessary evaluatioh histograms, which take time on large datasets 
     // // Code gleaned from original "EvaluateAllMethods()" function in tmva/tmva/src/Factory.cxx - AWB 31.01.17
     // if ( factX->fMethodsMap.empty() )
     //   std::cout << "factX->fMethodsMap is empty" << std::endl;
     
     // std::map<TString, std::vector<IMethod*>*>::iterator itrMap;
     // for (itrMap = factX->fMethodsMap.begin(); itrMap != factX->fMethodsMap.end(); itrMap++) {
       
     //   std::vector<IMethod*> *methods = itrMap->second;
     //   std::list<TString> datasets;
     //   Int_t nmeth_used[2] = {int(mlist.size()), 1};
       
     //   for (Int_t k = 0; k < 2; k++) {
     // 	 for (Int_t i = 0; i < nmeth_used[k]; i++) {
     // 	   MethodBase* theMethod = dynamic_cast<MethodBase*>((*methods)[i]);
     // 	   if (theMethod == 0) {
     // 	     std::cout << "For k = " << k << ", i = " << i << ", no valid method" << std::endl;
     // 	     continue;
     // 	   }
     // 	   TDirectory* RootBaseDir = (TDirectory*) out_file;
     // 	   RootBaseDir->cd( std::get<2>(factories.at(iFact)) );
     // 	   if ( std::find( datasets.begin(), datasets.end(), std::get<2>(factories.at(iFact)) ) == datasets.end() ) {
     // 	     theMethod->Data()->GetTree(Types::kTesting)->Write( "", TObject::kOverwrite );
     // 	     theMethod->Data()->GetTree(Types::kTraining)->Write( "", TObject::kOverwrite );
     // 	     datasets.push_back( std::get<2>(factories.at(iFact)) );
     // 	   }
     // 	 } // End loop: for (Int_t i = 0; i < nmeth_used[k]; i++)
     //   } // End loop: for (Int_t k = 0; k < 2; k++)
     // } // End loop: for (itrMap = factX->fMethodsMap.begin(); itrMap != factX->fMethodsMap.end(); itrMap++) 

     // --------------------------------------------------------------
     
   } // End loop: for (UInt_t iFact = 0; iFact < factories.size(); iFact++)
   
   // Save the output
   out_file->Close();

   std::cout << "==> Wrote root file: " << out_file->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   // delete factory;
   // delete dataloader;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVAGui( out_file_name );
}


int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   TMVAClassification_H2Mu(methodList);
   return 0;
}

