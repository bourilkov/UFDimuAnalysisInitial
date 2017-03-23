
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
#include "MuonSelection.h"
#include "EventSelection.h"
#include "CategorySelection.h"
#include "MuonCollectionCleaner.h"
#include "JetCollectionCleaner.h"
#include "EleCollectionCleaner.h"

// Extra tools - AWB 13.03.17
#include "TMVA_helper.h"

const int MAX_EVT    = 1000000;
const int MAX_SIG    =  500000;
const int MAX_BKG    =  500000;
const int REPORT_EVT =    1000;

const int MAX_TR_SIG = 10000000;
const int MAX_TR_BKG = 10000000;

const double PI = 3.14159265359;
const double BIT = 0.000001; // Tiny value or offset

const double LUMI = 36814; // pb-1 

using namespace TMVA;

void TMVAClassification_H2Mu ( TString myMethodList = "" ) {

   TMVA::Tools::Instance();

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
     for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
     mlist = gTools().SplitString( myMethodList, ',' );
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
   TString out_dir = "/afs/cern.ch/work/a/abrinke1/public/H2Mu/TMVA";
   // out_dir = ".";
   TString out_file_name;
   out_file_name.Form( "%s/TMVAClassification_H2Mu_17_03_23_test.root", out_dir.Data() );
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
   GetSamples(samples, "CERN_hiM", "SIGNAL" );
   GetSamples(samples, "CERN_hiM", "ZJets_MG");
   GetSamples(samples, "CERN_hiM", "tt_ll_MG");

   std::cout << std::endl << "\nGot the samples" << std::endl;

   // Tuple of sample and sample ID
   std::vector< std::tuple<Sample*, int> > sig_samps;
   std::vector< std::tuple<Sample*, int> > bkg_samps;
   std::vector< std::tuple<Sample*, int> > dat_samps;
   std::vector< std::tuple<Sample*, int> > all_samps;

   sig_samps.push_back( std::make_tuple(samples["H2Mu_gg"],     -1) );
   sig_samps.push_back( std::make_tuple(samples["H2Mu_VBF"],    -2) );
   sig_samps.push_back( std::make_tuple(samples["H2Mu_ZH"],     -3) );
   sig_samps.push_back( std::make_tuple(samples["H2Mu_WH_pos"], -4) );
   sig_samps.push_back( std::make_tuple(samples["H2Mu_WH_neg"], -5) );

   bkg_samps.push_back( std::make_tuple(samples["ZJets_MG"],              + 1) );
   bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_70_100"],    + 2) );
   bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_100_200"],   + 3) );
   bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_200_400"],   + 4) );
   bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_400_600"],   + 5) );
   bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_600_800"],   + 6) );
   bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_800_1200"],  + 7) );
   bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_1200_2500"], + 8) );
   bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_2500_inf"],  + 9) );
   bkg_samps.push_back( std::make_tuple(samples["tt_ll_MG"],              +10) );

   all_samps.insert( all_samps.end(), sig_samps.begin(), sig_samps.end() );
   all_samps.insert( all_samps.end(), bkg_samps.begin(), bkg_samps.end() );
   all_samps.insert( all_samps.end(), dat_samps.begin(), dat_samps.end() );

   // Get branches, set addresses
   // Tells the TTree that it should load the event information into samp->vars
   for (int iSamp = 0; iSamp < all_samps.size(); iSamp++) {
     std::get<0>(all_samps.at(iSamp))->setBranchAddresses(2);
     std::get<0>(all_samps.at(iSamp))->calculateNoriginal();
   }

   // Objects to help clean valid objects from the net collections in samp->vars
   JetCollectionCleaner   jetCollectionCleaner;
   MuonCollectionCleaner  muonCollectionCleaner;
   EleCollectionCleaner   eleCollectionCleaner;

   // Objects to cut events from the analysis
   // See selection/MuonSelection.h or selection/EventSelection.h
   // Choose from an available option or make your own implementing the interface
   // Use selection.evaluate(samp->vars) to see whether an event passes or fails the cuts
   Run2MuonSelectionCuts   run2MuonSelection;
   Run2EventSelectionCuts  run2EventSelectionMC;

   // Object to categorize the event appropriately.  Can make your own categorizer or use one of those in
   //   selection/CategorySelection.h.  Use one defined there or make your own implementing the interface.
   
   // Categorizer object has a map of the different categories in Categorizer.categoryMap<TString, Category>
   //   which maps the category name to the Category object
   // Then each category object can store histograms in category.histoMap<TString, TH1D*>
   // Use categorySelection.evaluate(s->vars) to see which category the event falls into
   CategorySelectionRun1 categorySelection;


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
   std::vector< std::tuple<TMVA::Factory*, TMVA::DataLoader*, TString, std::vector<TString>, std::vector<Double_t>, int, int, int> > factories;

   // factories.push_back( std::make_tuple( nullF, nullL, "f_muVars", var_names, var_vals, 
   // 					 0xffff, 0x0000, 0x0000) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_jetVars", var_names, var_vals, 
   // 					 0x0000, 0xffff, 0x0000) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_evtVars", var_names, var_vals, 
   // 					 0x0000, 0x0000, 0xfffe) );

   factories.push_back( std::make_tuple( nullF, nullL, "f_BASE", var_names, var_vals, 
					 0x001c, 0x0ff0, 0x0011) );
   factories.push_back( std::make_tuple( nullF, nullL, "f_Opt1", var_names, var_vals, 
					 0x045c, 0x0fff, 0x001e) );


   // Initialize factories and dataloaders
   for (UInt_t iFact = 0; iFact < factories.size(); iFact++) {
     std::get<0>(factories.at(iFact)) = new TMVA::Factory( std::get<2>(factories.at(iFact)), out_file, fact_set );
     std::get<1>(factories.at(iFact)) = new TMVA::DataLoader( std::get<2>(factories.at(iFact)) );
   }

   // Defined in interface/MVA_helper.h
   // TMVA_var(TString name, TString descr, TString unit, TString type, Double_t def_val)
   std::vector<TMVA_var> mu_vars;    // Muon input variables
   std::vector<TMVA_var> jet_vars;   // Jet input variables
   std::vector<TMVA_var> evt_vars;   // Global event / combined object input variables
   // std::vector<TMVA_var> in_vars;    // All input variables
   std::vector<TMVA_var> spec_vars;  // All spectator variables
   // std::vector<TMVA_var> all_vars;   // All variables
   
   /////////////////////////////////////////////////////////
   ///  Input variables: used in BDT to estimate the pT  ///
   /////////////////////////////////////////////////////////

   // Muon variables
   mu_vars.push_back( TMVA_var( "mu1_pt",     "p_{T}(#mu1)",       "GeV", 'F', -88 ) ); // 0x0001
   mu_vars.push_back( TMVA_var( "mu2_pt",     "p_{T}(#mu2)",       "GeV", 'F', -88 ) ); // 0x0002
   mu_vars.push_back( TMVA_var( "mu1_eta",    "#eta(#mu1)",           "", 'F', -88 ) ); // 0x0004
   mu_vars.push_back( TMVA_var( "mu2_eta",    "#eta(#mu2)",           "", 'F', -88 ) ); // 0x0008

   mu_vars.push_back( TMVA_var( "dimu_pt",    "p_{T}(#mu#mu)",     "GeV", 'F', -88 ) ); // 0x0010
   mu_vars.push_back( TMVA_var( "dimu_dMass", "#sigma M(#mu#mu)",  "GeV", 'F', -88 ) ); // 0x0020
   mu_vars.push_back( TMVA_var( "dimu_eta",   "#eta(#mu#mu)",         "", 'F', -88 ) ); // 0x0040  
   mu_vars.push_back( TMVA_var( "dimu_rapid", "rapid(#mu#mu)",        "", 'F', -88 ) ); // 0x0080 

   mu_vars.push_back( TMVA_var( "dimu_dR",    "dR(#mu#mu)",           "", 'F', -88 ) ); // 0x0100
   mu_vars.push_back( TMVA_var( "dimu_dEta",  "d#eta(#mu#mu)",        "", 'F', -88 ) ); // 0x0200
   mu_vars.push_back( TMVA_var( "dimu_dPhi",  "d#phi(#mu#mu)",        "", 'F', -88 ) ); // 0x0400  
   mu_vars.push_back( TMVA_var( "dimu_dPhiS", "d#phi*(#mu#mu)",       "", 'F', -88 ) ); // 0x0800 

   // Jet variables
   jet_vars.push_back( TMVA_var( "jet1_pt",     "p_{T}(jet1)",       "GeV", 'F', -88 ) ); // 0x0001
   jet_vars.push_back( TMVA_var( "jet2_pt",     "p_{T}(jet2)",       "GeV", 'F', -88 ) ); // 0x0002
   jet_vars.push_back( TMVA_var( "jet1_eta",    "#eta(jet1)",           "", 'F', -88 ) ); // 0x0004
   jet_vars.push_back( TMVA_var( "jet2_eta",    "#eta(jet2)",           "", 'F', -88 ) ); // 0x0008

   jet_vars.push_back( TMVA_var( "dijet1_mass", "1^{st} M(jj)",      "GeV", 'F', -88 ) ); // 0x0010
   jet_vars.push_back( TMVA_var( "dijet2_mass", "2^{nd} M(jj)",      "GeV", 'F', -88 ) ); // 0x0020
   jet_vars.push_back( TMVA_var( "dijet3_mass", "3^{rd} M(jj)",      "GeV", 'F', -88 ) ); // 0x0040
   jet_vars.push_back( TMVA_var( "dijet4_mass", "4^{th} M(jj)",      "GeV", 'F', -88 ) ); // 0x0080

   jet_vars.push_back( TMVA_var( "dijet1_dEta", "1^{st} d#eta(jj)",     "", 'F', -88 ) ); // 0x0100
   jet_vars.push_back( TMVA_var( "dijet2_dEta", "2^{nd} d#eta(jj)",     "", 'F', -88 ) ); // 0x0200
   jet_vars.push_back( TMVA_var( "dijet3_dEta", "3^{rd} d#eta(jj)",     "", 'F', -88 ) ); // 0x0400
   jet_vars.push_back( TMVA_var( "dijet4_dEta", "4^{th} d#eta(jj)",     "", 'F', -88 ) ); // 0x0800

   // Global event variables
   evt_vars.push_back( TMVA_var( "nJets",       "# of jets",            "", 'I', -88 ) ); // 0x0001
   evt_vars.push_back( TMVA_var( "nJetsCent",   "# of central jets",    "", 'I', -88 ) ); // 0x0002
   evt_vars.push_back( TMVA_var( "nJetsFwd",    "# of forward jets",    "", 'I', -88 ) ); // 0x0004
   evt_vars.push_back( TMVA_var( "nBMed",       "# of medium b-tags",   "", 'I', -88 ) ); // 0x0008

   evt_vars.push_back( TMVA_var( "MET",         "MET",               "GeV", 'F', -88 ) ); // 0x0010
   evt_vars.push_back( TMVA_var( "MHT",         "MHT",               "GeV", 'F', -88 ) ); // 0x0020
   evt_vars.push_back( TMVA_var( "MT_had",      "M_{T} of jets",     "GeV", 'F', -88 ) ); // 0x0040
   evt_vars.push_back( TMVA_var( "mass_had",    "Mass of jets",      "GeV", 'F', -88 ) ); // 0x0080


   /////////////////////////////////////////////////////////////////////////////
   ///  Spectator variables: not used in training, but saved in output tree  ///
   /////////////////////////////////////////////////////////////////////////////
   
   spec_vars.push_back( TMVA_var( "samp_ID",    "Sample ID",            "", 'I', -77 ) );
   spec_vars.push_back( TMVA_var( "samp_wgt",   "Sample weight",        "", 'F', -77 ) );
   spec_vars.push_back( TMVA_var( "LHE_HT",     "Sample weight",     "GeV", 'F', -77 ) );
   spec_vars.push_back( TMVA_var( "dimu_mass",  "mass(#mu#mu)",      "GeV", 'F', -77 ) );
   spec_vars.push_back( TMVA_var( "BASE_cat",   "BASELINE category",    "", 'I', -77 ) );


   // Fill each factory with the correct set of variables
   for (UInt_t iFact = 0; iFact < factories.size(); iFact++) {
     std::cout << "\n*** Factory " << std::get<2>(factories.at(iFact)) << " variables ***" << std::endl;
       
     std::cout << "*** Input muon variables ***" << std::endl;
     for (UInt_t i = 0; i < mu_vars.size(); i++) {
       if ( 0x1 & (std::get<5>(factories.at(iFact)) >> i) ) { // Hex bit mask for mu_vars
   	 TMVA_var v = mu_vars.at(i);
   	 std::cout << v.name << std::endl;
   	 std::get<1>(factories.at(iFact))->AddVariable( v.name, v.descr, v.unit, v.type ); // Add var to dataloader 
   	 std::get<3>(factories.at(iFact)).push_back( v.name );    // Add to vector of var names
   	 std::get<4>(factories.at(iFact)).push_back( v.def_val ); // Add to vector of var values
       }
     }

     std::cout << "*** Input jet variables ***" << std::endl;
     for (UInt_t i = 0; i < jet_vars.size(); i++) {
       if ( 0x1 & (std::get<6>(factories.at(iFact)) >> i) ) { // Hex bit mask for jet_vars
   	 TMVA_var v = jet_vars.at(i);
   	 std::cout << v.name << std::endl;
   	 std::get<1>(factories.at(iFact))->AddVariable( v.name, v.descr, v.unit, v.type ); // Add var to dataloader 
   	 std::get<3>(factories.at(iFact)).push_back( v.name );    // Add to vector of var names
   	 std::get<4>(factories.at(iFact)).push_back( v.def_val ); // Add to vector of var values
       }
     }

     std::cout << "*** Input event variables ***" << std::endl;
     for (UInt_t i = 0; i < evt_vars.size(); i++) {
       if ( 0x1 & (std::get<7>(factories.at(iFact)) >> i) ) { // Hex bit mask for evt_vars
   	 TMVA_var v = evt_vars.at(i);
   	 std::cout << v.name << std::endl;
   	 std::get<1>(factories.at(iFact))->AddVariable( v.name, v.descr, v.unit, v.type ); // Add var to dataloader 
   	 std::get<3>(factories.at(iFact)).push_back( v.name );    // Add to vector of var names
   	 std::get<4>(factories.at(iFact)).push_back( v.def_val ); // Add to vector of var values
       }
     }

     std::cout << "*** Spectator variables ***" << std::endl;
     for (UInt_t i = 0; i < spec_vars.size(); i++) {
       TMVA_var v = spec_vars.at(i);
       std::cout << v.name << std::endl;
       std::get<1>(factories.at(iFact))->AddSpectator( v.name, v.descr, v.unit, v.type );
       std::get<3>(factories.at(iFact)).push_back( v.name );
       std::get<4>(factories.at(iFact)).push_back( v.def_val );
     }
   } // End loop: for (UInt_t iFact = 0; iFact < factories.size(); iFact++)

   // in_vars.insert( in_vars.end(), mu_vars.begin(),  mu_vars.end()  );
   // in_vars.insert( in_vars.end(), jet_vars.begin(), jet_vars.end() );
   // in_vars.insert( in_vars.end(), evt_vars.begin(), evt_vars.end() );
   // assert( in_vars.size() > 0 );   // You need at least one input variable
   // // Order is important: input variables first, then specator
   // all_vars.insert( all_vars.end(), in_vars.begin(), in_vars.end() );
   // all_vars.insert( all_vars.end(), spec_vars.begin(), spec_vars.end() );


   std::cout << "\n******* About to loop over samples *******" << std::endl;
   UInt_t nEvt = 0;
   UInt_t nEvt_sig = 0;
   UInt_t nEvt_bkg = 0;
   UInt_t nTrain_sig = 0;
   UInt_t nTrain_bkg = 0;
   UInt_t nTest_sig  = 0;
   UInt_t nTest_bkg  = 0;

   for (int iSamp = 0; iSamp < all_samps.size(); iSamp++) {
     Sample* samp  = std::get<0>(all_samps.at(iSamp));
     int   samp_ID = std::get<1>(all_samps.at(iSamp));
     
     std::cout << "Looping over the " << samp->N << " events in sample " << samp->name << std::endl;
     for (UInt_t iEvt = 0; iEvt < samp->N; iEvt++) {
       if (nEvt > MAX_EVT) break;
       if ( samp_ID < 0 && nEvt_sig > ((MAX_SIG * (iSamp + 1))                    / sig_samps.size()) ) break;
       if ( samp_ID > 0 && nEvt_bkg > ((MAX_BKG * (iSamp + 1 - sig_samps.size())) / bkg_samps.size()) ) break;
       if ( (nEvt % REPORT_EVT) == 0 )
	 std::cout << "Looking at event " << nEvt << std::endl;

       
       ////////////////////////////////////////////////////////
       ///  Begin block mostly lifted from bin/example.cxx  ///
       ////////////////////////////////////////////////////////

       // Load info from the ttree into samp->vars
       // samp->branches.object (load info) <-> samp->vars.object (access info)
       samp->branches.muPairs->GetEntry(iEvt);

       if (samp->vars.muPairs->size() == 0)
	 continue;
       if ( samp->vars.muPairs->at(0).mass < 120 ||
	    samp->vars.muPairs->at(0).mass > 130 )
	 continue;

       // getEntry from lib/BranchSet.h allows us to set all of the branches at once
       samp->branches.getEntry(iEvt);

       // Reset the categories so we get the correct categorization for this event
       categorySelection.reset();

       // Set the Higgs dimuon candidate pair
       samp->vars.dimuCand = &(samp->vars.muPairs->at(0));
       
       // Clear vectors for the valid collections
       samp->vars.validMuons.clear();
       samp->vars.validExtraMuons.clear();
       samp->vars.validElectrons.clear();
       samp->vars.validJets.clear();
       samp->vars.validBJets.clear();

       // Load valid collections from samp->vars raw collections
       muonCollectionCleaner.getValidMuons(samp->vars, samp->vars.validMuons, samp->vars.validExtraMuons);
       jetCollectionCleaner.getValidJets(samp->vars, samp->vars.validJets, samp->vars.validBJets);
       eleCollectionCleaner.getValidElectrons(samp->vars, samp->vars.validElectrons);

       // Clean jets from muons
       CollectionCleaner::cleanByDR(samp->vars.validJets, samp->vars.validMuons, 0.4);

       // Can now categorize the event with our categorizer
       categorySelection.evaluate(samp->vars);
       // // See which categories the event fell into
       // std::cout << std::endl << "\nCategorizing event " << iEvt << " ..." << std::endl;
       // categorySelection.outputResults();

       // Look at each category
       Int_t BASE_cat = 0;
       for (auto &cat : categorySelection.categoryMap) {
	 // cat.first is the category name, cat.second is the category object
	 if (!cat.second.inCategory) continue;
	 if (cat.first == "c_01_Jet_Loose_EE") BASE_cat =  1;
	 if (cat.first == "c_01_Jet_Loose_OE") BASE_cat =  2;
	 if (cat.first == "c_01_Jet_Loose_BE") BASE_cat =  3;
	 if (cat.first == "c_01_Jet_Loose_OO") BASE_cat =  4;
	 if (cat.first == "c_01_Jet_Loose_BO") BASE_cat =  5;
	 if (cat.first == "c_01_Jet_Loose_BB") BASE_cat =  6;
	 if (cat.first == "c_01_Jet_Tight_EE") BASE_cat =  7;
	 if (cat.first == "c_01_Jet_Tight_OE") BASE_cat =  8;
	 if (cat.first == "c_01_Jet_Tight_BE") BASE_cat =  9;
	 if (cat.first == "c_01_Jet_Tight_OO") BASE_cat = 10;
	 if (cat.first == "c_01_Jet_Tight_BO") BASE_cat = 11;
	 if (cat.first == "c_01_Jet_Tight_BB") BASE_cat = 12;
	 if (cat.first == "c_2_Jet_GGF_Tight") BASE_cat = 13;
	 if (cat.first == "c_2_Jet_VBF_Loose") BASE_cat = 14;
	 if (cat.first == "c_2_Jet_VBF_Tight") BASE_cat = 15;
	 // std::cout << "In category " << cat.first << std::endl;
       }
       if (BASE_cat == 0) continue;

       //////////////////////////////////////////////////////
       ///  End block mostly lifted from bin/example.cxx  ///
       //////////////////////////////////////////////////////

       MuPairInfo& dimu = samp->vars.muPairs->at(0); 
       MuonInfo& mu1 = samp->vars.muons->at(dimu.iMu1);
       MuonInfo& mu2 = samp->vars.muons->at(dimu.iMu2);

       Double_t samp_wgt = samp->getWeight() * samp->getLumiScaleFactor(LUMI);
       if (samp_ID < 0 && 1.0 * MAX_SIG / sig_samps.size() < samp->N) 
	 samp_wgt *= (1.0 * samp->N * sig_samps.size() / MAX_SIG);
       if (samp_ID > 0 && 1.0 * MAX_BKG / bkg_samps.size() < samp->N) 
	 samp_wgt *= (1.0 * samp->N * bkg_samps.size() / MAX_BKG);
       Double_t LHE_HT   = samp->vars.lhe_ht;

       // std::cout << "Event " << iEvt << ", sample " << samp->name << " has weight " << samp->getWeight() 
       // 		 << ", lumi scale factor " << samp->getLumiScaleFactor(LUMI) << " (total = " << samp_wgt
       // 		 << "), and LHE HT = " << LHE_HT << std::endl;
       
       // std::cout << "Dimuon mass = " << dimu.mass << ", pT = " << dimu.pt << ", eta = " << dimu.eta << std::endl;
       // std::cout << "Muon 1 pT = " << mu1.pt << ", eta = " << mu1.eta << std::endl;
       // std::cout << "Muon 2 pT = " << mu2.pt << ", eta = " << mu2.eta << std::endl;

       Double_t mu1_pt  = mu1.pt;
       Double_t mu1_eta = mu1.eta;
       Double_t mu2_pt  = mu2.pt;
       Double_t mu2_eta = mu2.eta;

       Double_t dimu_pt    = dimu.pt;
       Double_t dimu_mass  = dimu.mass;
       Double_t dimu_dMass = dimu.massErr;
       Double_t dimu_eta   = dimu.eta;
       Double_t dimu_rapid = dimu.rapid;

       Double_t dimu_dR    = dimu.dR;
       Double_t dimu_dEta  = dimu.dEta;
       Double_t dimu_dPhi  = dimu.dPhi;
       Double_t dimu_dPhiS = dimu.dPhiStar;

       Double_t jet1_pt  = (samp->vars.jets->size() > 0) ? samp->vars.jets->at(0).pt  : 0;
       Double_t jet2_pt  = (samp->vars.jets->size() > 1) ? samp->vars.jets->at(1).pt  : 0;
       Double_t jet1_eta = (samp->vars.jets->size() > 0) ? samp->vars.jets->at(0).eta : 0;
       Double_t jet2_eta = (samp->vars.jets->size() > 1) ? samp->vars.jets->at(1).eta : 0;

       Double_t dijet1_mass = (samp->vars.jetPairs->size() > 0) ? samp->vars.jetPairs->at(0).mass : 0;
       Double_t dijet2_mass = (samp->vars.jetPairs->size() > 1) ? samp->vars.jetPairs->at(1).mass : 0;
       Double_t dijet3_mass = (samp->vars.jetPairs->size() > 2) ? samp->vars.jetPairs->at(2).mass : 0;
       Double_t dijet4_mass = (samp->vars.jetPairs->size() > 3) ? samp->vars.jetPairs->at(3).mass : 0;

       Double_t dijet1_dEta = (samp->vars.jetPairs->size() > 0) ? samp->vars.jetPairs->at(0).dEta : 0;
       Double_t dijet2_dEta = (samp->vars.jetPairs->size() > 1) ? samp->vars.jetPairs->at(1).dEta : 0;
       Double_t dijet3_dEta = (samp->vars.jetPairs->size() > 2) ? samp->vars.jetPairs->at(2).dEta : 0;
       Double_t dijet4_dEta = (samp->vars.jetPairs->size() > 3) ? samp->vars.jetPairs->at(3).dEta : 0;

       Double_t nJets     = samp->vars.nJets;
       Double_t nJetsCent = samp->vars.nJetsCent;
       Double_t nJetsFwd  = samp->vars.nJetsFwd;
       Double_t nBMed     = samp->vars.nBMed;

       Double_t MET      = samp->vars.met->pt;
       Double_t MHT      = samp->vars.mht->pt;
       Double_t MT_had   = samp->vars.mht->MT_had;
       Double_t mass_had = samp->vars.mht->mass_had;

	   
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
	   
	   // Muon variables
	   if ( vName == "mu1_pt" )
	     var_vals.at(iVar) = mu1_pt;
	   if ( vName == "mu2_pt" )
	     var_vals.at(iVar) = mu2_pt;
	   if ( vName == "mu1_eta" )
	     var_vals.at(iVar) = mu1_eta;
	   if ( vName == "mu2_eta" )
	     var_vals.at(iVar) = mu2_eta;

	   if ( vName == "dimu_pt" )
	     var_vals.at(iVar) = dimu_pt;
	   if ( vName == "dimu_dMass" )
	     var_vals.at(iVar) = dimu_dMass;
	   if ( vName == "dimu_eta" )
	     var_vals.at(iVar) = dimu_eta;
	   if ( vName == "dimu_rapid" )
	     var_vals.at(iVar) = dimu_rapid;

	   if ( vName == "dimu_dR" )
	     var_vals.at(iVar) = dimu_dR;
	   if ( vName == "dimu_dEta" )
	     var_vals.at(iVar) = dimu_dEta;
	   if ( vName == "dimu_dPhi" )
	     var_vals.at(iVar) = dimu_dPhi;
	   if ( vName == "dimu_dPhiS" )
	     var_vals.at(iVar) = dimu_dPhiS;

	   // Jet variables
	   if ( vName == "jet1_pt" )
	     var_vals.at(iVar) = jet1_pt;
	   if ( vName == "jet2_pt" )
	     var_vals.at(iVar) = jet2_pt;
	   if ( vName == "jet1_eta" )
	     var_vals.at(iVar) = jet1_eta;
	   if ( vName == "jet2_eta" )
	     var_vals.at(iVar) = jet2_eta;

	   if ( vName == "dijet1_mass" )
	     var_vals.at(iVar) = dijet1_mass;
	   if ( vName == "dijet2_mass" )
	     var_vals.at(iVar) = dijet2_mass;
	   if ( vName == "dijet3_mass" )
	     var_vals.at(iVar) = dijet3_mass;
	   if ( vName == "dijet4_mass" )
	     var_vals.at(iVar) = dijet4_mass;

	   if ( vName == "dijet1_dEta" )
	     var_vals.at(iVar) = dijet1_dEta;
	   if ( vName == "dijet2_dEta" )
	     var_vals.at(iVar) = dijet2_dEta;
	   if ( vName == "dijet3_dEta" )
	     var_vals.at(iVar) = dijet3_dEta;
	   if ( vName == "dijet4_dEta" )
	     var_vals.at(iVar) = dijet4_dEta;

	   // Global event variables
	   if ( vName == "nJets" )
	     var_vals.at(iVar) = nJets;
	   if ( vName == "nJetsCent" )
	     var_vals.at(iVar) = nJetsCent;
	   if ( vName == "nJetsFwd" )
	     var_vals.at(iVar) = nJetsFwd;
	   if ( vName == "nBMed" )
	     var_vals.at(iVar) = nBMed;

	   if ( vName == "MET" )
	     var_vals.at(iVar) = MET;
	   if ( vName == "MHT" )
	     var_vals.at(iVar) = MHT;
	   if ( vName == "MT_had" )
	     var_vals.at(iVar) = MT_had;
	   if ( vName == "mass_had" )
	     var_vals.at(iVar) = mass_had;


	   /////////////////////////////
	   ///  Spectator variables  ///
	   /////////////////////////////
	   
	   if ( vName == "samp_ID" )
	     var_vals.at(iVar) = samp_ID;
	   if ( vName == "samp_wgt" )
	     var_vals.at(iVar) = samp_wgt;
	   if ( vName == "LHE_HT" )
	     var_vals.at(iVar) = LHE_HT;
	   if ( vName == "dimu_mass" )
	     var_vals.at(iVar) = dimu_mass;
	   if ( vName == "BASE_cat" )
	     var_vals.at(iVar) = BASE_cat;
	   
	 } // End loop: for (UInt_t iVar = 0; iVar < var_names.size(); iVar++)
	     
	 // Unweighted signal and background
	 Double_t sig_evt_weight = 1.0;
	 Double_t bkg_evt_weight = 1.0;
	     
	 // Load values into event
	 if (samp_ID < 0) {
	   if ( (iEvt % 2) == 0 && nTrain_sig < MAX_TR_SIG ) {
	     std::get<1>(factories.at(iFact))->AddSignalTrainingEvent( var_vals, sig_evt_weight );
	     if (iFact == 0) nTrain_sig += 1;
	   } else {
	     std::get<1>(factories.at(iFact))->AddSignalTestEvent( var_vals, sig_evt_weight );
	     if (iFact == 0) nTest_sig += 1;
	   }
	 }
	 if (samp_ID > 0) {
	   if ( (iEvt % 2) == 0 && nTrain_bkg < MAX_TR_BKG ) {
	     std::get<1>(factories.at(iFact))->AddBackgroundTrainingEvent( var_vals, bkg_evt_weight );
	     if (iFact == 0) nTrain_bkg += 1;
	   } else {
	     std::get<1>(factories.at(iFact))->AddBackgroundTestEvent( var_vals, bkg_evt_weight );
	     if (iFact == 0) nTest_bkg += 1;
	   }
	 }


       } // End loop: for (UInt_t iFact = 0; iFact < factories.size(); iFact++) 
       
       nEvt += 1;
       if (samp_ID < 0) nEvt_sig += 1;
       if (samp_ID > 0) nEvt_bkg += 1;
     } // End loop: for (UInt_t iEvt = 0; iEvt < samp->GetEntries(); iEvt++)
   } // End loop: for (int iSamp = 0; iSamp < all_samps.size(); iSamp++)


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

