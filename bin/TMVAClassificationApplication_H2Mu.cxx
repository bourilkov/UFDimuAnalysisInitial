/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication_H2Mu                                          *
 *                                                                                *
 * This macro provides a simple example on how to use the trained regression MVAs *
 * within an analysis module                                                      *
 **********************************************************************************/

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
#include "TF1.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVAGui.h"
#include "TXMLEngine.h"

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
#include "EventTools.h"
#include "TMVA_helper.h"

//////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////

void loadFromXMLRecursive(TXMLEngine* xml, XMLNodePointer_t node, std::vector<TString>& tvars, std::vector<TString>& svars)
{
    TString node_name = xml->GetNodeName(node);
    //std::cout << "node: " << node_name << std::endl;
   
    // display attributes
    XMLAttrPointer_t attr = xml->GetFirstAttr(node);
    while (attr!=0) 
    {
        TString att_string = xml->GetAttrName(attr);
        TString val_string = xml->GetAttrValue(attr);
        if(node_name == "Variable" && att_string == "Label")
        {
            tvars.push_back(val_string);
            //printf("attr: \"%s\" value: \"%s\"\n", att_string.Data(), val_string.Data());
        }
        if(node_name == "Spectator" && att_string == "Label")
        {
            svars.push_back(val_string);
            //printf("attr: \"%s\" value: \"%s\"\n", att_string.Data(), val_string.Data());
        }
        attr = xml->GetNextAttr(attr);
    }

   // display content (if exists)
   //const char* content = xml->GetNodeContent(node);
   //if (content!=0)
   //    printf("cont: %s\n", content);
   
   // display all child nodes
   XMLNodePointer_t child = xml->GetChild(node);
   while (child!=0) 
   {
       loadFromXMLRecursive(xml, child, tvars, svars);
       child = xml->GetNext(child);
   }
}

//////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////

void getVarNames(TString filename, std::vector<TString>& tvars, std::vector<TString>& svars)
{
    std::cout << Form("/// Reading xmlfile from %s... \n\n", filename.Data());
    // First create the engine.
    TXMLEngine* xml = new TXMLEngine;

    // Now try to parse xml file.
    XMLDocPointer_t xmldoc = xml->ParseFile(filename);
    if (xmldoc==0)
    {    
        delete xml; 
        return;  
    }    

    // Get access to main node of the xml file.
    XMLNodePointer_t mainnode = xml->DocGetRootElement(xmldoc);
   
    // Recursively connect nodes together.
    loadFromXMLRecursive(xml, mainnode, tvars, svars);
   
    // Release memory before exit
    xml->FreeDoc(xmldoc);
    delete xml; 
}

//////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////

void TMVAClassificationApplication_H2Mu( TString myMethodList = "" ) 
{

   /////////////////////////////////////////////////////
   // Get a sample to test the classification on

   std::map<TString, Sample*> sampleMap;
   GetSamples(sampleMap, "UF", "ZJets_MG_incl");
   for(auto& i: sampleMap)
       printf("/// %s in Sample Map\n", i.second->name.Data());

   Sample* s = sampleMap["ZJets_MG"];
   printf("/// Using %s for training\n", s->name.Data());

   TMVA::Tools::Instance();

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication_H2Mu" << std::endl << std::endl;

   TString dir    = "classification/";
   TString weightfile = dir+"f_Opt1_BDTG_default.weights.xml";

   /////////////////////////////////////////////////////
   // Book training and spectator vars into reader

   std::vector<TString> tvars;
   std::map<TString, Float_t> tmap;

   std::vector<TString> svars;
   std::map<TString, Float_t> smap;
   getVarNames(weightfile, tvars, svars);

   // --- Create the Reader object
   // Vars need to have the same names as in the xml
   // and they need to be booked in the same order

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    
   Float_t x = -999;

   std::cout << "/// Booking feature variables ... " << std::endl << std::endl; 
   for(auto& var: tvars)
   {
       std::cout << var << ", ";
       tmap[var] = -999;
       reader->AddVariable(var, &tmap[var]);
   }
   std::cout << std::endl << std::endl;

   std::cout << "/// Listing feature variable map ... " << std::endl << std::endl; 
   for(auto& item: tmap)
   {
       printf("%s: %f\n", item.first.Data(), item.second);
   }

   // dumb... Spectator variables declared in the training have to be added to the reader, too
   std::cout << "/// Booking spectator variables ... " << std::endl << std::endl; 
   for(auto& var: svars)
   {
       std::cout << var << ", ";
       smap[var] = -999;
       reader->AddSpectator(var, &smap[var]);
   }
   std::cout << std::endl << std::endl;

   // Book BDT method into reader, now that the vars are set up
   TString methodName = "BDTG_default";
   reader->BookMVA( methodName, weightfile ); 
   
   /////////////////////////////////////////////////////
   // Loop over the events in the sample

   TStopwatch sw;
   sw.Start();

   // Objects to help with the cuts and selections
   JetCollectionCleaner      jetCollectionCleaner;
   MuonCollectionCleaner     muonCollectionCleaner;
   EleCollectionCleaner      eleCollectionCleaner;

   Run2MuonSelectionCuts  run2MuonSelection;
   Run2EventSelectionCuts run2EventSelection;

   // set some flags
   bool isData = s->sampleType.EqualTo("data");

   std::cout << Form("/// Processing %s \n", s->name.Data());
   s->setBranchAddresses(2);
   for(unsigned int i=0; i<s->N; i++)
   {
      std::cout << Form("  /// GetEntry muon info %s \n", s->name.Data());
      // only load essential information for the first set of cuts 
      s->branches.muPairs->GetEntry(i);
      s->branches.muons->GetEntry(i);

      // loop and find a good dimuon candidate
      if(s->vars.muPairs->size() != 1) continue;
     
      // the dimuon candidate and the muons that make up the pair
      std::cout << Form("  /// Set dimuon info %s \n", s->name.Data());
      MuPairInfo& dimu = s->vars.muPairs->at(0); 
      s->vars.dimuCand = &dimu;
      MuonInfo& mu1 = s->vars.muons->at(s->vars.dimuCand->iMu1);
      MuonInfo& mu2 = s->vars.muons->at(s->vars.dimuCand->iMu2);

      ////////////////////////////////////////////////////////////////////
      // CUTS  ----------------------------------------------------------
      ///////////////////////////////////////////////////////////////////
      
      std::cout << Form("  /// Cuts %s \n", s->name.Data());
      if(dimu.mass < 110 || dimu.mass > 160)
      {
          continue;
      }
      if(!run2EventSelection.evaluate(s->vars))
      {
          continue;
      }
      if(!mu1.isMediumID || !mu2.isMediumID)
      {
          continue;
      }
      if(!run2MuonSelection.evaluate(s->vars))
      {
          continue;
      }

      std::cout << Form("  /// Get remaining branches %s \n", s->name.Data());
      s->branches.getEntry(i);

      // clear vectors for the valid collections
      s->vars.validMuons.clear();
      s->vars.validExtraMuons.clear();
      s->vars.validElectrons.clear();
      s->vars.validJets.clear();
      s->vars.validBJets.clear();

      std::cout << Form("  /// Load valid collections %s \n", s->name.Data());
      // load valid collections from s->vars raw collections
      jetCollectionCleaner.getValidJets(s->vars, s->vars.validJets, s->vars.validBJets);
      muonCollectionCleaner.getValidMuons(s->vars, s->vars.validMuons, s->vars.validExtraMuons);
      eleCollectionCleaner.getValidElectrons(s->vars, s->vars.validElectrons);

      // Clean jets and electrons from muons, then clean remaining jets from remaining electrons
      CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validMuons, 0.4);
      CollectionCleaner::cleanByDR(s->vars.validElectrons, s->vars.validMuons, 0.4);
      CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validElectrons, 0.4);

      std::cout << Form("  /// Set map values %s \n", s->name.Data());
      for(auto& v: tmap)
      {
          tmap[v.first] = s->vars.getValue(v.first.Data());
          printf("  %s: %f\n", v.first.Data(), v.second);
      }

      std::cout << Form("  /// Get BDT Output %s \n", s->name.Data());
      Float_t val = (reader->EvaluateRegression(methodName))[0];
      printf("  !!! %d) BDT_Prediction: %f\n", i, val);
      break;
   }
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   delete reader;
    
   std::cout << "==> TMVAClassificationApplication_H2Mu is done!" << std::endl << std::endl;
}

//////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////

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
   TMVAClassificationApplication_H2Mu(methodList);
   return 0;
}
