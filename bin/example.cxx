#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "EventSelection.h"
#include "MuonSelection.h"
#include "CategorySelection.h"
#include "JetCollectionCleaner.h"
#include "MuonCollectionCleaner.h"
#include "EleCollectionCleaner.h"

#include "EventTools.h"
#include "PUTools.h"
#include "SignificanceMetrics.hxx"

#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TROOT.h"

#include "SampleDatabase.cxx"
#include "MuonInfo.h"
#include "EleInfo.h"

#include <sstream>
#include <map>
#include <vector>
#include <utility>
#include <unordered_map>

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

double function()
{
    return 1.0;
}

int main(int argc, char* argv[])
{

  std::cout << "Inside example.cxx" << std::endl;
  
    // save the errors for the histogram correctly so they depend upon 
    // the number used to fill originally rather than the scaling
    TH1::SetDefaultSumw2();
    std::map<TString, Sample*> samples;
    std::cout << "\nAbout to get the samples" << std::endl;
    
    GetSamples(samples, "UF", "DATA");

    std::cout << "\nGot the samples" << std::endl;
    Sample* s = samples["RunB"];

    std::cout << "Defined the sample" << std::endl;

   //////////////////////////////////////////////////////////
   // LOOP TCHAIN -------------------------------------------
   ///////////////////////////////////////////////////////// 

    // Get branches, set addresses
    s->setBranchAddresses(2);

    XMLCategorizer xmlc;
    xmlc.loadFromXML("xml/tree.xml");
    //std::cout << std::endl;
    //xmlc.outputResults();

    // Objects to help with the cuts and selections
    JetCollectionCleaner      jetCollectionCleaner;
    MuonCollectionCleaner     muonCollectionCleaner;
    EleCollectionCleaner      eleCollectionCleaner;

    Run2MuonSelectionCuts     run2MuonSelection;
    Run2EventSelectionCuts80X run2EventSelectionMC;

    //std::map<TString, double(*)()> m;
    //std::unordered_map<std::string, double(*)()> um;
    //um["f"] = &function;

    //std::cout << um["f"] << std::endl;

    for(unsigned int i=0; i<10; i++)
    {
       xmlc.reset();

       s->branches.recoDimuCands->GetEntry(i);
       s->branches.recoMuons->GetEntry(i);

       if(s->vars.recoDimuCands > 0)
       {
           PairInfo& dimu = s->vars.recoDimuCands->at(0); 
           s->vars.dimuCand = &dimu;
           MuonInfo& mu1 = s->vars.recoMuons->at(dimu.iMu1);
           MuonInfo& mu2 = s->vars.recoMuons->at(dimu.iMu2);

           // Load the rest of the information needed for run2 categories
           s->branches.jets->GetEntry(i);
           s->branches.mht->GetEntry(i);
           s->branches.nVertices->GetEntry(i);
           s->branches.recoElectrons->GetEntry(i);

           if(s->sampleType != "data")
           {
               //s->branches.gen_wgt->GetEntry(i);
               //s->branches.nPU->GetEntry(i);
               //s->branches.pu_wgt->GetEntry(i);
               //s->branches.eff_wgt->GetEntry(i);
           }

           // clear vectors for the valid collections
           s->vars.validMuons.clear();
           s->vars.validExtraMuons.clear();
           s->vars.validElectrons.clear();
           s->vars.validJets.clear();
           s->vars.validBJets.clear();

           // load valid collections from s->vars raw collections
           jetCollectionCleaner.getValidJets(s->vars, s->vars.validJets, s->vars.validBJets);
           muonCollectionCleaner.getValidMuons(s->vars, s->vars.validMuons, s->vars.validExtraMuons);
           eleCollectionCleaner.getValidElectrons(s->vars, s->vars.validElectrons);

           // Clean jets and electrons from muons, then clean remaining jets from remaining electrons
           CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validMuons, 0.3);
           CollectionCleaner::cleanByDR(s->vars.validElectrons, s->vars.validMuons, 0.3);
           CollectionCleaner::cleanByDR(s->vars.validJets, s->vars.validElectrons, 0.3);

           //std::cout << "dimu_pt   : " << dimu.pt_PF << std::endl;
           //std::cout << "dimu_pt_f : " << s->vars.dimu_pt() << std::endl;
           //std::cout << "dimu_pt_m : " << s->vars.getValue("dimu_pt") << std::endl;
           //std::cout << "dimu_pt_mt: " << s->vars.getValue(TString("dimu_pt").Data()) << std::endl;
                      
           std::cout << std::endl;
           xmlc.evaluate(s->vars);
           EventTools::outputEvent(s->vars);
           xmlc.outputResults();
           std::cout << std::endl;
       }

       //for(auto& dimu: (*s->vars.recoDimuCands))
       //{
       //    std::cout << i << " dimu: " << dimu.outputInfo() << std::endl;
       //}
       //std::cout << std::endl;

       //for(auto& m: (*s->vars.recoMuons))
       //{
       //    std::cout << i << " muon: " << m.outputInfo() << std::endl;
       //}
       //std::cout << std::endl;
    }

    // test to see if it racks up memory without this
    //delete s->vars.recoDimuCands;
    //delete s->vars.recoMuons;
}
