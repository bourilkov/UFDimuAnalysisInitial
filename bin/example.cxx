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

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{

  std::cout << "Inside example.cxx" << std::endl;
  
    // save the errors for the histogram correctly so they depend upon 
    // the number used to fill originally rather than the scaling
    TH1::SetDefaultSumw2();
    std::map<TString, Sample*> samples;
    std::cout << "\nAbout to get the samples" << std::endl;
    
    GetSamples(samples, "CERN", "RunB");
    GetSamples(samples, "CERN", "tt_ll_MG");
    
    // GetSamples(samples, "CERN", "DATA");
    // GetSamples(samples, "CERN", "SIGNAL");
    // GetSamples(samples, "CERN", "ZJets");
    // GetSamples(samples, "CERN", "ttbar");

    std::cout << "\nGot the samples" << std::endl;
    //Sample* s = samples["GGF_AWB"];
    // Sample* s = samples["DYJetsToLL_AWB"];
    //Sample* s = samples["DYJetsToLL_AWB_M100to200"];
    Sample* s = samples["RunB"];

    std::cout << "Defined the sample" << std::endl;

   //////////////////////////////////////////////////////////
   // LOOP TCHAIN -------------------------------------------
   ///////////////////////////////////////////////////////// 

    // Get branches, set addresses
    Int_t nMuons;
    TBranch* nMuonsBranch = s->chain->GetBranch("nMuons");
    TBranch* recoMuonsBranch = s->chain->GetBranch("muons");
    TBranch* recoDimuonsBranch = s->chain->GetBranch("pairs");
    recoMuonsBranch->SetAddress(&s->vars.recoMuons);
    recoDimuonsBranch->SetAddress(&s->vars.recoDimuCands);
    nMuonsBranch->SetAddress(&nMuons);

    for(unsigned int i=0; i<10; i++)
    {
       nMuonsBranch->GetEntry(i);
       recoDimuonsBranch->GetEntry(i);
       recoMuonsBranch->GetEntry(i);
       jetsBranch->GetEntry(i);
       mhtBranch->GetEntry(i);

       std::vector<TLorentzVector> muons;
       std::vector<TLorentzVector> jets;


       std::cout << "Muons before cleaning... " << std::endl;
       for(auto& m: (*s->vars.recoMuons))
       {
           std::cout << i << " muon: " << m.outputInfo() << std::endl;
           muons.push_back(m.get4vec());
       }
       std::cout << std::endl;

       std::cout << "Jets before cleaning... " << std::endl;
       for(auto& j: (*s->vars.jets))
       {
           std::cout << i << " jet: " << j.outputInfo() << std::endl;
           jets.push_back(j.get4vec());
       }

       std::cout << std::endl;

       c.cleanByDR(jets, muons, 0.3);

       std::cout << "Jets after cleaning... " << std::endl;
       for(auto& j: jets)
       {
           std::cout << i << " jet: " << ParticleTools::output4vecInfo(j) << std::endl;
       }

       std::cout << std::endl;

       std::cout << i << " recoMuons->size(): " << s->vars.recoMuons->size() << std::endl;
       std::cout << i << " jets->size(): " << s->vars.jets->size() << std::endl;
       std::cout << i << " mht.pt: " << s->vars.mht->pt << std::endl;
  
       for(auto& dimu: (*s->vars.recoDimuCands))
       {
           s->vars.dimuCand = &dimu;
           std::cout << i << " dimu.mass_      : " << s->vars.dimuCand->mass << std::endl;
           std::cout << i << " pt1: " << s->vars.recoMuons->at(s->vars.dimuCand->iMu1).pt << std::endl;
           std::cout << i << " pt2: " << s->vars.recoMuons->at(s->vars.dimuCand->iMu2).pt << std::endl;
       }

       std::cout << std::endl;
    }

    // test to see if it racks up memory without this
    //delete s->vars.recoDimuCands;
    //delete s->vars.recoMuons;
}
