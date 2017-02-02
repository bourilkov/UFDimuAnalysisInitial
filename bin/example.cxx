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

    TBranch* recoMuonsBranch = s->chain->GetBranch("muons");
    std::cout << "Got muons branch" << std::endl;
    TBranch* recoDimuonsBranch = s->chain->GetBranch("pairs");
    std::cout << "Got Branches." << std::endl;
    recoMuonsBranch->SetAddress(&s->vars.recoMuons);

    recoDimuonsBranch->SetAddress(&s->vars.recoDimuCands);

    Int_t nMuons;
    TBranch* nMuonsBranch = s->chain->GetBranch("nMuons");
    nMuonsBranch->SetAddress(&nMuons);

    for(unsigned int i=0; i<s->N; i++)
    {
       nMuonsBranch->GetEntry(i);
       if(i%1000==0) std::cout << i << " nMuons          : " << nMuons << std::endl;
       recoDimuonsBranch->GetEntry(i);
       recoMuonsBranch->GetEntry(i);
       if(i%1000==0) std::cout << i << " recoMuons->size(): " << s->vars.recoMuons->size() << std::endl;

       if (i > 100000) break;
  
       for(auto& dimu: (*s->vars.recoDimuCands))
       {
           s->vars.dimuCand = &dimu;
           if(i%1000==0) std::cout << i << " dimu.mass       : " << dimu.mass << std::endl;
           if(i%1000==0) std::cout << i << " dimu.mass_      : " << s->vars.dimuCand->mass << std::endl;
           if(i%1000==0) std::cout << i << " dimu.address    : " << &dimu << std::endl;
           if(i%1000==0) std::cout << i << " dimu.address_p  : " << s->vars.dimuCand << std::endl;

           if(i%1000==0) std::cout << i << " pt1: " << s->vars.recoMuons->at(s->vars.dimuCand->iMu1).pt << std::endl;
           if(i%1000==0) std::cout << i << " pt2: " << s->vars.recoMuons->at(s->vars.dimuCand->iMu2).pt << std::endl;
       }
    }

    // test to see if it racks up memory without this
    //delete s->vars.recoDimuCands;
    //delete s->vars.recoMuons;
}
