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
    
    GetSamples(samples, "UF", "DATA");

    std::cout << "\nGot the samples" << std::endl;
    Sample* s = samples["RunB"];

    std::cout << "Defined the sample" << std::endl;

   //////////////////////////////////////////////////////////
   // LOOP TCHAIN -------------------------------------------
   ///////////////////////////////////////////////////////// 

    // Get branches, set addresses
    s->setBranchAddresses(1);

    for(unsigned int i=0; i<10; i++)
    {
       s->branches.recoDimuCands->GetEntry(i);
       s->branches.recoMuons->GetEntry(i);

       float pto = s->vars.recoMuons->at(0).pt;
       //float ptb = s->branches.recoMuons->GetLeaf("pt")->GetValue(0);
       float ptb = 0;

       std::cout << "CHECK GET BRANCH: " << pto << ", " << ptb << std::endl;

       if(s->vars.recoDimuCands > 0)
       {
           PairInfo& dimu = s->vars.recoDimuCands->at(0); 
           std::cout << "before0 dimu_pt   : " << dimu.pt << std::endl;
           std::cout << "before0 dimu_pt_R : " << dimu.pt_Roch << std::endl;
           std::cout << "before1 dimu_pt   : " << s->vars.recoDimuCands->at(0).pt << std::endl;
           std::cout << "before1 dimu_pt_R : " << s->vars.recoDimuCands->at(0).pt_Roch << std::endl;
           std::cout << std::endl;
           dimu.pt = dimu.pt_Roch;
           std::cout << "after0 dimu_pt    : " << dimu.pt << std::endl;
           std::cout << "after1 dimu_pt    : " << s->vars.recoDimuCands->at(0).pt << std::endl;
           std::cout << std::endl;
       }

       for(auto& dimu: (*s->vars.recoDimuCands))
       {
           std::cout << i << " dimu: " << dimu.outputInfo() << std::endl;
       }
       std::cout << std::endl;

       for(auto& m: (*s->vars.recoMuons))
       {
           std::cout << i << " muon: " << m.outputInfo() << std::endl;
       }
       std::cout << std::endl;
    }

    // test to see if it racks up memory without this
    //delete s->vars.recoDimuCands;
    //delete s->vars.recoMuons;
}
