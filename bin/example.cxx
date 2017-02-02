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
    gROOT->SetBatch();
    // save the errors for the histogram correctly so they depend upon 
    // the number used to fill originally rather than the scaling
    TH1::SetDefaultSumw2();
    std::map<TString, Sample*> samples;
    getSamples(36814, samples);
    //Sample* s = samples["GGF_AWB"];
    Sample* s = samples["WminusH"];
    //Sample* s = samples["DYJetsToLL_AWB_M100to200"];

   //////////////////////////////////////////////////////////
   // LOOP TTREE -------------------------------------------
   ///////////////////////////////////////////////////////// 

    TBranch* jetsBranch = s->tree->GetBranch("jets");
    TBranch* mhtBranch = s->tree->GetBranch("mht");
    TBranch* recoMuonsBranch = s->tree->GetBranch("muons");
    TBranch* recoDimuonsBranch = s->tree->GetBranch("pairs");
    TBranch* nMuonsBranch = s->tree->GetBranch("nMuons");
    std::cout << "Got Branches." << std::endl;

    Int_t nMuons;
    jetsBranch->SetAddress(&s->vars.jets);
    mhtBranch->SetAddress(&s->vars.mht);
    recoMuonsBranch->SetAddress(&s->vars.recoMuons);
    recoDimuonsBranch->SetAddress(&s->vars.recoDimuCands);
    nMuonsBranch->SetAddress(&nMuons);
    std::cout << "Set Branches." << std::endl;


    for(unsigned int i=0; i<100; i++)
    {
       std::cout << i << " get nMuons" << std::endl;
       nMuonsBranch->GetEntry(i);
       std::cout << i << " get recoDimuons" << std::endl;
       recoDimuonsBranch->GetEntry(i);
       std::cout << i << " get recoMuons" << std::endl;
       recoMuonsBranch->GetEntry(i);
       std::cout << i << " get recoJets" << std::endl;
       jetsBranch->GetEntry(i);
       std::cout << i << " get mht" << std::endl;
       mhtBranch->GetEntry(i);
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

       for(auto& j: (*s->vars.jets))
           std::cout << i << " jet: " << j.outputInfo() << std::endl;

       std::cout << std::endl;
    }

    // test to see if it racks up memory without this
    //delete s->vars.recoDimuCands;
    //delete s->vars.recoMuons;
}
