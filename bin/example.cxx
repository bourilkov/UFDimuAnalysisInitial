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
    // save the errors for the histogram correctly so they depend upon 
    // the number used to fill originally rather than the scaling
    TH1::SetDefaultSumw2();
    std::map<TString, Sample*> samples;
    getSamples(33598, samples);
    Sample* s = samples["GGF_AWB"];

   //////////////////////////////////////////////////////////
   // LOOP TTREE -------------------------------------------
   ///////////////////////////////////////////////////////// 

    TBranch* recoMuonsBranch = s->tree->GetBranch("muons");
    TBranch* recoDimuonsBranch = s->tree->GetBranch("pairs");
    std::cout << "Got Branches." << std::endl;
    recoMuonsBranch->SetAddress(&s->vars.recoMuons);

    recoDimuonsBranch->SetAddress(&s->vars.recoDimuCands);

    Int_t nMuons;
    TBranch* nMuonsBranch = s->tree->GetBranch("nMuons");
    nMuonsBranch->SetAddress(&nMuons);

    for(unsigned int i=0; i<10; i++)
    {
       nMuonsBranch->GetEntry(i);
       std::cout << i << " nMuons          : " << nMuons << std::endl;
       recoDimuonsBranch->GetEntry(i);
       recoMuonsBranch->GetEntry(i);
       std::cout << i << " recoMuons->size(): " << s->vars.recoMuons->size() << std::endl;
  
       for(auto& dimu: (*s->vars.recoDimuCands))
       {
           s->vars.dimuCand = &dimu;
           std::cout << i << " dimu.mass       : " << dimu.mass << std::endl;
           std::cout << i << " dimu.mass_      : " << s->vars.dimuCand->mass << std::endl;
           std::cout << i << " dimu.address    : " << &dimu << std::endl;
           std::cout << i << " dimu.address_p  : " << s->vars.dimuCand << std::endl;

           std::cout << i << " pt1: " << s->vars.recoMuons->at(s->vars.dimuCand->iMu1).pt << std::endl;
           std::cout << i << " pt2: " << s->vars.recoMuons->at(s->vars.dimuCand->iMu2).pt << std::endl;
       }
    }
}
