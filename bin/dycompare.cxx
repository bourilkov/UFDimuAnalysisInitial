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
  std::map<TString, Sample*> sampleMap;
  std::vector<Sample*> sampleVec;
  getSamples(36814, sampleMap);
  //Sample* s = samples["GGF_AWB"];
  sampleVec.push_back(sampleMap["DYJetsToLL_AWB"]);
  //sampleVec.push_back(sampleMap["DYJetsToLL_AWBi_madgraph"]);
  sampleVec.push_back(sampleMap["DYJetsToLL_AWB_M100to200"]);
  sampleVec.push_back(sampleMap["DYJetsToLL_AWB_M100to200_SpringPU"]);

  // store histograms
  std::map<TString, TH1F*> histMap;

  for(auto& s: sampleVec)
  {
    std::cout << "/// Looping over " << s->name << std::endl;

    TBranch* recoMuonsBranch = s->tree->GetBranch("muons");
    TBranch* recoDimuonsBranch = s->tree->GetBranch("pairs");
    TBranch* genWeightBranch = s->tree->GetBranch("GEN_wgt");

    recoMuonsBranch->SetAddress(&s->vars.recoMuons);
    recoDimuonsBranch->SetAddress(&s->vars.recoDimuCands);
    genWeightBranch->SetAddress(&s->vars.gen_wgt);

    std::cout << "Got Branches." << std::endl;

    // maps sample to histogram
    histMap[s->name] = new TH1F(s->name, s->name, 110, 100, 200);

    for(unsigned int i=0; i<s->N; i++)
    {
       recoDimuonsBranch->GetEntry(i);

       if(s->vars.recoDimuCands->size() == 0) continue;

       recoMuonsBranch->GetEntry(i);
       genWeightBranch->GetEntry(i);
  
       for(auto& dimu: (*s->vars.recoDimuCands))
       {
           s->vars.dimuCand = &dimu;
           if(s->vars.dimuCand->mass_PF > 200 || s->vars.dimuCand->mass_PF < 110) continue;
           histMap[s->name]->Fill(s->vars.dimuCand->mass_PF, s->vars.gen_wgt);
       }
    }

    // test to see if it racks up memory without this
    //delete s->vars.recoDimuCands;
    //delete s->vars.recoMuons;
  }
  TFile* f = new TFile("compare_drell_yan.root", "RECREATE");
  f->cd();

  for(auto& h: histMap)
      h.second->Write();

  f->Write();
  f->Close();
}
