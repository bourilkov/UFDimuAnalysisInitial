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
    TFile* f = new TFile("HToMuMu_NTuple_1.root");
    TTree* t = (TTree*) f->Get("dimuons/tree");
    std::cout << t->GetEntries() << std::endl;

   //////////////////////////////////////////////////////////
   // LOOP TTREE -------------------------------------------
   ///////////////////////////////////////////////////////// 

    TBranch* recoMuonsBranch = t->GetBranch("muons");
    std::cout << "Got Branch." << std::endl;
    std::vector<MuonInfo>* recoMuons = 0;
    recoMuonsBranch->SetAddress(&recoMuons);

    for(unsigned int i=0; i<10; i++)
    {
       recoMuonsBranch->GetEntry(i);
       std::cout << i << " pt: " << recoMuons->at(0).pt << std::endl;
    }
}
