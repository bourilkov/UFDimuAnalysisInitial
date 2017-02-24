#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "CutSet.h"
#include "Cut.h"
#include "SelectionCuts.h"
#include "CategorySelection.h"
#include "JetSelectionTools.h"
#include "EventTools.cxx"

#include "TLorentzVector.h"

#include <sstream>
#include <map>
#include <vector>
#include <utility>


//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
    }   

    std::map<std::string, Sample*> samples;

    //CutSet cuts();
    //Cut cut();  // Abstract class, can't initialize
  
    ///////////////////////////////////////////////////////////////////
    // SAMPLES---------------------------------------------------------
    ///////////////////////////////////////////////////////////////////
 
    //TString datafilename = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data_from_json/25ns/golden/stage_1_singleMuon_RunCD_GOLDEN_ALL.root");
    //TString datafilename = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data_from_json/25ns/golden/stage_1_singleMuon_RunDBoth_MINIAOD_GOLDEN_ALL.root");
    TString datafilename = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data_from_json/25ns/golden/revised/stage_1_singleMuon_RunDBoth_MINIAOD_GOLDEN_ALL.root");
    Sample* datasample = new Sample(datafilename, "Data", "data");
    datasample->lumi = 2169;
    //datasample->pileupfile = "/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data_from_json/25ns/golden/pileup/PUCalib_Golden_71mb.root";
    datasample->pileupfile = "/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data_from_json/25ns/golden/pileup/old/pileup_data_71mb.root";
    samples["Data"] = datasample;

    ///////////////////////////////////////////////////////////////////
    // PREPROCESSING---------------------------------------------------
    ///////////////////////////////////////////////////////////////////

    // Loop through all of the samples to do some pre-processing
    std::cout << std::endl;
    std::cout << "======== Preprocess the samples... " << std::endl;
    std::cout << std::endl;

    //makePUHistos(samples);

    for(auto const &i : samples)
    {
        // Output some info about the current file
        std::cout << "  /// Looping over " << i.second->name << std::endl;
        std::cout << std::endl;
        std::cout << "    sample name:       " << i.second->name << std::endl;
        std::cout << "    sample file:       " << i.second->filename << std::endl;
        std::cout << "    pileup file:       " << i.second->pileupfile << std::endl;
        std::cout << "    nOriginal:         " << i.second->nOriginal << std::endl;
        std::cout << "    N:                 " << i.second->N << std::endl;
        std::cout << "    nOriginalWeighted: " << i.second->nOriginalWeighted << std::endl;
        std::cout << std::endl;

        if(!i.second->sampleType.Contains("data"))
        {
            // Pileup reweighting
            std::cout << "    +++ PU Reweighting " << i.second->name << "..."  << std::endl;
            std::cout << std::endl;

            i.second->lumiWeights = new reweight::LumiReWeighting(i.second->pileupfile.Data(), samples["Data"]->pileupfile.Data(), "pileup", "pileup");
            std::cout << "        " << i.first << "->lumiWeights: " << i.second->lumiWeights << std::endl;
            std::cout << std::endl;
        }
    }

    ///////////////////////////////////////////////////////////////////
    // Cut and Categorize ---------------------------------------------
    ///////////////////////////////////////////////////////////////////
    
    // Objects to help with the cuts and selections
    TightMuonIdCuts tightMuonId;
    JetSelectionTools jetSelectionTools;
    CategorySelection categorySelection;
    SynchEventSelectionCuts synchEventSelection;
    SynchMuonSelectionCuts synchMuonSelection;

    //synchEventSelection.cutset.cuts[2].on = false;

    // Keep track of the events are in each category
    std::vector<std::pair<int,long long int>> VBFTightVec;
    std::vector<std::pair<int,long long int>> GGFTightVec;
    std::vector<std::pair<int,long long int>> VBFLooseVec;
    std::vector<std::pair<int,long long int>> loose01Vec;
    std::vector<std::pair<int,long long int>> tight01Vec;

    // Keep track of the events from the cross-check list
    std::vector<std::pair<int,long long int>> xmissingVec;
    std::vector<std::pair<int,long long int>> xVBFTightVec;
    std::vector<std::pair<int,long long int>> xGGFTightVec;
    std::vector<std::pair<int,long long int>> xVBFLooseVec;
    std::vector<std::pair<int,long long int>> xloose01Vec;
    std::vector<std::pair<int,long long int>> xtight01Vec;


    // Output information of the first numToPrint events
    int numToPrint = 0;
    int numCounted = 0;

    // std::pair(run, event)
    // The list of events Adrian gave us to look at
    std::vector<std::pair<int,long long int>> eventsToCheck;
    loadEventsFromFile("txt/andrew/xMissing.txt", eventsToCheck);

    // Loop through the events in the data sample
    Sample* s = samples["Data"];
    for(unsigned int i=0; i<s->nOriginal; i++)
    //for(unsigned int i=0; i<10; i++)
    {
       s->getEntry(i); 
       s->vars.validJets = std::vector<TLorentzVector>();
       jetSelectionTools.getValidJetsdR(s->vars, s->vars.validJets);
       std::pair<int,long long int> e(s->vars.eventInfo.run, s->vars.eventInfo.event); // create a pair that identifies the event uniquely

       if(eventInVector(e, eventsToCheck) || i < numToPrint) // Adrian gave a list of events to look at for synch purposes
         outputEvent(s->vars);
       
       // The event is cut
       //if(!tightMuonId.evaluate(s->vars))
       //{ 
       //    if(eventInVector(e, eventsToCheck)) xmissingVec.push_back(e); 
       //    if(eventInVector(e, eventsToCheck))
       //    { 
       //        std::cout << "   !!! Fail Tight Muon ID !!!" << std::endl; 
       //        tightMuonId.cutset.outputFailedCuts();
       //    }
       //    continue; 
       //}
       if(!s->vars.reco1.isTightMuon || !s->vars.reco2.isTightMuon)
       {
           if(eventInVector(e, eventsToCheck))
           { 
               std::cout << "   !!! Fail Tight Muon ID !!!" << std::endl; 
           }
           continue;
       }
       if(!synchEventSelection.evaluate(s->vars))
       { 
           if(eventInVector(e, eventsToCheck)) xmissingVec.push_back(e);
           if(eventInVector(e, eventsToCheck)) 
           { 
               std::cout << "   !!! Fail Event Selection !!!" << std::endl; 
               synchEventSelection.cutset.outputFailedCuts();
           }
           continue; 
       }
       if(!synchMuonSelection.evaluate(s->vars)) 
       {
           if(eventInVector(e, eventsToCheck)) xmissingVec.push_back(e);  
           if(eventInVector(e, eventsToCheck)) 
           { 
               std::cout << "   !!! Fail Muon Selection !!!" << std::endl; 
               synchMuonSelection.cutset.outputFailedCuts();
           }
           continue; 
       }
   
       // The event passes the cuts, see which category it falls into
       categorySelection.evaluate(s->vars);
       if(categorySelection.isVBFTight) 
       {
           VBFTightVec.push_back(e);
           if(eventInVector(e, eventsToCheck)) xVBFTightVec.push_back(e);
       }
       else if(categorySelection.isGGFTight)
       {
           GGFTightVec.push_back(e);
           if(eventInVector(e, eventsToCheck)) xGGFTightVec.push_back(e);
       }
       else if(categorySelection.isVBFLoose)
       {
           VBFLooseVec.push_back(e);
           if(eventInVector(e, eventsToCheck)) xVBFLooseVec.push_back(e);
       }
       else if (categorySelection.isTight01)
       {
           tight01Vec.push_back(e);
           if(eventInVector(e, eventsToCheck)) xtight01Vec.push_back(e);
       }
       else if (categorySelection.isLoose01)
       {
           loose01Vec.push_back(e);
           if(eventInVector(e, eventsToCheck)) xloose01Vec.push_back(e);
       }

       if(eventInVector(e, eventsToCheck) || i < numToPrint) // Adrian gave a list of events to look at for synch purposes
         outputCategory(s->vars, categorySelection);

       // Reset the flags in preparation for the next event
       categorySelection.reset();
    }

    // ////////////////////////////////////////////////////////////////////////////
    // ========= Total Counts =====================================================
    // ////////////////////////////////////////////////////////////////////////////
    
    std::cout << std::endl;
    std::cout << "=========== Category Counts ============" << std::endl;
    std::cout << "VBFTight: " << VBFTightVec.size() << std::endl;
    std::cout << "GGFTight: " << GGFTightVec.size() << std::endl;
    std::cout << "VBFLoose: " << VBFLooseVec.size() << std::endl;
    std::cout << "10Tight: " << tight01Vec.size() <<  std::endl;
    std::cout << "10Loose: " << loose01Vec.size() <<  std::endl;
    std::cout << std::endl;

    int dimutotal = VBFTightVec.size() + GGFTightVec.size() + VBFLooseVec.size() + tight01Vec.size() + loose01Vec.size();
    std::cout << "total: " << dimutotal << std::endl;
    std::cout << std::endl;

    std::cout << "=========== xCheck Category Counts ============" << std::endl;
    std::cout << "VBFTight: " << xVBFTightVec.size() << std::endl;
    std::cout << "GGFTight: " << xGGFTightVec.size() << std::endl;
    std::cout << "VBFLoose: " << xVBFLooseVec.size() << std::endl;
    std::cout << "10Tight:  " << xtight01Vec.size() <<  std::endl;
    std::cout << "10Loose:  " << xloose01Vec.size() <<  std::endl;
    std::cout << "Cut:      " << xmissingVec.size() <<  std::endl;

    int xdimutotal = xVBFTightVec.size() + xGGFTightVec.size() + xVBFLooseVec.size() + xtight01Vec.size() + xloose01Vec.size() + xmissingVec.size();
    std::cout << std::endl;
    std::cout << "total: " << xdimutotal << std::endl;
    std::cout << "to check total: " << eventsToCheck.size() << std::endl;

    std::cout << std::endl;
    //outputEventsToFile(VBFTightVec, "txt/andrew/VBFTight.txt");
    //outputEventsToFile(GGFTightVec, "txt/andrew/GGFTight.txt");
    //outputEventsToFile(VBFLooseVec, "txt/andrew/VBFLoose.txt");
    //outputEventsToFile(tight01Vec, "txt/andrew/10Tight.txt");
    //outputEventsToFile(loose01Vec, "txt/andrew/10Loose.txt");
    //outputEventsToFile(xmissingVec, "txt/andrew/xMissing.txt");


    // Plot some things if we wish
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

    return 0;
}
