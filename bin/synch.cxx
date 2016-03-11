#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "CutSet.h"
#include "Cut.h"
#include "SelectionCuts.h"
#include "CategorySelection.h"
#include "JetSelectionTools.h"

#include "TLorentzVector.h"

#include <sstream>
#include <map>
#include <vector>
#include <utility>

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

bool sameRunAndEvent(std::pair<int,int> a, std::pair<int,int> b)
{
   if(a.first == b.first && a.second == b.second) return true;
   else return false;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

bool eventInVector(std::pair<int,int> e, std::vector<std::pair<int,int>> events)
{
    for(unsigned int i=0; i<events.size(); i++)
    {
        if(sameRunAndEvent(e, events[i])) return true;
    }
    return false;
}

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
 
    TString datafilename = TString("/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data_from_json/25ns/golden/stage_1_singleMuon_RunCD_GOLDEN_ALL.root");
    Sample* datasample = new Sample(datafilename, "Data", "data");
    datasample->lumi = 2260;
    datasample->pileupfile = "/cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1/data_from_json/25ns/golden/pileup/PUCalib_Golden_71mb.root";
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
    // Do Stuff -------------------------------------------------------
    ///////////////////////////////////////////////////////////////////
    
    // Objects to help with the cuts and selections
    TightMuonIdCuts tightMuonId;
    JetSelectionTools jetSelectionTools;
    CategorySelection categorySelection;
    SynchEventSelectionCuts synchEventSelection;
    SynchMuonSelectionCuts synchMuonSelection;

    // Count the number of events in each category
    int loose01 = 0;
    int tight01 = 0;
    int VBFTight = 0;
    int GGFTight = 0;
    int VBFLoose = 0;
    int n2jets = 0;

    // Output information of the first numToPrint events
    int numToPrint = 0;

    // std::pair(run, event)
    // The list of events Adrian gave us to look at
    std::vector<std::pair<int,int>> eventsToCheck;
    eventsToCheck.push_back( std::pair<int,int>(257751, 618680179) );
    //eventsToCheck.push_back( std::pair<int,int>(257645, 309094265) );
    //eventsToCheck.push_back( std::pair<int,int>(256843, 570071294) );
    //eventsToCheck.push_back( std::pair<int,int>(258158, 979723330) );
    //eventsToCheck.push_back( std::pair<int,int>(257751, 618680179) );
    //eventsToCheck.push_back( std::pair<int,int>(256843, 322474699) );

    // Loop through the events in the data sample
    Sample* s = samples["Data"];
    for(unsigned int i=0; i<s->nOriginal; i++)
    //for(unsigned int i=0; i<10; i++)
    {
       s->getEntry(i); 
       s->vars.validJets = std::vector<TLorentzVector>();
       jetSelectionTools.getValidJets(s->vars.jets, s->vars.validJets);
       std::pair<int,int> e(s->vars.eventInfo.run, s->vars.eventInfo.event); // create a pair that identifies the event uniquely

       if(eventInVector(e, eventsToCheck) || i < numToPrint) // Adrian gave a list of events to look at for synch purposes
       //if(i < numToPrint) 
       {
           std::cout << std::endl;
           std::cout << "========= RUN: " << e.first << ", EVENT: " << e.second << " =============" << std::endl;
           std::cout << "  Count #: " << i << std::endl;
           std::cout << std::endl;
           std::cout << "  recoCandMass: " << s->vars.recoCandMass << std::endl;
           std::cout << "  recoCandPt: " << s->vars.recoCandPt << std::endl;
           std::cout << std::endl;
           std::cout << "  reco1.pt: " << s->vars.reco1.pt << std::endl;
           std::cout << "  reco1.phi: " << s->vars.reco1.phi << std::endl;
           std::cout << "  reco1.eta: " << s->vars.reco1.eta << std::endl;
           std::cout << std::endl;
           std::cout << "  reco2.pt: " << s->vars.reco2.pt << std::endl;
           std::cout << "  reco2.phi: " << s->vars.reco2.phi << std::endl;
           std::cout << "  reco2.eta: " << s->vars.reco2.eta << std::endl;
           std::cout << std::endl;
           std::cout << "  nJets: " << s->vars.jets.nJets << std::endl;
           std::cout << "  nValidJets: " << s->vars.validJets.size() << std::endl;
           std::cout << "  MET: " << s->vars.met.pt << std::endl;
           std::cout << std::endl;

           for(unsigned int j=0; j<s->vars.jets.nJets && j<10; j++)
           {
               std::cout << "  jet" << j << " pt: " <<  s->vars.jets.pt[j] << std::endl;
               std::cout << "  jet" << j << " phi: " << s->vars.jets.phi[j] << std::endl;
               std::cout << "  jet" << j << " eta: " << s->vars.jets.eta[j] << std::endl;
               std::cout << std::endl;

           }

           std::cout << std::endl;

           for(unsigned int j=0; j<s->vars.validJets.size(); j++)
           {
               std::cout << "  validjet" << j << " pt: " <<  s->vars.validJets[j].Pt() << std::endl;
               std::cout << "  validjet" << j << " phi: " << s->vars.validJets[j].Phi() << std::endl;
               std::cout << "  validjet" << j << " eta: " << s->vars.validJets[j].Eta() << std::endl;
               std::cout << std::endl;

           }
           std::cout << std::endl;
       }

       if(!tightMuonId.evaluate(s->vars)) continue;
       if(!synchEventSelection.evaluate(s->vars)) continue;
       if(!synchMuonSelection.evaluate(s->vars)) continue;

       categorySelection.evaluate(s->vars);
       if(categorySelection.isVBFTight) VBFTight++;
       else if(categorySelection.isGGFTight) GGFTight++;
       else if(categorySelection.isVBFLoose) VBFLoose++;
       else if (categorySelection.isTight01) tight01++;
       else if (categorySelection.isLoose01) loose01++;

       categorySelection.reset();
    }

    // ////////////////////////////////////////////////////////////////////////////
    // ========= Total Counts =====================================================
    // ////////////////////////////////////////////////////////////////////////////
    
    std::cout << std::endl;
    std::cout << "=========== Category Counts ============" << std::endl;
    std::cout << "VBFTight: " << VBFTight << std::endl;
    std::cout << "GGFTight: " << GGFTight << std::endl;
    std::cout << "VBFLoose: " << VBFLoose << std::endl;
    std::cout << "10Tight: " << tight01 << std::endl;
    std::cout << "10Loose: " << loose01 << std::endl;

    int dimutotal = VBFTight + GGFTight + VBFLoose + tight01 + loose01;
    std::cout << "total: " << dimutotal << std::endl;

    DiMuPlottingSystem* dps = new DiMuPlottingSystem();

    return 0;
}
