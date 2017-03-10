///////////////////////////////////////////////////////////////////////////
//                           EventTools.cxx                              //
//=======================================================================//
//                                                                       //
//        Miscellaneous tools: output event info to terminal             //
//        output run,event from vector to CSV, output vars and values    //
//        from map to CSV, etc                                           //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <vector>
#include <utility> // pair
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>

#include "../bin/SignificanceMetrics.hxx"
#include "EventTools.h"
#include "ParticleTools.h"
#include "DiMuPlottingSystem.h"
#include "CategorySelection.h"
#include "JetCollectionCleaner.h"

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

void EventTools::loadEventsFromFile(TString filename, std::vector<std::pair<int, long long int>>& v)
{
// Used for the synchronization exercise. Run, event were saved to different csv files.
// Load the run, event from the csv into the vector.

    std::ifstream file(filename);
    std::string run;
    std::string event;
    int runn; 
    long long int eventn;
    while(std::getline(file, run, ','))
    {
        std::getline(file, event);
        std::stringstream ssrun;
        std::stringstream ssevent;
        ssrun << run;
        ssrun >> runn;
        ssevent << event;
        ssevent >> eventn; 

        v.push_back(std::pair<int, long long int>(runn, eventn)); 
    } 
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

void EventTools::outputEventsToFile(std::vector< std::pair<int,long long int> >& v, TString filename)
{
// used for the synchronization exercise
// output the run, event info from the vector into a csv file

    std::cout << "Exporting events to " << filename << " ..." << std::endl;
    std::ofstream file(filename, std::ofstream::out);
    for(auto const &e : v)
    {   
        file << e.first << ", " << e.second << std::endl;
    }   
    file.close();  
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

TString EventTools::outputMapKeysCSV(std::map<TString,double>& map)
{
// Given a map output the keys to a csv string, used to output the variable names
// then later the variable values in csv to create a dataset useable by other
// machine learning algorithms outside of ROOT
    TString out = "";
    for(auto const &c : map)
    {   
        out+=","+c.first;
    }
    out = out(1,out.Length());
    return out;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

TString EventTools::outputMapValuesCSV(std::map<TString,double>& map)
{
// used to output the values of the map to csv string
// easy to export the var values to use this dataset outside ROOT
    TString out = "";
    for(auto const &c : map)
    {   
        std::stringstream ss;
        ss << c.second;
        out+=","+ss.str();
    }
    out = out(1,out.Length());
    return out;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

bool EventTools::sameRunAndEvent(std::pair<int,long long int> a, std::pair<int,long long int> b)
{
// are the run and event number the same for both?
   if(a.first == b.first && a.second == b.second) return true; 
   else return false;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

bool EventTools::eventInVector(std::pair<int,long long int> e, std::vector<std::pair<int,long long int>> events)
{
// see if the same run, event number is in the vector 
    for(unsigned int i=0; i<events.size(); i++)
    {
        if(sameRunAndEvent(e, events[i])) return true;
    }
    return false;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

void EventTools::outputEvent(VarSet& vars)
{
// output all of the information for an event, should make this more modular
         JetCollectionCleaner js;
         if(vars.eventInfo !=0) 
             std::cout << "========= RUN: " << vars.eventInfo->run << ", EVENT: " << vars.eventInfo->event << " =============" << std::endl;
         std::cout << std::endl;
         std::cout << "  dimu_mass: " << vars.dimuCand->mass << std::endl;
         std::cout << "  dimu_mass_PF: " << vars.dimuCand->mass_PF << std::endl;
         std::cout << "  dimu_pt: " << vars.dimuCand->pt << std::endl;
         std::cout << std::endl;
         std::cout << "  mu1: " << vars.recoMuons->at(vars.dimuCand->iMu1).outputInfo() << std::endl;
         std::cout << "  mu2: " << vars.recoMuons->at(vars.dimuCand->iMu2).outputInfo() << std::endl;
         std::cout << std::endl;
         std::cout << "  nJets: " << vars.jets->size() << std::endl;
         std::cout << "  nValidJets: " << vars.validJets.size() << std::endl;
         std::cout << "  nValidBJets: " << vars.validBJets.size() << std::endl;
         std::cout << std::endl;
         std::cout << "  nValidMuons: " << vars.validMuons.size() << std::endl;
         std::cout << "  nValidExtraMuons: " << vars.validExtraMuons.size() << std::endl;
         std::cout << "  nValidElectrons: " << vars.validElectrons.size() << std::endl;
         std::cout << "  nValidExtraLeptons: " << vars.validElectrons.size() + vars.validExtraMuons.size() << std::endl;
         std::cout << std::endl;
         std::cout << "  MET: " << vars.mht->pt << std::endl;
         std::cout << std::endl;
         //std::cout << "  nGenJets: " << vars.genJets.nJets << std::endl;
         //std::cout << "  nValidGenJets: " << vars.validGenJets.size() << std::endl;
         //std::cout << std::endl;

         for(unsigned int j=0; j<vars.jets->size(); j++)
         {
             std::cout << "  jet" << j << ": " <<  vars.jets->at(j).outputInfo() << std::endl;
             if(j==vars.jets->size()-1) std::cout << std::endl;
         }

         for(unsigned int j=0; j<vars.validJets.size(); j++)
         {
             std::cout << "  validjet" << j << ": " << ParticleTools::output4vecInfo(vars.validJets[j]) << std::endl;
             if(j==vars.validJets.size()-1) std::cout << std::endl;
         }

         for(unsigned int j=0; j<vars.validBJets.size(); j++)
         {
             std::cout << "  validBjet" << j <<  ": " << ParticleTools::output4vecInfo(vars.validBJets[j]) << std::endl;
             if(j==vars.validBJets.size()-1) std::cout << std::endl;
         }

         for(unsigned int j=0; j<vars.validMuons.size(); j++)
         {
             std::cout << "  validMuon" << j <<  ": " << ParticleTools::output4vecInfo(vars.validMuons[j]) << std::endl;
             if(j==vars.validMuons.size()-1) std::cout << std::endl;
         }

         for(unsigned int j=0; j<vars.validExtraMuons.size(); j++)
         {
             std::cout << "  validExtraMuon" << j <<  ": " << ParticleTools::output4vecInfo(vars.validExtraMuons[j]) << std::endl;
             if(j==vars.validExtraMuons.size()-1) std::cout << std::endl;
         }

         for(unsigned int j=0; j<vars.validElectrons.size(); j++)
         {
             std::cout << "  validElectron" << j <<  ": " << ParticleTools::output4vecInfo(vars.validElectrons[j]) << std::endl;
             if(j==vars.validElectrons.size()-1) std::cout << std::endl;
         }
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

void EventTools::outputEvent(VarSet& vars, Categorizer& categorizer)
{
// output the event and its categorization
         // output standard information about the event
         outputEvent(vars);
         categorizer.outputResults();
}
