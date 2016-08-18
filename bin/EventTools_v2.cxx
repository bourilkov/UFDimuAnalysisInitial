// EventTools.cxx

#ifndef EVENTTOOLS
#define EVENTTOOLS

#include <vector>
#include <utility> // pair
#include <map>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>

#include "TString.h"
#include "VarSet.h"
#include "CategorySelection_v2.h"
#include "JetSelectionTools.h"

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

void loadEventsFromFile(TString filename, std::vector<std::pair<int, long long int>>& v)
{
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

void outputEventsToFile(std::vector<std::pair<int,long long int>>& v, TString filename)
{
    std::cout << "  /// Exporting events to " << filename << " ..." << std::endl;
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

bool sameRunAndEvent(std::pair<int,long long int> a, std::pair<int,long long int> b)
{
   if(a.first == b.first && a.second == b.second) return true; 
   else return false;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

bool eventInVector(std::pair<int,long long int> e, std::vector<std::pair<int,long long int>> events)
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

void outputEvent(VarSet& vars)
{
         JetSelectionTools js;
         std::cout << "========= RUN: " << vars.eventInfo.run << ", EVENT: " << vars.eventInfo.event << " =============" << std::endl;
         std::cout << std::endl;
         std::cout << "  recoCandMass: " << vars.recoCandMass << std::endl;
         std::cout << "  recoCandMassPF: " << vars.recoCandMassPF << std::endl;
         std::cout << "  recoCandPt: " << vars.recoCandPt << std::endl;
         std::cout << std::endl;
         std::cout << "  reco1.pt: " << vars.reco1.pt << std::endl;
         std::cout << "  reco1.phi: " << vars.reco1.phi << std::endl;
         std::cout << "  reco1.eta: " << vars.reco1.eta << std::endl;
         std::cout << "  reco1.isHltMatched: " << (vars.reco1.isHltMatched[0] || vars.reco2.isHltMatched[1]) << std::endl;
         std::cout << std::endl;
         std::cout << "  reco2.pt: " << vars.reco2.pt << std::endl;
         std::cout << "  reco2.phi: " << vars.reco2.phi << std::endl;
         std::cout << "  reco2.eta: " << vars.reco2.eta << std::endl;
         std::cout << "  reco2.isHltMatched: " << (vars.reco2.isHltMatched[0] || vars.reco2.isHltMatched[1]) << std::endl;
         std::cout << std::endl;
         std::cout << "  nJets: " << vars.jets.nJets << std::endl;
         std::cout << "  nValidJets: " << vars.validJets.size() << std::endl;
         std::cout << "  MET: " << vars.met.pt << std::endl;
         std::cout << std::endl;

         for(unsigned int j=0; j<vars.jets.nJets && j<10; j++)
         {
             std::cout << "  jet" << j << " pt: " <<  vars.jets.pt[j] << std::endl;
             std::cout << "  jet" << j << " phi: " << vars.jets.phi[j] << std::endl;
             std::cout << "  jet" << j << " eta: " << vars.jets.eta[j] << std::endl;
             std::cout << "  jet" << j << " dR1: " << js.dR(vars.jets.eta[j], vars.jets.phi[j], vars.reco1.eta, vars.reco1.phi) << std::endl;
             std::cout << "  jet" << j << " dR2: " << js.dR(vars.jets.eta[j], vars.jets.phi[j], vars.reco2.eta, vars.reco2.phi) << std::endl;
             std::cout << std::endl;

         }

         std::cout << std::endl;

         for(unsigned int j=0; j<vars.validJets.size(); j++)
         {
             std::cout << "  validjet" << j << " pt: " <<  vars.validJets[j].Pt() << std::endl;
             std::cout << "  validjet" << j << " phi: " << vars.validJets[j].Phi() << std::endl;
             std::cout << "  validjet" << j << " eta: " << vars.validJets[j].Eta() << std::endl;
             std::cout << "  validjet" << j << " dR1: " << js.dR(vars.validJets[j].Eta(), vars.validJets[j].Phi(), vars.reco1.eta, vars.reco1.phi) << std::endl;
             std::cout << "  validjet" << j << " dR2: " << js.dR(vars.validJets[j].Eta(), vars.validJets[j].Phi(), vars.reco2.eta, vars.reco2.phi) << std::endl;
             std::cout << std::endl;

         }
         std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

void outputEvent(VarSet& vars, Categorizer& categorizer)
{
         // output standard information about the event
         outputEvent(vars);
         categorizer.outputResults();
}

#endif
