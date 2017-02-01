// EventTools.cxx

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

void EventTools::outputEventsToFile(std::vector<std::pair<int,long long int>>& v, TString filename)
{
// used for the synchronization exercise
// output the run, event info from the vector into a csv file

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

void EventTools::outputCategoryCountsToFile(Categorizer& categories, TString filename)
{
// Given a categorizer, output the expected net signal, bkg to a csv

    std::cout << "  /// Exporting events to " << filename << " ..." << std::endl;
    std::ofstream file(filename, std::ofstream::out);

    TString fieldnames = "Category, Significance, Signal_Net, Bkg_Net";

    // Loop through signal samples to get names
    TIter nextSig(categories.categoryMap["ALL"].signalList);
    TObject* object = 0;

    while ((object = nextSig()))
    {
        TH1D* h = (TH1D*) object;
        TString samplename = h->GetName();
        samplename = samplename.ReplaceAll("ALL_", ""); // remove "ALL_" and we are left with the sample name
        fieldnames += ", " + samplename; 
    }
    
    // Loop through bkg samples to get names
    TIter nextBkg(categories.categoryMap["ALL"].bkgList);
    object = 0;

    while ((object = nextBkg()))
    {
        TH1D* h = (TH1D*) object;
        TString samplename = h->GetName();
        samplename = samplename.ReplaceAll("ALL_", ""); // remove "ALL_" and we are left with the sample name
        fieldnames += ", " + samplename; 
    }

    // put the titles of the csv fields on the first line
    file << fieldnames.Data() << std::endl;

    // now get the info for each category    
    for(auto const &c : categories.categoryMap)
    {   
        TH1D* hsig = DiMuPlottingSystem::addHists(c.second.signalList, c.first+"_Net_Signal", c.first+"_Net_Signal");
        TH1D* hbkg = DiMuPlottingSystem::addHists(c.second.bkgList, c.first+"_Net_Bkg", c.first+"_Net_Bkg");

        double sintegral = hsig->Integral(0, hsig->GetSize());
        double bintegral = hbkg->Integral(0, hbkg->GetSize());

        PoissonSignificance poisson0(0);
        double significance = poisson0.significance(sintegral, bintegral);

        file << c.first << ", " << significance << ", " << sintegral << ", " << bintegral;

        // Loop through signal samples to get names
        TIter nextSig2(c.second.signalList);
        object = 0;

        while ((object = nextSig2()))
        {
            TH1D* h = (TH1D*) object;
            file << ", " << h->Integral(0, h->GetSize()); 
        }
        
        // Loop through bkg samples to get names
        TIter nextBkg2(c.second.bkgList);
        object = 0;

        while ((object = nextBkg2()))
        {
            TH1D* h = (TH1D*) object;
            file << ", " << h->Integral(0, h->GetSize()); 
        }
        file << std::endl;

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
         std::cout << "========= RUN: " << vars.eventInfo.run << ", EVENT: " << vars.eventInfo.event << " =============" << std::endl;
         std::cout << std::endl;
         std::cout << "  dimuCand->mass: " << vars.dimuCand->mass << std::endl;
         std::cout << "  dimuCand->mass_PF: " << vars.dimuCand->mass_PF << std::endl;
         std::cout << "  dimuCand->pt: " << vars.dimuCand->pt << std::endl;
         std::cout << std::endl;
         std::cout << "  recoMuons->at(0): " << vars.recoMuons->at(0).outputInfo() << std::endl;
         std::cout << std::endl;
         std::cout << "  recoMuons->at(1): " << vars.recoMuons->at(1).outputInfo() << std::endl;
         std::cout << std::endl;
         std::cout << std::endl;
         std::cout << "  nJets: " << vars.jets->size() << std::endl;
         std::cout << "  nValidJets: " << vars.validJets.size() << std::endl;
         std::cout << "  nValidBJets: " << vars.validBJets.size() << std::endl;
         std::cout << std::endl;
         std::cout << "  nGenJets: " << vars.genJets.nJets << std::endl;
         std::cout << "  nValidGenJets: " << vars.validGenJets.size() << std::endl;
         std::cout << std::endl;
         std::cout << "  nValidMuons: " << vars.validMuons.size() << std::endl;
         std::cout << "  nValidExtraMuons: " << vars.validExtraMuons.size() << std::endl;
         std::cout << "  nValidElectrons: " << vars.validElectrons.size() << std::endl;
         std::cout << "  nValidExtraLeptons: " << vars.validElectrons.size() + vars.validExtraMuons.size() << std::endl;
         std::cout << std::endl;
         std::cout << "  MET: " << vars.mht->pt << std::endl;
         std::cout << std::endl;

         for(unsigned int j=0; j<vars.jets->size(); j++)
         {
             std::cout << "  jet" << j << ": " <<  vars.jets->at(j).outputInfo() << std::endl;
         }

         std::cout << std::endl;

         for(unsigned int j=0; j<vars.validJets.size(); j++)
         {
             std::cout << "  validjet" << j <<  ParticleTools::output4vecInfo(vars.validJets[j]) << std::endl;
             std::cout << "  validBjet" << j << " dR1: " << vars.validJets[j].DeltaR(vars.recoMuons->at(0).get4vec()) << std::endl;
             std::cout << "  validBjet" << j << " dR2: " << vars.validJets[j].DeltaR(vars.recoMuons->at(1).get4vec()) << std::endl;
             std::cout << std::endl;

         }
         std::cout << std::endl;

         for(unsigned int j=0; j<vars.validBJets.size(); j++)
         {
             std::cout << "  validBjet" << j <<  ParticleTools::output4vecInfo(vars.validBJets[j]) << std::endl;
             std::cout << "  validBjet" << j << " dR1: " << vars.validBJets[j].DeltaR(vars.recoMuons->at(0).get4vec()) << std::endl;
             std::cout << "  validBjet" << j << " dR2: " << vars.validBJets[j].DeltaR(vars.recoMuons->at(1).get4vec()) << std::endl;
             std::cout << std::endl;

         }
         std::cout << std::endl;

         for(unsigned int j=0; j<vars.validMuons.size(); j++)
         {
             std::cout << "  validMuon" << j <<  ParticleTools::output4vecInfo(vars.validMuons[j]) << std::endl;
         }
         std::cout << std::endl;
         for(unsigned int j=0; j<vars.validExtraMuons.size(); j++)
         {
             std::cout << "  validExtraMuon" << j <<  ParticleTools::output4vecInfo(vars.validExtraMuons[j]) << std::endl;
         }
         std::cout << std::endl;
         for(unsigned int j=0; j<vars.validElectrons.size(); j++)
         {
             std::cout << "  validEletron" << j <<  ParticleTools::output4vecInfo(vars.validElectrons[j]) << std::endl;
         }
         std::cout << std::endl;


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

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

void EventTools::cleanByDR(std::vector<TLorentzVector>& v1, std::vector<TLorentzVector>& v2, float dRmin)
{
// clean 4vecs in v1 by dR based upon 4vecs in v2, should probably move to CollectionCleaners
// Not debugged, I have no idea if this works right at the moment

    for(unsigned int i=0; i<v1.size(); i++)
    {
        for(unsigned int j=0; j<v2.size(); j++)
        {
            if(!(v1[i].DeltaR(v2[j]) > dRmin)) v1.erase(v1.begin()+i); 
            i--;
        }
    }

}
