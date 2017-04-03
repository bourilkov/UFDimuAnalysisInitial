///////////////////////////////////////////////////////////////////////////
//                           EventTools.h                                //
//=======================================================================//
//                                                                       //
//        Miscellaneous tools: output event info to terminal             //
//        output run,event from vector to CSV, output vars and values    //
//        from map to CSV, etc                                           //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef EVENTTOOLS
#define EVENTTOOLS

#include <vector>
#include <utility> // pair
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>

#include "TString.h"
#include "CategorySelection.h"
#include "JetCollectionCleaner.h"

class EventTools
{
    public: 
        EventTools(){};
        ~EventTools(){};

        static void loadEventsFromFile(TString filename, std::vector<std::pair<int, long long int>>& v);
        static void outputEventsToFile(std::vector<std::pair<int,long long int>>& v, TString filename);
        static TString outputMapKeysCSV(std::map<TString,double>& map);
        static TString outputMapValuesCSV(std::map<TString,double>& map);
        static bool sameRunAndEvent(std::pair<int,long long int> a, std::pair<int,long long int> b);
        static bool eventInVector(std::pair<int,long long int> e, std::vector<std::pair<int,long long int>> events);
        static void outputEvent(VarSet& vars);
        static void outputEvent(VarSet& vars, Categorizer& categorizer);
};
#endif
