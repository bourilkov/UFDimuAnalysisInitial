///////////////////////////////////////////////////////////////////////////
//                           TMVATools.h                                //
//=======================================================================//
//                                                                       //
//        Miscellaneous tools to help with TMVA Classification.          //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef TMVATOOLS
#define TMVATOOLS

#include <vector>
#include <map>
#include <iostream>

#include "TXMLEngine.h"
#include "TMVA/Reader.h"
#include "TString.h"
#include "VarSet.h"

class TMVATools
{
    public: 
        TMVATools(){};
        ~TMVATools(){};

        static void getNames(TString filename, TString node_name, TString node_att, std::vector<TString>& names);
        static void getNamesRecursive(TXMLEngine* xml, XMLNodePointer_t node, TString node_name, TString node_att, std::vector<TString>& names);
        static void getVarNames(TString filename, std::vector<TString>& tvars, std::vector<TString>& svars);
        static void getClassNames(TString filename, std::vector<TString>& classes);
        
        static TMVA::Reader* bookVars(TString methodName, TString weightfile, std::map<TString, Float_t>& tmap, std::map<TString, Float_t>& smap);
        static float getClassifierScore(TMVA::Reader* reader, TString methodName, std::map<TString, Float_t>& tmap, VarSet& varset);
        static std::vector<float> getMulticlassScores(TMVA::Reader* reader, TString methodName, std::map<TString, Float_t>& tmap, VarSet& varset);
};
#endif
