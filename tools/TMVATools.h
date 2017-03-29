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

        static void getVarNamesRecursive(TXMLEngine* xml, XMLNodePointer_t node, std::vector<TString>& tvars, std::vector<TString>& svars);
        static void getVarNames(TString filename, std::vector<TString>& tvars, std::vector<TString>& svars);
        
        static TMVA::Reader* bookVars(TString methodName, TString weightfile, std::map<TString, Float_t>& tmap, std::map<TString, Float_t>& smap);
        static float getClassifierScore(TMVA::Reader* reader, TString methodName, std::map<TString, Float_t>& tmap, VarSet& varset);
};
#endif
