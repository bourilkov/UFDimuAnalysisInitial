/////////////////////////////////////////////////////////////////////////////
//                           listXMLNodes.cxx                              //
//=========================================================================//
//                                                                         //
//       List the nodes created by our decision tree autocategorizer.      //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

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
#include <unordered_map>

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    TString xmlfile; 

    for(int i=1; i<argc; i++) 
    {    
        std::stringstream ss;  
        ss << argv[i];
        if(i==1) xmlfile = TString(ss.str().c_str());
    }    

   CategorySelectionBDT xmlc(xmlfile);
   xmlc.outputResults();

}
