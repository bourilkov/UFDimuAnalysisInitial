// Sample.h

#ifndef ADD_SAMPLE
#define ADD_SAMPLE

#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TEntryList.h"
#include "TROOT.h"

#include "LumiReweightingStandAlone.h"
#include "VarSet.h"
#include "BranchSet.h"

class Sample
{
    public:
        Sample();
	Sample(std::vector<TString> infilenames, TString name, TString sampleType="NONE");
        ~Sample();

        TString name;
        TString filename;
	std::vector<TString> filenames;
        TString treename;
	TChain* chain;
        reweight::LumiReWeighting* lumiWeights;  // Information for pileup reweighting 

        TString dir;           // DAS directory
        TString pileupfile;    // used for pile up reweighting through lumiWeights
        TString sampleType;    // "data", "signal", "background"

        int plotColor;         // the color used when plotting the sample
        int nOriginal;         // the number of events run over to get this sample
        int nOriginalWeighted; // the number of original events run over to get this sample accounting for genWeights
        int N;                 // the number of events in the sample

        float xsec;            // xsec in pb
        float lumi;            // the luminosity of the data or effective luminosity of MC

        VarSet vars;           // all of the variables from the ttree
        BranchSet branches;    // all of the branches from the ttree, which we link to the vars

        int getEntry(int i);                    // load the ith event from the ttree into vars
        int getEntry(int i, TEntryList* list);  // load ith event from the list into vars
                                                // the ith event in the list maps to the jth tree entry

        void calculateNoriginal();                    // calculate nOriginal and nOriginalWeighted
        void setBranchAddresses(int whichCategories=1); // link the values in the tree to vars
        double getWeight();         // get the weight for the histogram based upon the pileup weight and the MC gen weight

        // get the scale factor for the MC histogram based upon the number of events, the data luminosity, and the xsec for the process 
        float getLumiScaleFactor(float luminosity); 

    protected:
	std::vector<TFile*> files; // files with the ttree

};

#endif
