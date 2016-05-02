// DiMuPlottingSystem.h

#ifndef ADD_DIMUPLOTTINGSYSTEM
#define ADD_DIMUPLOTTINGSYSTEM

#include "TMath.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TList.h"
#include "TEntryList.h"
#include "TEventList.h"
#include "TTree.h"
#include "TCut.h"
#include "THStack.h"
#include "TPaveStats.h"

#include "Sample.h"
#include <fstream>

class DiMuPlottingSystem
{
    public:
        // ====================================================
        // Constructors/Destructors ---------------------------
        // ====================================================

        DiMuPlottingSystem();
        ~DiMuPlottingSystem();

        // ====================================================
        // Variables ------------------------------------------
        // ====================================================

        // ====================================================
        // Functions-------------------------------------------
        // ====================================================

        THStack* stackMCHistosAndData(TList* list, TString title, TString xaxistitle, TString yaxistitle);
        TCanvas* overlay(TList* list, TString name, TString title, TString xaxistitle, TString yaxistitle);
        TCanvas* stackedHistogramsAndRatio(TList* list, TString name, TString title, TString xaxistitle, TString yaxistitle);
        TH1F* addHists(TList* list, TString name);

        void arrangeStatBox(TCanvas* c, int i);
        void arrangeLegend(TCanvas* c, int i);

};

#endif
