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
        DiMuPlottingSystem(Sample* sample, TEntryList* list);
        ~DiMuPlottingSystem();

        // ====================================================
        // Variables ------------------------------------------
        // ====================================================
        TEntryList* list;
        Sample* sample;

        // ====================================================
        // Functions-------------------------------------------
        // ====================================================
        void initialize(Sample* sample, TEntryList* list);

        TH1F* hist1D(TString name, TString xtitle, TString title, float& varToPlot, Int_t nbins, Int_t xmin, Int_t xmax);
        TH1F* hist1D(TString name, TString xtitle, TString title, int& varToPlot, Int_t nbins, Int_t xmin, Int_t xmax);
  
        TH1F* genWeightedHist1D(TString name, TString xtitle, TString title, float& varToPlot, Int_t nbins, Int_t xmin, Int_t xmax);
        TH1F* genWeightedHist1D(TString name, TString xtitle, TString title, int& varToPlot, Int_t nbins, Int_t xmin, Int_t xmax);

        TH1F* gPUWeightedHist1D(TString name, TString xtitle, TString title, float& varToPlot, Int_t nbins, Int_t xmin, Int_t xmax);
        TH1F* gPUWeightedHist1D(TString name, TString xtitle, TString title, int& varToPlot, Int_t nbins, Int_t xmin, Int_t xmax);

        Double_t scaleByXsec(Float_t luminosity);
        void scaleHistByXsec(TH1F* hist, Float_t luminosity);

        void arrangeStatBox(TCanvas* c, int i);
        void arrangeLegend(TCanvas* c, int i);

        THStack* overlay(TList* list, TString title, TString xaxistitle, TString yaxistitle);
        TCanvas* stackedHistogramsAndRatio(TList* list, TString title, TString xaxistitle, TString yaxistitle);
        TH1F* addHists(TList* list, TString name);
};

#endif
