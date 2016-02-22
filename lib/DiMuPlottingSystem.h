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
        TEntryList* list;      // the list of events to associate with the sample, presumably the events that passed some set of cuts
        Sample* sample;        // the sample to plot from
        float reductionFactor; // reduce the number of events in each plot by this factor N->N/reductionFactor
                               // useful when developing code or debugging so that you don't have to take as much time to see results

        // ====================================================
        // Functions-------------------------------------------
        // ====================================================
        void initialize(Sample* sample, TEntryList* list);

        TH1F* hist1D(TString name, TString xtitle, TString title, float& varToPlot, Int_t nbins, Int_t xmin, Int_t xmax);
        TH1F* hist1D(TString name, TString xtitle, TString title, int& varToPlot, Int_t nbins, Int_t xmin, Int_t xmax);

        TH1F* weightedHist1D(TString name, TString xtitle, TString title, float& varToPlot, Int_t nbins, Int_t xmin, Int_t xmax);
        TH1F* weightedHist1D(TString name, TString xtitle, TString title, int& varToPlot, Int_t nbins, Int_t xmin, Int_t xmax);

        Double_t getScaleFactor(Float_t luminosity);
        void scaleHistByXsec(TH1F* hist, Float_t luminosity);

        void arrangeStatBox(TCanvas* c, int i);
        void arrangeLegend(TCanvas* c, int i);

        THStack* overlay(TList* list, TString title, TString xaxistitle, TString yaxistitle);
        TCanvas* stackedHistogramsAndRatio(TList* list, TString title, TString xaxistitle, TString yaxistitle);
        TH1F* addHists(TList* list, TString name);
};

#endif
