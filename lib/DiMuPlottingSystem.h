///////////////////////////////////////////////////////////////////////////
// ======================================================================//
// DiMuPlottingSystem.h                                                  //
// ======================================================================//
// Plotting utilities for the analysis.                                  //
// Compare MC stack to Data w/ ratio plot, add histograms together,      //
// overlay tgraphs, rebin histos based on rato plot error, etc           //
// ======================================================================//
///////////////////////////////////////////////////////////////////////////

#ifndef ADD_DIMUPLOTTINGSYSTEM
#define ADD_DIMUPLOTTINGSYSTEM

#include "TMath.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1D.h"
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

        static THStack* stackComparison(TList* list, TString title, TString xaxistitle, TString yaxistitle, bool log = true, bool stats = false, int legend = 0);
        static TCanvas* overlay(TList* list, TString name, TString title, TString xaxistitle, TString yaxistitle, bool log = true);

        static float ratioError2(float numerator, float numeratorError2, float denominator, float denominatorError2);
        static void  getBinningForRatio(TH1D* numerator, TH1D* denominator, std::vector<Double_t>& newBins, float maxPercentError = 0.1);

        static TCanvas* stackedHistogramsAndRatio(TList* list, TString name, TString title, TString xaxistitle, TString yaxistitle, bool rebin = false, bool fit = true, 
                                           TString ratiotitle = "Data/MC", bool log = true, bool stats = false, int legend = 0);
        static TH1D* addHists(TList* list, TString name, TString title);

        static void arrangeStatBox(TCanvas* c, int i);
        static void arrangeLegend(TCanvas* c, int i);
};

#endif
