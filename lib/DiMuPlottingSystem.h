// DiMuPlottingSystem.h

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

struct VoigtFitInfo 
{
// Datastructure to keep track of fit information for the ZCalibration

  Float_t hmean = -999;
  Float_t hmean_err = -999;
  Float_t hrms = -999;
  Float_t hrms_err = -999;

  Float_t gamma = 0.08399;

  Float_t vmean = -999;
  Float_t vmean_err = -999;
  Float_t vsigma = -999;
  Float_t vsigma_err = -999;
  Float_t vgamma = -999;
  Float_t vgamma_err = -999;

  Float_t x = -999;
  Float_t x_err = -999;

  TF1* fit = 0;
};


class ZCalibration
{
// Fit the z peak with a voigtian in different ranges of the x variable

    public:
        // ====================================================
        // Constructors/Destructors ---------------------------
        // ====================================================

        ZCalibration();
        ZCalibration(TString xname, TString massname, Float_t fitsig, Float_t massmin, 
                     Float_t massmax, Int_t massbins, Float_t xmin, Float_t xmax, Int_t xbins);
        ZCalibration(TString xname, TString massname, Float_t fitsig, Float_t massmin, 
                     Float_t massmax, Int_t massbins, std::vector<Float_t> binning);
        ~ZCalibration(){};

        // ====================================================
        // Variables ------------------------------------------
        // ====================================================
        TString xname;
        TString massname;

        bool variableBinning = false;

        Float_t xmin;
        Float_t xmax;
        Int_t xbins;

        Float_t massmin;
        Float_t massmax;
        Int_t massbins;

        Float_t fitsig;

        std::vector<Float_t> binning;
        std::vector<TH1D*> histos;
        std::vector<VoigtFitInfo> vfis;

        TGraphErrors* mean_vs_x = 0;
        TGraphErrors* resolution_vs_x = 0;

        // ====================================================
        // Functions-------------------------------------------
        // ====================================================

        void init();
        void init(std::vector<Float_t>& binning);
        void fill(Float_t xvalue, Float_t massvalue);
        void fill(Float_t xvalue, Float_t massvalue, Double_t weight);
        int whichTH1D(Float_t xvalue);
        VoigtFitInfo fit(TH1D* inhist, Float_t x, Float_t x_err);
        void fit();
        void plot();
};

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
