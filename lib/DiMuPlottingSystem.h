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
        ZCalibration(TString xname, Float_t fitsig, Float_t massmin, Float_t massmax, Int_t massbins, Float_t xmin, Float_t xmax, Int_t xbins);
        ~ZCalibration();

        // ====================================================
        // Variables ------------------------------------------
        // ====================================================
        TString xname;

        Float_t xmin;
        Float_t xmax;
        Int_t xbins;

        Float_t massmin;
        Float_t massmax;
        Int_t massbins;

        Float_t fitsig;

        std::vector<Float_t> binning;
        std::vector<TH1F*> histos;
        std::vector<VoigtFitInfo> vfis;
        TGraphErrors* massvsx;

        // ====================================================
        // Functions-------------------------------------------
        // ====================================================

        void init();
        void fill(Float_t xvalue, Float_t massvalue);
        int whichTH1F(Float_t xvalue);
        VoigtFitInfo fit(TH1F* inhist, Float_t x, Float_t x_err);
        void fit();
        TGraphErrors* plot();
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

        THStack* stackMCHistosAndData(TList* list, TString title, TString xaxistitle, TString yaxistitle);
        TCanvas* overlay(TList* list, TString name, TString title, TString xaxistitle, TString yaxistitle);
        TCanvas* stackedHistogramsAndRatio(TList* list, TString name, TString title, TString xaxistitle, TString yaxistitle, TString ratiotitle = "Data/MC");
        TH1F* addHists(TList* list, TString name);

        void arrangeStatBox(TCanvas* c, int i);
        void arrangeLegend(TCanvas* c, int i);

};

#endif
