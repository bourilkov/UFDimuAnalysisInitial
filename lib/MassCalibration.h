///////////////////////////////////////////////////////////////////////////
// ======================================================================//
// MassCalibration.h                                                        //
// ======================================================================//
// Plot Z Mass Resolution vs mu+/- eta, mu+/- phi, etc.                  //
// Creates mass histograms in different bins of some x variable.         //
// The object can then fit the histograms with a Voigtian and plot       //
// the voigt fit mean/resolution vs x.                                   //
// ======================================================================//
///////////////////////////////////////////////////////////////////////////

#ifndef ADD_ZCALIBRATION
#define ADD_ZCALIBRATION

#include "TMath.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TList.h"
#include "TTree.h"
#include "TPaveStats.h"

#include "Sample.h"
#include <fstream>

///////////////////////////////////////////////////////////////////////////
// -----------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////

// Datastructure to keep track of the histogram mean and rms as well as
// the fit mean and rms for each histogram, one for each xbin
struct VoigtFitInfo 
{
// Datastructure to keep track of fit information for the MassCalibration

  Float_t hmean = -999;
  Float_t hmean_err = -999;
  Float_t hrms = -999;
  Float_t hrms_err = -999;

  Float_t gamma = 0.0838;

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

///////////////////////////////////////////////////////////////////////////
// -----------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////

class MassCalibration
{
// Fit the z peak with a voigtian in different ranges of the x variable

    public:
        // ====================================================
        // Constructors/Destructors ---------------------------
        // ====================================================

        MassCalibration();

        // For a nonvariable binning scheme, all bins with the same width
        MassCalibration(TString xname, TString massname, Float_t fitsig, Float_t massmin, 
                     Float_t massmax, Int_t massbins, Float_t xmin, Float_t xmax, Int_t xbins, Float_t voigt_gamma = 2.5);

        // For a variable binning scheme
        MassCalibration(TString xname, TString massname, Float_t fitsig, Float_t massmin, 
                     Float_t massmax, Int_t massbins, std::vector<Float_t> binning, Float_t voigt_gamma = 2.5);
        ~MassCalibration(){};

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

        // The number of sigma-withs on the histogram to fit the voigtian
        // fitsig=1 means fit the voigtian to a 1 sigma fit window around the mean
        Float_t fitsig;
        Float_t voigt_gamma = 2.5;

        std::vector<Float_t> binning;
        std::vector<TH1D*> histos;
        std::vector<VoigtFitInfo> vfis;

        TGraphErrors* mean_vs_x = 0;
        TGraphErrors* resolution_vs_x = 0;

        // ====================================================
        // Functions-------------------------------------------
        // ====================================================

        // Create and initialize histograms for every x bin
        void init();
        void init(std::vector<Float_t>& binning);

        // Fill the correct x-bin histogram with the massvalue
        void fill(Float_t xvalue, Float_t massvalue);
        void fill(Float_t xvalue, Float_t massvalue, Double_t weight);
        int whichTH1D(Float_t xvalue);

        // Fit a single histogram with a voigtian
        VoigtFitInfo fit(TH1D* inhist, Float_t x, Float_t x_err);

        // Fit all of the histograms with voigtians
        // Save fit info (mean/resolution) to vfis
        void fit();
       
        // create mean_vs_x and resolution_vs_x TGraphs 
        void plot();
};

#endif
