// SignificanceMetrics.cxx

#ifndef SIGMETRICS
#define SIGMETRICS

#include "TMath.h"
#include "TH1F.h"
#include <vector>
#include <utility>
#include <cmath>

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

class SignificanceMetric
{

    public:

        // the significance is different depending on the metric, so make this abstract
        virtual double significance(double signal, double background, double unc) = 0;

        // if we have a suitable significance function, then this function is determined
        std::pair<double, double> significanceGivenCut(TH1F* signal, TH1F* background, double unc, int cutbin, bool ismin)
        {
            double sIntegral = 0;
            double bIntegral = 0;
            double cutvalue = signal->GetBinCenter(cutbin);
            // the cut is a minimum bound we want the integral from the cut to the end
            if(ismin)
            {
                sIntegral = signal->Integral(cutbin, signal->GetSize()-1);
                bIntegral = background->Integral(cutbin, signal->GetSize()-1);
            }
            // the cut is a maximum bound we want the integral from beginning to the cut
            else
            {
                sIntegral = signal->Integral(0, cutbin);
                bIntegral = background->Integral(0, cutbin);
            }
            return std::pair<double,double>(cutvalue, significance(sIntegral, bIntegral, unc) );
        }

        // if we have a suitable significance function, then this function is determined as well
        void significanceVsCut(std::vector<std::pair<double,double>>& svec, TH1F* signal, TH1F* background, double unc, bool ismin)
        {
            for(int i=0; i<signal->GetSize()-1; i++)
            {
                svec.push_back(significanceGivenCut(signal, background, unc, i, ismin)); 
            }     
        }
};

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

class AsimovSignificance : public SignificanceMetric
{
    public:

        double significance(double signal, double background, double unc)
        {
            if(unc == 0 && background == 0) return 0;
            if(background == 0 && signal == 0) return 0;
            
            double varb = background*unc*background*unc; 
            double tot = signal + background;
        
            // return the simple case for zero uncertainty
            if(unc == 0) return std::sqrt(2*(tot*std::log(1+signal/background) - signal));

            // return the full calculation when there is an uncertainty
            return std::sqrt(2*(tot*std::log((tot*(varb+background))/((background*background)+tot*varb))-(1/unc/unc)*std::log(1.+(varb*signal)/(background*(background+varb)))));
        }
};

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

class PoissonSignificance : public SignificanceMetric
{
    public:

        double significance(double signal, double background, double unc)
        {
            if(background == 0) return 0;
            return signal/TMath::Sqrt(background + unc*unc*background*background);
        }
};

#endif
