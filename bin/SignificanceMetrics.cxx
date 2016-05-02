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

        double unc_0 = 0;  // the percent uncertainty for a given number of backbround events nbg_0
        double nbg_0 = 0;  // the number of background events when the percent uncertainty unc_0 was calculated
        double unc = 0;    // the percent uncertainty for the current number of bg events

        SignificanceMetric(double nbg_0, double unc_0)
        {
            this->nbg_0 = nbg_0;
            this->unc_0 = unc_0;
        }

        void setUncertainty(double background)
        {
        // The uncertainty appears to scale with the sqrt(B).
        // we use the uncertainty and nbackground from the run1 paper or some other benchmark and then scale from there
        // depending on how much background there is in the current calculation
            unc = unc_0 * TMath::Sqrt(nbg_0/background);
        }


        // the significance is different depending on the metric, so make this abstract
        virtual double significance(double signal, double background) = 0;

        // if we have a suitable significance function, then this function is determined
        std::pair<double, double> significanceGivenCut(TH1F* signal, TH1F* background, int cutbin, bool ismin)
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
            return std::pair<double,double>(cutvalue, significance(sIntegral, bIntegral) );
        }

        // if we have a suitable significance function, then this function is determined as well
        void significanceVsCut(std::vector<std::pair<double,double>>& svec, TH1F* signal, TH1F* background, bool ismin)
        {
            for(int i=0; i<signal->GetSize()-1; i++)
            {
                svec.push_back(significanceGivenCut(signal, background, i, ismin)); 
            }     
        }

        static void outputSignificanceVsCut(std::vector<std::pair<double,double>>& svec)
        {
            for(unsigned int i=0; i<svec.size(); i++)
            {   
                std::cout << svec[i].first << ", " << svec[i].second << std::endl;
            }   
            std::cout << std::endl;
        }

        static TGraph* makeTGraph(std::vector<std::pair<double,double>>& svec, TString name, TString title, TString xtitle, TString ytitle)
        {
            TGraph* graph = new TGraph();
            graph->SetName(name);
            graph->SetTitle(title);

            for(unsigned int i=0; i<svec.size(); i++)
            {   
                graph->SetPoint(i, svec[i].first, svec[i].second);
            }   

            graph->GetXaxis()->SetTitle(xtitle);
            graph->GetYaxis()->SetTitle(ytitle);

            return graph;
        }
};

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

class AsimovSignificance : public SignificanceMetric
{
    public:

        AsimovSignificance(double nbg_0, double unc_0) : SignificanceMetric(nbg_0, unc_0){}

        double significance(double signal, double background)
        {
            setUncertainty(background);

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

        PoissonSignificance(double nbg_0, double unc_0) : SignificanceMetric(nbg_0, unc_0){}

        double significance(double signal, double background)
        {
            setUncertainty(background);
            if(background == 0) return 0;
            return signal/TMath::Sqrt(background + unc*unc*background*background);
        }
};

#endif
