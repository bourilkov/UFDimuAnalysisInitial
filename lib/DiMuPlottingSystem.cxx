///////////////////////////////////////////////////////////////////////////
//                             DiMuPlottingSystem.cxx                  //
//=======================================================================//
//                                                                       //
//        This class draws all the plots we will need for analysis.      //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "DiMuPlottingSystem.h"

#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TList.h"
#include "TProfile.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLatex.h"

#include <vector>
#include <sstream>
#include <cmath>
#include <iostream>


///////////////////////////////////////////////////////////////////////////
// _______________________Constructor/Destructor_________________________//
///////////////////////////////////////////////////////////////////////////

DiMuPlottingSystem::DiMuPlottingSystem() {
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

DiMuPlottingSystem::~DiMuPlottingSystem() {
}

///////////////////////////////////////////////////////////////////////////
// _______________________Other Functions________________________________//
///////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::arrangeStatBox(TCanvas* c, int i)
{
  // Align the stat box better.
  // Input i determines how big the stat box should be and where it should go

  c->Update(); // Need this line to get the stats box
  TPaveStats *st = (TPaveStats*)c->GetPrimitive("stats");
  st->SetX1NDC(0.68);
  st->SetY1NDC(0.64);
  st->SetX2NDC(0.88);
  st->SetY2NDC(0.88);
  st->Draw("same");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::arrangeLegend(TCanvas* c, int i)
{
  // Align the legend better.
  // Input i determines how big the legend should be and where it should go

  // Might need to get TLegend from the pad instead of the canvas

  c->Update(); // Need this line to get the stats box
  TPaveStats *st = (TPaveStats*)c->GetPrimitive("stats");
  st->SetX1NDC(0.68);
  st->SetY1NDC(0.64);
  st->SetX2NDC(0.88);
  st->SetY2NDC(0.88);
  st->Draw("same");
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

TCanvas* DiMuPlottingSystem::overlay(TList* ilist, TString name, TString title, TString xaxistitle, TString yaxistitle, bool log)
{
// Overlay itmes in the list. 


// Wide bottom left
//  TLegend* l = new TLegend(1, 0.15, 0.7, 1.25, "", "brNDC");
 
  // Square-ish top right
  //TLegend* l = new TLegend(0.68, 0.56, 0.88, 0.87, "", "brNDC");

  // Square-ish top left
  //TLegend* l = new TLegend(0.13, 0.56, 0.33, 0.88, "", "brNDC");
  
  // Wide rectangle top right
  TLegend* l = new TLegend(0.846, 0.854, 1.0, 1.0, "", "brNDC");
  TCanvas* c = new TCanvas(name);
  c->SetTitle(title);
  c->SetGridx(1);
  c->SetGridy(1);

  TIter next(ilist);
  TObject* object = 0;
  int i=0;

  std::vector<int> colors = {1, 2, 4, 7, 8, 36, 50, 30, 9, 29, 3, 42, 98, 62, 74, 20, 29, 32, 49, 12, 3, 91};

  float max = -9999;
  float min = 9999;

  TMultiGraph* multigraph = 0;
  THStack* stack = 0;
  float minimum = 999999999;  // minimum of all the histograms
  float minMax = 999999999;   // minimum of the largest values of the histograms

  while ((object = next()))
  {
      //TH1D* hist = (TH1D*) object;
      TGraph* graph = 0;
      TH1D* hist = 0;

      if(object->InheritsFrom("TH1D"))
      {
          hist = (TH1D*) object;
          hist->SetLineColor(colors[i]);
          hist->SetLineWidth(2);
          //hist->SetFillColor(colors[i]);
          
          std::cout << hist->GetName() << std::endl;

          if(object == ilist->First())
          {
              minimum = hist->GetMinimum();
              minMax = hist->GetMaximum();
              stack = new THStack();
              stack->SetTitle(title+"_"+xaxistitle+"_"+yaxistitle);
          }

          if(hist->GetMinimum() < minimum) minimum = hist->GetMinimum();
          if(hist->GetMaximum() < minMax) minMax = hist->GetMaximum();

          stack->Add(hist);
          
          TString legend_entry = TString(hist->GetTitle());
          l->AddEntry(hist, legend_entry, "f");

          if(object == ilist->Last())
          {   
               //stack->Draw();
              stack->Draw("nostack,l");
              stack->GetXaxis()->SetTitle(xaxistitle);
              stack->GetYaxis()->SetTitle(yaxistitle);
              gPad->Modified();
              if(log) 
              {
                  gPad->SetLogy(1);
                  stack->SetMaximum(stack->GetMaximum()*1.5);
                  if(minimum > 0) stack->SetMinimum(minimum*0.5);
                  else stack->SetMinimum(minMax*10e-3);
              }
          }
      }
      if(object->InheritsFrom("TGraph"))
      {
          if(object == ilist->First())
          {
              multigraph = new TMultiGraph();
              multigraph->SetTitle(title);
          }

          graph = (TGraph*) object;
          graph->SetLineColor(colors[i]);
          graph->SetLineWidth(3);

          multigraph->Add(graph);

          TString legend_entry = TString(graph->GetTitle());
          l->AddEntry(graph, legend_entry, "l");

          if(object == ilist->Last())
          {
              multigraph->Draw("a");
              multigraph->GetXaxis()->SetTitle(xaxistitle);
              multigraph->GetYaxis()->SetTitle(yaxistitle);
          }
      }

      i++;
  }
    
  l->Draw();
  return c;
}
//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

THStack* DiMuPlottingSystem::stackComparison(TList* ilist, TString title, TString xaxistitle, TString yaxistitle, bool log, bool stats, int legend)
{
// Creates a THStack of the histograms in the list. Overlays the data ontop of this without adding it to the stack.
// Assumes data is in the last ilist location so that it appears on top of the other histograms.


// Wide bottom left
//  TLegend* l = new TLegend(1, 0.15, 0.7, 1.25, "", "brNDC");
 
  // Square-ish top right
  //TLegend* l = new TLegend(0.68, 0.56, 0.88, 0.87, "", "brNDC");

  // Wide rectangle top right
  TLegend* l = 0;
  if(legend == 0) l = new TLegend(0.71, 0.82, 1.0, 1.0, "", "brNDC");
  else if(legend == 1) l = new TLegend(0.313, 0.491, 0.883, 0.886, "", "brNDC");
  else l = new TLegend(0.71, 0.82, 1.0, 1.0, "", "brNDC"); 

  THStack* stack = new THStack();
  stack->SetTitle(title);

  // Square-ish top left
  //TLegend* l = new TLegend(0.13, 0.56, 0.33, 0.88, "", "brNDC");

  TIter next(ilist);
  TObject* object = 0;
  int i=0;
  float minimum = 999999999;
  float minMax = 999999999;

  std::vector<int> colors = {58, 2, 8, 36, 91, 46, 50, 30, 9, 29, 3, 42, 98, 62, 74, 20, 29, 32, 49, 12, 3, 91};

  while ((object = next()))
  {
      TH1D* hist = (TH1D*) object;
      //std::cout << Form("%d: Adding %s to the stack \n", i, hist->GetName());
      hist->SetStats(0);
      hist->SetFillColor(colors[i]);
      hist->SetLineColor(colors[i]);

      // Print name + num events in legend
      TString legend_entry = TString(hist->GetTitle());
      //legend_entry.Form("%s %d", v[i]->GetName(), (int)v[i]->GetEntries());
      //std::cout << v[i]->GetName() << ": " << v[i]->GetEntries() << std::endl;
      //std::cout << legend_entry << std::endl;
      //if((i+1)==10) v[i]->SetFillColor(50);
      //v[i]->SetLineWidth(3);
 
      if(object == ilist->First()) 
      {
          minimum = hist->GetMinimum();
          minMax = hist->GetMaximum();
      }

      if(hist->GetMinimum() < minimum) minimum = hist->GetMinimum();
      if(hist->GetMaximum() < minMax) minMax = hist->GetMaximum();

      // Assuming data is in the last location
      if(object == ilist->Last())
      {
          hist->SetMarkerStyle(20);
          hist->SetLineColor(1);
          hist->SetFillColor(0);

          stack->Draw("hist");
          stack->GetXaxis()->SetTitle(xaxistitle);
          stack->GetYaxis()->SetTitle(yaxistitle);
          
          // overlay the mc error band
          TH1D* sum = (TH1D*) stack->GetStack()->Last()->Clone();
          sum->SetFillStyle(3001);
          sum->SetLineColor(12);
          sum->SetFillColor(12);
          sum->Draw("E2SAME");
          l->AddEntry(sum, "MC Error Band", "f");

          gPad->Modified();

          // use y log scale if the flag is set
          if(log) 
          {
              gPad->SetLogy(1);
              stack->SetMaximum(stack->GetMaximum()*10);
              if(minimum > 0) stack->SetMinimum(minimum*0.1);
              else stack->SetMinimum(minMax*10e-3);
          }

          // only draw the stat box if the flag is set
          if(stats) hist->SetStats(1);
          //hist->SetMinimum(10e-4);
	  hist->Draw("epsames");

          // put the stat box in a decent location
          if(stats)
          {
              gPad->Update();
              TPaveStats *st=(TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
              // legend == 0 default alignment for stat box
              st->SetX1NDC(0.52);
              st->SetY1NDC(0.83);
              st->SetX2NDC(0.70);
              st->SetY2NDC(0.93);

              // different legend setting alignments not yet implemented
              gPad->Modified();
          }
          l->AddEntry(hist, legend_entry, "p");
      }
      else
      {
        stack->Add(hist);
        l->AddEntry(hist, legend_entry, "f");
      }

      i++;
  }
    
  l->Draw();

  return stack;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

TH1D* DiMuPlottingSystem::addHists(TList* ilist, TString name, TString title)
{
// Add all of the histograms in the list into a single histogram

    TH1D* htotal = (TH1D*)ilist->First()->Clone(name);
    htotal->SetTitle(title);

    TIter next(ilist);
    TObject* object = 0;

    while ((object = next()))
    {
        TH1D* h = (TH1D*) object;
        if(object != ilist->First()) htotal->Add(h);
    }
    return htotal;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

float DiMuPlottingSystem::ratioError2(float numerator, float numeratorError2, float denominator, float denominatorError2)
{
// Propogation of errors for N_numerator/N_denominator
// Returns error squared, takes in numerator error squared and denom err squared as well as N_num and N_denom
  if (numerator   <= 0) return 0;
  if (denominator <= 0) return 0;

  const float ratio = numerator/denominator;

  return  ratio*ratio*( numeratorError2/numerator/numerator + denominatorError2/denominator/denominator );
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::getBinningForRatio(TH1D* numerator, TH1D* denominator, std::vector<Double_t>& newBins, float maxPercentError)
{
// The ratio plots are a bit crazy with huge errors sometimes, so we want to rebin with variable binning
// such that the error is always low in each of the ratio plot bins

    TString numName = numerator->GetName();
    bool isMassData = false;
    if(numName.Contains("dimu_mass") && numName.Contains("Data")) isMassData = true;

    // check if the two histos are binned the same way
    bool compatible = numerator->GetNbinsX() == denominator->GetNbinsX() && numerator->GetBinLowEdge(1) == denominator->GetBinLowEdge(1)
                      && numerator->GetBinLowEdge(numerator->GetNbinsX()) == denominator->GetBinLowEdge(denominator->GetNbinsX());
    if(!compatible)
    {
        std::cout << "!!! ERROR: DiMuPlottingSystem::getBinningForRatio numerator and denominator histograms are not binned the same!" << std::endl;
        return;
    }

    if(numerator->Integral() <= 0 || denominator->Integral() <=0)
    {
        std::cout << "!!! ERROR: DiMuPlottingSystem::getBinningForRatio numerator or denominator integral <= 0" << std::endl;
        return;
    }

    std::vector<double> errVec = std::vector<double>();

    // first low edge is the lowest edge by default
    newBins.push_back(numerator->GetBinLowEdge(1));
    double sumNum = 0;
    double sumErr2Num = 0;
    double sumDenom = 0;
    double sumErr2Denom = 0;
    bool lastBinIsCutBin = 0;

    // Collect bins together until the error is low enough for the corresponding ratio plot bin.
    // Once the error is low enough call that the new bin, move on and repeat.
    for(unsigned int i=1; i<=numerator->GetNbinsX(); i++)
    {
        sumNum += numerator->GetBinContent(i);
        sumErr2Num += numerator->GetBinError(i)*numerator->GetBinError(i);

        sumDenom += denominator->GetBinContent(i);
        sumErr2Denom += denominator->GetBinError(i)*denominator->GetBinError(i);

        // squared ratio plot error
        // 0 if sumNum or sumDenom is 0
        float ratioErr2 = ratioError2(sumNum, sumErr2Num, sumDenom, sumErr2Denom); 

        float percentError = (ratioErr2>0)?TMath::Sqrt(ratioErr2)/(sumNum/sumDenom):0;

        // If the bins suddenly drop to zero or go up from zero then this is probably due to some cut that was applied
        // and we dont' want to combine the 0 bins on one side of the cut with the bins after the cut. It messes up the comparison plots.
        // it's better to see exactly where the cut caused everything to drop to zero.
        if(i!=numerator->GetNbinsX() && TMath::Abs(numerator->GetBinContent(i)-numerator->GetBinContent(i+1)) > 10 
          && (numerator->GetBinContent(i) == 0 || numerator->GetBinContent(i+1) == 0))
        {
            newBins.push_back(numerator->GetBinLowEdge(i)+numerator->GetBinWidth(i));
            errVec.push_back(percentError);

            sumNum = 0;
            sumErr2Num = 0;

            sumDenom = 0;
            sumErr2Denom = 0;
 
            lastBinIsCutBin = true;
        }
        // we have the minimum error necessary create the bin in the vector
        // or just make a bin if we are in the blinded region of the mass spectrum
        else if(percentError !=0 && percentError <= maxPercentError || (isMassData && numerator->GetBinLowEdge(i) >= 110 && numerator->GetBinLowEdge(i) < 140))
        {
            newBins.push_back(numerator->GetBinLowEdge(i)+numerator->GetBinWidth(i));
            errVec.push_back(percentError);

            sumNum = 0;
            sumErr2Num = 0;

            sumDenom = 0;
            sumErr2Denom = 0;

            lastBinIsCutBin = false;
        }
        // we got to the end of the histogram and the last bins can't add up
        // to get the error low enough, so we merge these last bins with the 
        // last bin set up in the new variable binning scheme.
        if(i==numerator->GetNbinsX() && (percentError > maxPercentError || percentError == 0))
        {
            // get rid of the last bin added, but make sure there is at least one bin. Don't remove the low edge of the zero-th bin.
            // la la we need at least a value for the beginning and end, size 1 vector doesn't make sense and rebinning will fail.
            if(newBins.size() > 1 && !lastBinIsCutBin) newBins.pop_back(); 
            // replace with the end bin value in the numerator histo
            newBins.push_back(numerator->GetBinLowEdge(i)+numerator->GetBinWidth(i));
        }
    }
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

TCanvas* DiMuPlottingSystem::stackedHistogramsAndRatio(TList* ilist, TString name, TString title, TString xaxistitle, TString yaxistitle, bool rebin, 
                                                       bool fit, TString ratiotitle, bool log, bool stats, int legend)
{
  
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // first let's see if the numerator and denom are valid to see if this will work out
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    // data usually
    TH1D* numerator = (TH1D*) ilist->Last();
    // The rest of the list is usually the MC used to make the stack to compare to data
    TList* denomlist = (TList*) ilist->Clone("blah");
    // get rid of the data histogram from the list
    denomlist->RemoveLast();
    // Add up the MC histograms
    TH1D* denominator = addHists(denomlist, name+"_add", name+"_add");

    // Don't worry about drawing the stack. Just return a reasonable TCanvas.
    // If we were to return the canvas after creating the stack and ratio with bad info then TCanvas::SaveAs() gives a seg fault
    if(numerator->Integral() <= 0 || denominator->Integral() <= 0)
    {
        std::cout << "!!! " << std::endl;
        std::cout << "!!! ERROR: DiMuPlottingSystem::stackedHistogramsAndRatio numerator or denominator for stack has <= 0 integral." << std::endl;
        std::cout << "!!! Returning canvas with only the numerator drawn for "<< name << "." << std::endl;
        std::cout << "!!! " << std::endl;
        std::cout << std::endl;

        TCanvas* c = new TCanvas(name);
        c->SetTitle(title);
        numerator->Draw();
        return c;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // OK the numerator and denominator are fine, let's move on to making the stack and ratio plot
    ///////////////////////////////////////////////////////////////////////////////////////////////

    // rebin the histograms by combining bins so that the error in each ratio plot bin is low enough
    // Leave the histograms from the input list alone, and make the stack the same way
    // but make the ratio plot with the rebinned numerator and denominator. 
    std::vector<Double_t> newBins;
    if(rebin)
    {
        // now put the variable binning into newBins
        getBinningForRatio(numerator, denominator, newBins, 0.20);

        // use the rebinnedList for the stack and ratio plot
        if(newBins.size() > 2)
        {
            numerator = (TH1D*) numerator->Rebin(newBins.size()-1, numerator->GetName()+TString("_"), &newBins[0]);
            denominator = (TH1D*) denominator->Rebin(newBins.size()-1, numerator->GetName()+TString("_"), &newBins[0]);
        }
        else
        {
            std::cout << "!!! " << std::endl;
            std::cout << "!!! ERROR: DiMuPlottingSystem::stackedHistogramsAndRatio rebinning vector has <=2 entries. Cannot rebin. Using original binning." << std::endl;
            std::cout << "!!! " << std::endl;
        }
        std::cout << "newBins for " << name << std::endl;
        for(unsigned int i=0; i< newBins.size(); i++)
        {
            std::cout << newBins[i] << ", ";
        }
        std::cout << std::endl << std::endl;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Put the stack in the upper pad of the canvas
    ///////////////////////////////////////////////////////////////////////////////////////////////

    // Define the Canvas
    TCanvas *c = new TCanvas(name);
    c->SetTitle(title);
    c->SetCanvasSize(800,800);

    // Need access to this to set the axes and titles and whatnot
    // also need access to the last one, but we already have numerator = ilist->last
    TH1D* first = (TH1D*) ilist->At(0);

    Double_t scale = 1;

     // Upper pad
    TPad *pad1 = new TPad(TString("pad1")+first->GetName(), "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0); 
    //pad1->SetGridx();         // Vertical grid
    //pad1->SetGridy();         // horizontal grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad

    // Draw stacked histograms here
    first->SetStats(0);          // No statistics on upper plot
    THStack* stack = stackComparison(ilist, title, xaxistitle, yaxistitle, log, stats, legend);

    // Not sure how to work with this
    // Do not draw the Y axis label on the upper plot and redraw a small
    // axis instead, in order to avoid the first label (0) to be clipped.
    //v[0]->GetYaxis()->SetLabelSize(0.);
    //TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
    //axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    //axis->SetLabelSize(15);
    //axis->Draw();

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Put the ratio plot in the lower pad of the canvas
    ///////////////////////////////////////////////////////////////////////////////////////////////

    // Lower pad
    c->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad(TString("pad2_")+name, "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy(); // horizontal grid
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad

    // Define the ratio plot
    // Clone the data histogram from the last location in the vector
    // The remainder of the vector consits of MC samples
    TH1D* hratio = (TH1D*)numerator->Clone("hratio");

    // Output the overall scale discrepancy between the stack and the hist of interest
    scale = hratio->Integral()/denominator->Integral();
    std::cout << std::endl;
    std::cout << "########## Numerator integral for " << numerator->GetName() << "  : " << numerator->Integral() << std::endl;
    std::cout << "########## Denominator integral for " << denominator->GetName() << ": " << denominator->Integral() << std::endl;
    std::cout << "$$$$$$$$$$ Scale factor for MC stack: " << scale << std::endl;
    //hadd->Scale(scale);

    hratio->SetLineColor(kBlack);
    hratio->SetMinimum(0.58);  // Define Y ..
    hratio->SetMaximum(1.42); // .. range
    hratio->Sumw2();
    hratio->SetStats(1);

    // make the ratio plot
    hratio->Divide(denominator);

    // continue setting up the plot to draw
    hratio->SetMarkerStyle(20);
    hratio->Draw("ep");       // Draw the ratio plot

    // I don't know what this is doing ...
    //last->GetYaxis()->SetTitleSize(20);
    //last->SetTitleFont(43);
    //last->GetYaxis()->SetTitleOffset(1.55);

    // Ratio plot (hratio) settings
    hratio->SetTitle(""); // Remove the ratio title

    // Y axis ratio plot settings
    hratio->GetYaxis()->SetTitle(ratiotitle);
    hratio->GetYaxis()->SetNdivisions(505);
    hratio->GetYaxis()->SetTitleSize(20);
    hratio->GetYaxis()->SetTitleFont(43);
    hratio->GetYaxis()->SetTitleOffset(1.55);
    hratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hratio->GetYaxis()->SetLabelSize(15);

    // X axis ratio plot settings
    hratio->GetXaxis()->SetTitle(xaxistitle);
    hratio->GetXaxis()->SetTitleSize(20);
    hratio->GetXaxis()->SetTitleFont(43);
    hratio->GetXaxis()->SetTitleOffset(4.);
    hratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hratio->GetXaxis()->SetLabelSize(15);

    if(!fit) hratio->SetStats(0);

    if(fit && hratio->Integral() != 0)
    {
        // Display fit info on canvas, mess with stat box styles.
        //gStyle->SetStatStyle(0);
        //gStyle->SetStatBorderSize(0);
        gStyle->SetOptFit(0111);

        // Fit things.
        TString fitname = TString("linear_fit_")+name;
        TF1* fit = new TF1(fitname, "[1]*x+[0]", hratio->GetXaxis()->GetXmin(), hratio->GetXaxis()->GetXmax());
        fit->SetParNames("Fit_Offset", "Fit_Slope");
        
        // Reasonable initial guesses for the fit parameters
        fit->SetParameters(hratio->GetMean(), 0.0001);
        fit->SetParLimits(0, 0.01, 100);
        //fit->SetParLimits(1, 0.00001, 0.1);

        // Sometimes we have to fit the data a few times before the fit converges
        bool converged = 0;
        int ntries = 0;
        
        // Make sure the fit converges.
        while(!converged) 
        {
          if(ntries >= 50) break;
          std::cout << Form("==== Fitting %s ==== \n", name.Data());
          hratio->Fit(fitname);

          TString sconverge = gMinuit->fCstatu.Data();
          converged = sconverge.Contains(TString("CONVERGED"));
          ntries++; 
        }


        // Seg fault if we do this stat box editing earlier
        // needs to be here
        gStyle->SetOptStat("n");
        TPaveStats* st = (TPaveStats*) hratio->GetListOfFunctions()->FindObject("stats");
        st->SetName("mystats");

        // legend == 0 default
        // directly overlays stats from pad1 if the user drew pad1 stats... 
        // should ifnd better alignment later
        st->SetX1NDC(0.52);
        st->SetY1NDC(0.83);
        st->SetX2NDC(0.70);
        st->SetY2NDC(0.93);

        // for fewz mass plots
        if(legend == 1)
        {
            st->SetX1NDC(0.425);
            st->SetY1NDC(0.280);
            st->SetX2NDC(0.882);
            st->SetY2NDC(0.458);
        }

        float slope_value = fit->GetParameter(1);
        float slope_err = fit->GetParError(1);
        float slope_percent_err = abs(slope_err/slope_value*100);

        // Remove the offset line
        TList* listOfLines = st->GetListOfLines();
        TText* title = st->GetLineWith("hratio");
        title->SetText(0.5, 0.5, "Ratio Info");
        //TText* tconst_s = st->GetLineWith("Slope");
        //listOfLines->Remove(tconst_o);
        //listOfLines->Remove(tconst_s);
        st->AddText(Form("Mean = %5.3f", scale));

        hratio->SetStats(0);
        TPaveStats* st2 = (TPaveStats*)st->Clone();
        st->Delete();
        gPad->Modified();
        pad1->cd();
        st2->Draw();
        gPad->Modified();
    }
    return c;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

ZCalibration::ZCalibration()
{
    massname = "mass";
    xname = "phi_plus";
    massmin = 86.2;
    massmax = 96.2;
    massbins = 50;
    xmin = -3.14;
    xmax = 3.14;
    xbins = 25;
    fitsig = 1;

    histos = std::vector<TH1D*>();
    binning = std::vector<Float_t>();
    vfis = std::vector<VoigtFitInfo>();
    init();
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

ZCalibration::ZCalibration(TString xname, TString massname, Float_t fitsig, Float_t massmin, 
                           Float_t massmax, Int_t massbins, Float_t xmin, Float_t xmax, Int_t xbins)
{
    this->xname = xname;
    this->massname = massname;
    this->massmin = massmin;
    this->massmax = massmax;
    this->massbins = massbins;
    this->xmin = xmin;
    this->xmax = xmax;
    this->xbins = xbins;
 
    this->fitsig = fitsig;

    histos = std::vector<TH1D*>();
    binning = std::vector<Float_t>();
    vfis = std::vector<VoigtFitInfo>();
    init();
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void ZCalibration::init()
{
    TString basename = "Z_Peak_Calibration_";
    basename+=massname+"_"+xname;

    binning.push_back(xmin);
    Float_t interval = (xmax - xmin)/xbins;
    Float_t upper_bin_edge = xmin;

    // Set up the bin edges 
    for(unsigned int i=0; i<xbins; i++)
    {
        upper_bin_edge += interval;
        binning.push_back(upper_bin_edge);
    }

    // print binning for debugging
    //std::cout << "printing binning..." << std::endl;
    //for(unsigned int i=0; i<binning.size(); i++)
    //{
    //    std::cout << binning[i] << std::endl;
    //}

    // Set up the TH1Ds
    for(unsigned int i=0; i<binning.size()-1; i++)
    {
        TString range = Form("_%5.2f_to_%5.2f", binning[i], binning[i+1]);
        TString title = basename+range;
        title.ReplaceAll(" ","");
        title.ReplaceAll(".","p");
        histos.push_back(new TH1D(title, title, massbins, massmin, massmax));
        histos[i]->GetXaxis()->SetTitle(massname);
    }
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

int ZCalibration::whichTH1D(Float_t xvalue)
{
    Float_t interval = (xmax - xmin)/xbins;
    int bin = (xvalue-xmin)/interval;

    if(xvalue == xmin) bin = 0;
    if(xvalue == xmax) bin = binning.size()-2;
    if(xvalue < xmax && bin > binning.size()-2) bin = binning.size()-2;

    if(xvalue < xmin || xvalue > xmax) bin = -9999;

    return bin;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void ZCalibration::fill(Float_t xvalue, Float_t massvalue)
{
    int h = whichTH1D(xvalue);
    if(h>=histos.size())
    {
        std::cout << Form("whichHist: %d, histos.size(): %d, binning.size(): %d, xmin: %f, xmax: %12.10f, xvalue: %12.10f, massvalue: %f \n", 
                           h, histos.size(), binning.size(), xmin, xmax, xvalue, massvalue); 
    }
    if(h>=0) 
        histos[h]->Fill(massvalue);
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void ZCalibration::fill(Float_t xvalue, Float_t massvalue, double weight)
{
    int h = whichTH1D(xvalue);
    if(h>=histos.size())
    {
        std::cout << Form("whichHist: %d, histos.size(): %d, binning.size(): %d, xmin: %f, xmax: %12.10f, xvalue: %12.10f, massvalue: %f \n", 
                           h, histos.size(), binning.size(), xmin, xmax, xvalue, massvalue); 
    }
    if(h>=0) 
        histos[h]->Fill(massvalue, weight);
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

VoigtFitInfo ZCalibration::fit(TH1D* inhist, Float_t x, Float_t x_err)
{
    VoigtFitInfo vfi;
    vfi.x = x;
    vfi.x_err = x_err;
    
    // Grab the results we are interested in.
    vfi.hmean = inhist->GetMean();
    vfi.hmean_err = inhist->GetMeanError();
    
    vfi.hrms = inhist->GetRMS();
    vfi.hrms_err = inhist->GetRMSError();
    
    // Display fit info on canvas.
    gStyle->SetOptFit(0011);
    // Fit things.
    TString fitname = TString("fit_")+massname+"_"+xname+Form("_%5.2f", x);
    fitname.ReplaceAll(" ", "");
    fitname.ReplaceAll(".", "p");
    TF1* fit = new TF1(fitname, "[0]*TMath::Voigt(x - [1], [2], [3])", 87, 95);
    fit->SetParNames("Constant", "Mean", "Sigma", "Gamma");
    
    // Reasonable initial guesses for the fit parameters
    fit->SetParameters(inhist->GetMaximum(), inhist->GetMean(), inhist->GetRMS(), vfi.gamma);
    fit->SetParLimits(2, 0.001, 2 * inhist->GetRMS());

    // Fix the intrinsic width to the theoretical value
    fit->FixParameter(3, 0.083985);

    // Sometimes we have to fit the data a few times before the fit converges
    bool converged = 0;
    int ntries = 0;
    
    // Make sure the fit converges.
    while(!converged) 
    {
      if(ntries >= 50) break;
      //std::cout << Form("==== Fitting %s ==== \n", inhist->GetTitle());
      fit->SetParameter(1, inhist->GetRMS());
      inhist->Fit(fitname, "q");

      // Fit to a width of fitsig sigmas
      inhist->Fit(fitname,"q","", fit->GetParameter(1) - fitsig*fit->GetParameter(2), fit->GetParameter(1) + fitsig*fit->GetParameter(2));
      TString sconverge = gMinuit->fCstatu.Data();
      converged = sconverge.Contains(TString("CONVERGED"));
      ntries++; 
    }

    vfi.fit = fit;
    vfi.vmean = fit->GetParameter(1);
    vfi.vmean_err = fit->GetParError(1);
    vfi.vsigma = fit->GetParameter(2);
    vfi.vsigma_err = fit->GetParError(2);
    vfi.vgamma = fit->GetParameter(3);
    vfi.vgamma_err = fit->GetParError(3);

    return vfi;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void ZCalibration::fit()
{
    for(unsigned int h=0; h<histos.size(); h++)
    {
        VoigtFitInfo vfi = fit(histos[h], (binning[h+1] + binning[h])/2, binning[h+1] - binning[h]);
        vfis.push_back(vfi);
    }
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void ZCalibration::plot()
{
    fit();
    mean_vs_x = new TGraphErrors();
    resolution_vs_x = new TGraphErrors();
    for(unsigned int i=0; i<vfis.size(); i++)
    {
        mean_vs_x->SetPoint(i, vfis[i].x, vfis[i].vmean);
        mean_vs_x->SetPointError(i, vfis[i].x_err, vfis[i].vmean_err);

        resolution_vs_x->SetPoint(i, vfis[i].x, vfis[i].vsigma);
        resolution_vs_x->SetPointError(i, vfis[i].x_err, vfis[i].vsigma_err);
    }
    mean_vs_x->SetTitle(massname+"_fit_mean_vs_"+xname);
    mean_vs_x->SetName("mean_"+massname+"_"+xname);
    mean_vs_x->GetXaxis()->SetTitle(xname);
    mean_vs_x->GetYaxis()->SetTitle(massname+"_mean");

    resolution_vs_x->SetTitle(massname+"_fit_resolution_vs_"+xname);
    resolution_vs_x->SetName("resolution_"+massname+"_"+xname);
    resolution_vs_x->GetXaxis()->SetTitle(xname);
    resolution_vs_x->GetYaxis()->SetTitle(massname+"_resolution");
}
