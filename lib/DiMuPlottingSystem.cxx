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

#include "Sample.h"
#include "DiMuPlottingSystem.h"

#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TH1F.h"
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

DiMuPlottingSystem::DiMuPlottingSystem(Sample* isample)
{
    sample = isample;
    list = 0;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

DiMuPlottingSystem::DiMuPlottingSystem(Sample* isample, TEntryList* ilist)
{
    sample = isample;
    list = ilist;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

DiMuPlottingSystem::~DiMuPlottingSystem() {
  // free pointed to memory!
  if (list != 0) {
    delete list;
  }
  if (sample !=0) {
    delete sample;
  }
}

///////////////////////////////////////////////////////////////////////////
// _______________________Other Functions________________________________//
///////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::initialize(Sample* isample, TEntryList* ilist)
{
    sample = isample;
    list = ilist;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TEntryList* DiMuPlottingSystem::getEntryList()
{
// Return the entry list currently associated with the sample.
    return list;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::setEntryList(TEntryList* ilist) {
// Associate some entrylist with the current sample. 
  list = ilist;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

int DiMuPlottingSystem::getEntry(int i)
{
  // Get the appropriate entry in the collection. If we applied cuts
  // we only want entries that passed the cuts, those in the entrylist.
  // This is useful if we want to iterate only over those events that passed
  // the cuts, which is usually the case.

    if(list != 0)
    {
        return sample->getEntry(i, list);
    } 
    else
    {
        return sample->getEntry(i);
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

int DiMuPlottingSystem::getNumEvents()
{
// Get the total number of events in the collection. If cuts have been applied
// the events we care about are in the TEntryList associated with the sample.

    if(list != 0)
    {
        return list->GetN();
    } 
    else
    {
        return sample->tree->GetEntries();
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TH1F* DiMuPlottingSystem::hist1D(TString name, TString xtitle, TString title, float& varToPlot, int nbins, int xmin, int xmax)
{
// Return a 1D histogram of varToPlot. Don't weight the events at all.

    TH1F* hist = new TH1F(name, title, nbins, xmin, xmax);
    //hist->Sumw2();
    hist->GetXaxis()->SetTitle(title);
    int sum_events = 0;

    for(int i=0; i<getNumEvents()/reductionFactor; i++)
    {
        getEntry(i);
        sum_events+=1;
        hist->Fill(varToPlot);
    }

    std::cout << std::endl;
    std::cout << "sum_events: " << sum_events << std::endl;
    std::cout << "hist " << hist->GetName() << " numEntries: " << hist->GetEntries() << std::endl;
    std::cout << "hist " << hist->GetName() << " Integral: " << hist->Integral() << std::endl;
    std::cout << std::endl;
    hist->SetDirectory(0);
    
    return hist;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TH1F* DiMuPlottingSystem::hist1D(TString name, TString xtitle, TString title, int& varToPlot, int nbins, int xmin, int xmax)
{
// Return a 1D histogram of varToPlot. Don't weight the events at all.

    TH1F* hist = new TH1F(name, title, nbins, xmin, xmax);
    //hist->Sumw2();
    hist->GetXaxis()->SetTitle(title);
    int sum_events = 0;

    for(int i=0; i<getNumEvents()/reductionFactor; i++)
    {
        getEntry(i);
        sum_events+=1;
        hist->Fill(varToPlot);
    }

    std::cout << std::endl;
    std::cout << "sum_events: " << sum_events << std::endl;
    std::cout << "hist " << hist->GetName() << " numEntries: " << hist->GetEntries() << std::endl;
    std::cout << "hist " << hist->GetName() << " Integral: " << hist->Integral() << std::endl;
    std::cout << std::endl;
    hist->SetDirectory(0);
    
    return hist;
}
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TH1F* DiMuPlottingSystem::weightedHist1D(TString name, TString xtitle, TString title, float& varToPlot, int nbins, int xmin, int xmax)
{
// Return a 1D histogram of varToPlot where each entry is weighted by the nnlo generated weights (+1/-1) and the PU weights via lumiWeights.

    // Data samples don't have to account for pu weighting or gen weighting
    if(sample->sampleType.Contains("data")) return hist1D(name, xtitle, title, varToPlot, nbins, xmin, xmax); 

    TH1F* hist = new TH1F(name, title, nbins, xmin, xmax);
    //hist->Sumw2();
    hist->GetXaxis()->SetTitle(title);
    int sum_weights = 0;
    int sum_events = 0;


    for(int i=0; i<getNumEvents()/reductionFactor; i++)
    {
        getEntry(i);

        sum_weights+=sample->lumiWeights->weight(sample->vars.nPU)*sample->vars.genWeight;
        sum_events+=1;

        if(i < 20) std::cout << i << ": genWeight = " << sample->vars.genWeight << std::endl;
        if(i < 20) std::cout << i << ": varToPlot = " << varToPlot << std::endl;
        if(i < 20) std::cout << i << ": reco1.pt = " << sample->vars.reco1.pt << std::endl;
        if(i < 20) std::cout << i << ": reco1.eta = " << sample->vars.reco1.eta << std::endl;
        if(i < 20) std::cout << i << ": reco2.pt = " << sample->vars.reco2.pt << std::endl;
        if(i < 20) std::cout << i << ": reco2.eta = " << sample->vars.reco2.eta << std::endl;
        if(i < 20) std::cout << i << ": vertexInfo.nVertices = " << sample->vars.vertices.nVertices << std::endl;
        if(i < 20) std::cout << i << ": nPU = " << sample->vars.nPU << std::endl;
        if(i < 20) std::cout << i << ": lumiWeight = " << sample->lumiWeights->weight(sample->vars.nPU) << std::endl;

        hist->Fill(varToPlot, 1.0*sample->vars.genWeight*sample->lumiWeights->weight(sample->vars.nPU));
    }

    std::cout << std::endl;
    std::cout << "sum_events: " << sum_events << std::endl;
    std::cout << "sum_weights: " << sum_weights << std::endl;
    std::cout << "hist " << hist->GetName() << " numEntries: " << hist->GetEntries() << std::endl;
    std::cout << "hist " << hist->GetName() << " Integral: " << hist->Integral() << std::endl;
    std::cout << "hist " << hist->GetName() << " N/Integral: " << hist->GetEntries()/hist->Integral() << std::endl;
    std::cout << std::endl;
    hist->SetDirectory(0);
    
    return hist;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TH1F* DiMuPlottingSystem::weightedHist1D(TString name, TString xtitle, TString title, int& varToPlot, int nbins, int xmin, int xmax)
{
// Return a 1D histogram of varToPlot where each entry is weighted by the nnlo generated weights (+1/-1) and the PU weights via lumiWeights.

    // Data samples don't have to account for pu weighting or gen weighting
    if(sample->sampleType.Contains("data")) return hist1D(name, xtitle, title, varToPlot, nbins, xmin, xmax); 

    TH1F* hist = new TH1F(name, title, nbins, xmin, xmax);
    //hist->Sumw2();
    hist->GetXaxis()->SetTitle(title);
    int sum_weights = 0;
    int sum_events = 0;

    for(int i=0; i<getNumEvents()/reductionFactor; i++)
    {
        getEntry(i);

        sum_weights+=sample->lumiWeights->weight(sample->vars.nPU)*sample->vars.genWeight;
        sum_events+=1;

        if(i < 20) std::cout << i << ": genWeight = " << sample->vars.genWeight << std::endl;
        if(i < 20) std::cout << i << ": varToPlot = " << varToPlot << std::endl;
        if(i < 20) std::cout << i << ": reco1.pt = " << sample->vars.reco1.pt << std::endl;
        if(i < 20) std::cout << i << ": reco1.eta = " << sample->vars.reco1.eta << std::endl;
        if(i < 20) std::cout << i << ": reco2.pt = " << sample->vars.reco2.pt << std::endl;
        if(i < 20) std::cout << i << ": reco2.eta = " << sample->vars.reco2.eta << std::endl;
        if(i < 20) std::cout << i << ": vertexInfo.nVertices = " << sample->vars.vertices.nVertices << std::endl;
        if(i < 20) std::cout << i << ": nPU = " << sample->vars.nPU << std::endl;
        if(i < 20) std::cout << i << ": lumiWeight = " << sample->lumiWeights->weight(sample->vars.nPU) << std::endl;

        hist->Fill(varToPlot, 1.0*sample->vars.genWeight*sample->lumiWeights->weight(sample->vars.nPU));
    }

    std::cout << std::endl;
    std::cout << "sum_events: " << sum_events << std::endl;
    std::cout << "sum_weights: " << sum_weights << std::endl;
    std::cout << "hist " << hist->GetName() << " numEntries: " << hist->GetEntries() << std::endl;
    std::cout << "hist " << hist->GetName() << " Integral: " << hist->Integral() << std::endl;
    std::cout << "hist " << hist->GetName() << " N/Integral: " << hist->GetEntries()/hist->Integral() << std::endl;
    std::cout << std::endl;
    hist->SetDirectory(0);
    
    return hist;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////


double DiMuPlottingSystem::getScaleFactor(float luminosity)
{
// Calculate the scale factor according to the integrated luminosity of the data and the cross section of the sample.
    if(sample->sampleType.Contains("data")) return 1;
    Double_t scalefactor = luminosity*sample->xsec/(sample->nOriginalWeighted/reductionFactor);
    return scalefactor;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::scaleHistByXsec(TH1F* hist, float luminosity)
{
// Scale the histogram according to the cross section and luminosity
    hist->Scale(getScaleFactor(luminosity));
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

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

THStack* DiMuPlottingSystem::overlay(TList* ilist, TString title, TString xaxistitle, TString yaxistitle)
{
// Assumes data is in the last ilist location so that it appears on top of the other histograms.


// Wide bottom left
//  TLegend* l = new TLegend(1, 0.15, 0.7, 1.25, "", "brNDC");
 
  // Square-ish top right
  //TLegend* l = new TLegend(0.68, 0.56, 0.88, 0.87, "", "brNDC");

  // Wide rectangle top right
  TLegend* l = new TLegend(0.47, 0.56, 0.88, 0.88, "", "brNDC");
  THStack* stack = new THStack();
  stack->SetTitle(title);

  // Square-ish top left
  //TLegend* l = new TLegend(0.13, 0.56, 0.33, 0.88, "", "brNDC");

  TIter next(ilist);
  TObject* object = 0;
  int i=0;

  std::vector<int> colors;
  colors.push_back(2);
  colors.push_back(4);
  colors.push_back(6);
  colors.push_back(7);
  colors.push_back(36);
  colors.push_back(50);
  colors.push_back(30);
  colors.push_back(9);
  colors.push_back(29);
  colors.push_back(3);

  while ((object = next()))
  {
      TH1F* hist = (TH1F*) object;
      hist->SetStats(0);
      hist->SetFillColor(colors[i]);
      hist->SetLineColor(colors[i]);

      // Print name + num events in legend
      TString legend_entry = TString(hist->GetName());
      //legend_entry.Form("%s %d", v[i]->GetName(), (int)v[i]->GetEntries());
      //std::cout << v[i]->GetName() << ": " << v[i]->GetEntries() << std::endl;
      //std::cout << legend_entry << std::endl;
      //if((i+1)==10) v[i]->SetFillColor(50);
      //v[i]->SetLineWidth(3);
 

      // Assuming data is in the last location
      if(object == ilist->Last())
      {
          hist->SetMarkerStyle(20);
          hist->SetLineColor(1);
          hist->SetFillColor(0);

          stack->Draw("hist");
          stack->GetXaxis()->SetTitle(xaxistitle);
          stack->GetYaxis()->SetTitle(yaxistitle);
          gPad->Modified();

          hist->SetStats(1);
	  hist->Draw("epsames");
          gPad->Update();
          TPaveStats *st=(TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
          st->SetX1NDC(0.69);
          st->SetY1NDC(0.37);
          st->SetX2NDC(0.88);
          st->SetY2NDC(0.53);
          gPad->Modified();
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

TH1F* DiMuPlottingSystem::addHists(TList* ilist, TString name)
{
// Add all of the histograms in the list into a single histogram

    TH1F* htotal = (TH1F*)ilist->First()->Clone(name);

    TIter next(ilist);
    TObject* object = 0;

    while ((object = next()))
    {
        TH1F* h = (TH1F*) object;
        if(object != ilist->First()) htotal->Add(h);
    }
    return htotal;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

TCanvas* DiMuPlottingSystem::stackedHistogramsAndRatio(TList* ilist, TString title, TString xaxistitle, TString yaxistitle)
{
    // Define two gaussian histograms. Note the X and Y title are defined
    // at booking time using the convention "Hist_title ; X_title ; Y_title"

    // Define the Canvas
    TCanvas *c = new TCanvas();
    c->SetCanvasSize(800,800);

    TH1F* first = (TH1F*) ilist->At(0);
    TH1F* last = (TH1F*) ilist->Last();

    // force the histograms to match
    // Right now this only scales the data to match the major player DYJetsToLL assumed to be in the first location
    Double_t scale = last->Integral()/first->Integral();
    std::cout << "########## Scale factor for " << first->GetName() << ": " << scale << std::endl;
    // Uncomment the next line to implement the hard scaling
    //first->Scale(scale);

     // Upper pad
    TPad *pad1 = new TPad(TString("pad1")+first->GetName(), "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0); 
    //pad1->SetGridx();         // Vertical grid
    //pad1->SetGridy();         // horizontal grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad

    // Draw stacked histograms here
    first->SetStats(0);          // No statistics on upper plot
    THStack* stack = overlay(ilist, title, xaxistitle, yaxistitle);

    // Not sure how to work with this
    // Do not draw the Y axis label on the upper plot and redraw a small
    // axis instead, in order to avoid the first label (0) to be clipped.
    //v[0]->GetYaxis()->SetLabelSize(0.);
    //TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
    //axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    //axis->SetLabelSize(15);
    //axis->Draw();

    // Lower pad
    c->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad(TString("pad2")+first->GetName(), "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy(); // horizontal grid
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad

    // Define the ratio plot
    // Clone the data histogram from the last location in the vector
    // The remainder of the vector consits of MC samples
    TH1F* hratio = (TH1F*)last->Clone("hratio");

    // We want a ratio to the sum of the MC in order to compare data to the net contribution from all MC, BG+SIG
    //TList* sublist = (TList*)list->Clone();
    //sublist->Remove(sublist->Last());
    //TH1F* hadd = addHists(sublist, "mctotalhist"); 
    
    // Create the ratio plot an easier way using THStack::GetStack()->Last()
    TH1F* hadd = (TH1F*)stack->GetStack()->Last();

    // scale the histograms to match
    scale = hratio->Integral()/hadd->Integral();
    std::cout << "########## Scale factor for MC stack: " << scale << std::endl;
    //hadd->Scale(scale);

    hratio->SetLineColor(kBlack);
    hratio->SetMinimum(0.58);  // Define Y ..
    hratio->SetMaximum(1.42); // .. range
    hratio->Sumw2();
    hratio->SetStats(0);      // No statistics on lower plot
    hratio->Divide(hadd);
    hratio->SetMarkerStyle(20);
    hratio->Draw("ep");       // Draw the ratio plot

    // Y axis v[v.size()-1] plot settings
    last->GetYaxis()->SetTitleSize(20);
    last->SetTitleFont(43);
    last->GetYaxis()->SetTitleOffset(1.55);

    // Ratio plot (hratio) settings
    hratio->SetTitle(""); // Remove the ratio title

    // Y axis ratio plot settings
    hratio->GetYaxis()->SetTitle("Data/MC");
    hratio->GetYaxis()->SetNdivisions(505);
    hratio->GetYaxis()->SetTitleSize(20);
    hratio->GetYaxis()->SetTitleFont(43);
    hratio->GetYaxis()->SetTitleOffset(1.55);
    hratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hratio->GetYaxis()->SetLabelSize(15);

    // X axis ratio plot settings
    hratio->GetXaxis()->SetTitleSize(20);
    hratio->GetXaxis()->SetTitleFont(43);
    hratio->GetXaxis()->SetTitleOffset(4.);
    hratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hratio->GetXaxis()->SetLabelSize(15);

    return c;
}
