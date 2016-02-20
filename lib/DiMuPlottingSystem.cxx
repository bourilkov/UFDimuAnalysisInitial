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

#include "DataFormats.h"
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

DiMuPlottingSystem::DiMuPlottingSystem(TString infilen, TString dps_name)
{
    infilename = infilen;
    dpsname = dps_name;
    treename = TString("dimuons/tree");
    infile = new TFile(infilename);
    tree = (TTree*)infile->Get(treename);
    if(tree == 0)
    {
        treename = TString("tree");
        tree = (TTree*)infile->Get(treename);
    }

    elist = 0;
    lumiWeights = 0;
    
    setBranchAddresses();
//    setNumEvents();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

DiMuPlottingSystem::DiMuPlottingSystem(TString infilen, TString dps_name, TString treen)
{
    dpsname = dps_name;
    infilename = infilen;
    treename = treen;
    infile = new TFile(infilename);
    tree = (TTree*)infile->Get(treename);
    elist = 0;
    lumiWeights = 0;

    setBranchAddresses();
//    setNumEvents();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

DiMuPlottingSystem::~DiMuPlottingSystem() {
  // free pointed to memory!
  if (tree != 0) {
    delete tree;
  }
  if (infile !=0) {
    delete infile;
  }
}

///////////////////////////////////////////////////////////////////////////
// _______________________Other Functions________________________________//
///////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::setBranchAddresses()
{
    tree->SetBranchAddress("recoCandMass", &recoCandMass);
    tree->SetBranchAddress("recoCandMassPrime", &recoCandMassPrime);
    tree->SetBranchAddress("nPU", &nPU);
    tree->SetBranchAddress("recoCandPt", &recoCandPt);
    tree->SetBranchAddress("vertexInfo", &vertexInfo);
    tree->SetBranchAddress("reco1", &reco1);
    tree->SetBranchAddress("reco2", &reco2);
    tree->SetBranchAddress("genM1ZpreFSR", &genM1ZpreFSR);
    tree->SetBranchAddress("genM2ZpreFSR", &genM2ZpreFSR);
    tree->SetBranchAddress("genM1ZpostFSR", &genM1ZpostFSR);
    tree->SetBranchAddress("genM2ZpostFSR", &genM2ZpostFSR);
    tree->SetBranchAddress("genZpreFSR", &genZpreFSR);
    tree->SetBranchAddress("genZpostFSR", &genZpostFSR);
    tree->SetBranchAddress("genWeight", &genWeight);
    tree->SetBranchAddress("pfJets", &rawJets);
    tree->SetBranchAddress("met", &met);
    tree->SetBranchAddress("recoCandMassPF", &recoCandMassPF);
    tree->SetBranchAddress("recoCandPtPF", &recoCandPtPF);
}
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::setNumEvents()
{
// Get the number of original events from the meta data tree
    TTree* metadata = (TTree*)infile->Get("dimuons/metadata");
    metadata->Draw("sumEventWeights>>eweights_"+dpsname);
    TH1F* sumEventWeightsHist = (TH1F*) gDirectory->Get("eweights_"+dpsname); 

    // There are many ttrees combined so the histogram has a numEvents entry for each
    // file. The total number is equal to the average number times the total number of entries.
    NoriginalWeighted = sumEventWeightsHist->GetEntries()*sumEventWeightsHist->GetMean();

    metadata->Draw("originalNumEvents>>nevents_"+dpsname);
    TH1F* nEventsHist = (TH1F*) gDirectory->Get("nevents_"+dpsname); 
    Noriginal = nEventsHist->GetEntries()*nEventsHist->GetMean();
    if(metadata !=0) delete metadata;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TEntryList DiMuPlottingSystem::getEntryList(TCut cuts)
{
// Apply cuts to form an entrylist then return this entrylist.
    tree->Draw(">>l" + dpsname, cuts, "entrylist");
    TEntryList* list = (TEntryList*)gDirectory->Get("l" + dpsname);
    return *list;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::applyCuts(TCut cuts) {
  // Apply whatever cuts are input and make sure the tree only uses these events.
  tree -> Draw(">>l"+dpsname, cuts, "entrylist");
  elist = (TEntryList*) gDirectory -> Get("l"+dpsname);
  tree -> SetEntryList(elist);
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::applyRun1Cuts() {
  // Apply the cuts from run1
  tree -> Draw(">>l"+dpsname, run1cuts, "entrylist");
  elist = (TEntryList*) gDirectory -> Get("l"+dpsname);
  tree -> SetEntryList(elist);
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::resetTree()
{
  // Reset the tree as if no cuts were ever applied.

  tree->SetEntryList(0);
  delete elist;
  elist = 0;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::setEntryList(TEntryList* list)
{
  // Reset the tree as if no cuts were ever applied.

  tree->SetEntryList(list);
  elist = list;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Int_t DiMuPlottingSystem::getEntry(Int_t i)
{
  // Get the appropriate entry in the collection. If we applied cuts
  // we only want entries that passed the cuts, those in the entrylist.
  // This is useful if we want to iterate only over those events that passed
  // the cuts, which is usually the case.

  if( elist != 0)
    {
      tree -> GetEntry(elist->GetEntry(i));
      return elist -> GetEntry(i);
    } 
  else
    {
      tree -> GetEntry(i);
      return i;
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Int_t DiMuPlottingSystem::getNumEvents()
{
  // We want to determine the number of events in the collection we care about.
  // If we apply cuts we only care about the set that passes the cuts.

  if( elist != 0)
    {
      return elist->GetN();
    } 
  else
    {
      return tree->GetEntries();
    }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TGraphAsymmErrors* DiMuPlottingSystem::drawEfficiency(TH1F* pass, TH1F* total, TString title, TString xtitle)
{
  // Given numerator and denominator histograms, draw the efficiency.

  Float_t npass  = pass->GetEntries();
  Float_t ntotal = total->GetEntries();
    
  std::cout << std::endl;
  std::cout << title << std::endl;
  std::cout << "npass  : " << npass << std::endl;
  std::cout << "ntotal : " << ntotal << std::endl; 

  // Plot the efficiency results from pass/total
  TCanvas* c = new TCanvas();
  c->SetName(title);
  TGraphAsymmErrors* eff = new TGraphAsymmErrors(pass, total, "cp");
  c->SetGridx(kTRUE);
  c->SetGridy(kTRUE);
  eff->SetTitle(title);
  eff->SetName(title);
  eff->GetXaxis()->SetTitle(xtitle);
  eff->GetYaxis()->SetTitle("Efficiency");
  eff->Draw("AL");
  return eff;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool DiMuPlottingSystem::isGenMatched(_MuonInfo reco, Float_t dRcutoff)
{
  // Determine whether the reco muon reconstructs a generated muon by determining
  // whether dR is less than the dRcutoff value.

  // We check to see if there is an acceptable dR value with a same signed preFSR muon
  // or a same signed postFSR muon.

  // Used by efficiency plots

  _TrackInfo preGen1 = genM1ZpreFSR;
  _TrackInfo preGen2 = genM2ZpreFSR;
  _TrackInfo postGen1 = genM1ZpostFSR;
  _TrackInfo postGen2 = genM2ZpostFSR;


  Float_t dEta[4] = {reco.eta - preGen1.eta, reco.eta - preGen2.eta, reco.eta - postGen1.eta, reco.eta - postGen2.eta};
  Float_t dPhi[4] = {reco.phi - preGen1.phi, reco.phi - preGen2.phi, reco.phi - postGen1.phi, reco.phi - postGen2.phi};
  // Determine whether the charges between the reco and the candidate gen muon match +1 for match -1 for mismatch
  Int_t charge[4] =  {reco.charge*preGen1.charge, reco.charge*preGen2.charge, reco.charge*postGen1.charge, reco.charge*postGen2.charge};

  // Calculate the dR values between the reco and all the candidate gen muons
  std::vector<Float_t> deltaRvec = std::vector<Float_t>();
  for(unsigned int i=0; i<4; i++)
    {
      Float_t deltaRsq = dEta[i]*dEta[i] + dPhi[i]*dPhi[i];
      Float_t deltaR = TMath::Sqrt(deltaRsq);
      deltaRvec.push_back(deltaR);
    }

  for(int i=3; i>=0; i--)
    {
      // The reco candidate must have a small dR value with a same signed gen muon
      // If this criteria is met then the reco muon reconstructs a gen muon.
      if(deltaRvec[i] < dRcutoff && charge[i] ==1)
        {
	  if(genM1ZpreFSR.pt == -999 && genM1ZpostFSR.pt == -999) std::cout << "$$$$$$$$$ -999 and dR passes" << std::endl;
	  return true;
        }
    }

  return false;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

bool DiMuPlottingSystem::fillRecoGivenECI(_MuonInfo reco, EffCutInfo eci)
{
  // Used when plotting the efficiencies. Tells us whether to put the event into the histo.
 
  // If we require Gen Matching and the reco does not match a gen then return false.
  if(eci.match && !isGenMatched(reco, eci.dRcutoff)) return false;

  // If we require the muon to pass the trigger and the muon does not pass the trigger return false.
  if(eci.trigger && !reco.isHltMatched[eci.triggernum]) return false;

  // Otherwise return true.
  return true;

}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::plotEfficiency(EffCutInfo numerator, EffCutInfo denominator, TString title, TString varToPlotAgainst, Int_t bins, Float_t min, Float_t max)
{
  // Make some efficiency plots. Had to make a special function due to trigger matching and dR calculations. 

  // This gives us the list of entries for the events that pass the numerator cuts.
  TEntryList nlist = getEntryList(numerator.cuts);

  // This gives us the list of entries for the events that pass the denominator cuts.
  TEntryList dlist = getEntryList(denominator.cuts);

  Int_t n = nlist.Next();
  Int_t d = dlist.Next();

  // Print out some debugging info.
  std::cout << "&n: " << &n << std::endl;
  std::cout << "&d: " << &d << std::endl;
  std::cout << "n size: " << nlist.GetN() << std::endl;
  std::cout << "d size: " << dlist.GetN() << std::endl;

  TString name    = TString(title);

  // Initialize with default values.
  Float_t fill1 = -999;
  Float_t fill2 = -999;

  std::cout << varToPlotAgainst << std::endl;
  std::cout << std::endl;

  // Get the number of events in the collection. If cuts have been applied it only
  // returns the number of events that pass the cuts.
  Int_t nentries = getNumEvents();
  std::cout << nentries << std::endl;

  // The numerator histogram
  TH1F* pass  = new TH1F("pass", "pass", bins, min, max);

    // The denominator histogram
    TH1F* total = new TH1F("total", "total", bins, min, max);

    for(Int_t i=0; i<nentries; i++) 
    {    
        // Get the entry using our DiMuPlottingSystem. If cuts have been applied
        // it gets the ith entry that passes the cuts. This might be buggy when cuts are applied and should be tested.
        // Should return the number of the tree entry and check as we iterate.
        // This will get the ith entry in the list but if cuts have been applied then the ith entry that passes the cuts
        // is not the ith entry in the tree which the n,d values refer to.
        // Hopefully this is fixed by returning tree_entry_id.
        Int_t tree_entry_id = getEntry(i);

        // If the number of the next numerator event is not the current event
        // and the next denominator event is not the current event then go to the next event.
        if(tree_entry_id!=n && tree_entry_id!=d) continue;

        // Determine which variable to bin the histogram in based upon user input
        if(varToPlotAgainst.EqualTo("Reco_Pt"))
        {
	  fill1 = reco1.pt;
	  fill2 = reco2.pt;
        }
      else if(varToPlotAgainst.EqualTo("Reco_Eta")) 
        {
	  fill1 = reco1.eta;
	  fill2 = reco2.eta;
        }
      else if(varToPlotAgainst.EqualTo("N_Primary_Vertices"))
        {
	  fill1 = vertexInfo.nVertices;
	  fill2 = vertexInfo.nVertices;
        }
      else
        {
	  std::cout << varToPlotAgainst << " is not recognized." << std::endl;
	  return;
        }
        // This event is an event that passes the denominator cuts
        if(tree_entry_id==d)
        {
	  // Check to see if the muons pass the trigger and gen muon matching criteria for the numerator.
	  if(fillRecoGivenECI(reco1, denominator)) total->Fill(fill1);
	  if(fillRecoGivenECI(reco2, denominator)) total->Fill(fill2);
	  d = dlist.Next();
        }
        // This event is an event that passes the numerator cuts
        if(tree_entry_id==n)
        {
	  // Check to see if the muons pass the trigger and gen muon matching criteria for the denominator.
	  if(fillRecoGivenECI(reco1, numerator)) pass->Fill(fill1);
	  if(fillRecoGivenECI(reco2, numerator))  pass->Fill(fill2);
	  n = nlist.Next();
        }

    }

  drawEfficiency(pass, total, title, varToPlotAgainst);
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

GaussFitInfo DiMuPlottingSystem::fitGaussHisto(TH1F* inhist, Float_t fitsig, Float_t x, Float_t x_err)
{
  // Given a 1D histogram, fit it with a Gaussian. Fit the gaussian to fitsig deviations.
  // Fill the fit info into the GaussFitInfo data structure and return it.

  GaussFitInfo gfi;
  gfi.x = x;
  gfi.x_err = x_err;

  // Grab the results we are interested in. 
  gfi.hmean = inhist->GetMean();
  gfi.hmean_err = inhist->GetMeanError();

  gfi.hrms = inhist->GetRMS();
  gfi.hrms_err = inhist->GetRMSError();

  // Display fit info on canvas.
  gStyle->SetOptFit(0011);
  // Fit things.
  TF1* fit = new TF1("fit", "[0]*TMath::Exp(-0.5*(((x-[1])/[2])^2))", -20, 20);

  fit->SetParName(0, "Constant");
  fit->SetParName(1, "Mean");
  fit->SetParName(2, "Sigma");

  fit->SetParameter(0, inhist->GetMaximum());
  fit->SetParameter(1, inhist->GetMean());
  fit->SetParameter(2, inhist->GetRMS());

  fit->SetParLimits(0, 1, 1*inhist->GetMaximum());
  //fit->SetParLimits(1, -2*inhist->GetMean(), 2*inhist->GetMean());
  fit->SetParLimits(2, 0.001, 1*inhist->GetRMS());

  bool converged = 0;
  int ntries = 0;

  // Make sure the fit converges.
  while(!converged)
    {
      if(ntries >= 50) break;
      std::cout << "==== " << inhist->GetTitle() << " ====" << std::endl;
      fit->SetParameter(2, inhist->GetRMS());
      inhist->Fit("fit");
      // Fit within a certain sigma range
      inhist->Fit("fit","","", fit->GetParameter(1)-fitsig*fit->GetParameter(2), fit->GetParameter(1)+fitsig*fit->GetParameter(2));
      TString sconverge = gMinuit->fCstatu.Data();
      converged = sconverge.Contains(TString("CONVERGED"));
      ntries++;
    }

  gfi.gmean = fit->GetParameter(1);
  gfi.gmean_err = fit->GetParError(1);
  gfi.gsigma = fit->GetParameter(2);
  gfi.gsigma_err = fit->GetParError(2);

  return gfi;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

VoigtFitInfo DiMuPlottingSystem::fitVoigtHisto(TH1F* inhist, Float_t fitsig, Float_t x, Float_t x_err) {
  // Given a 1D histogram, fit it with a Voigt profile. Fit the gaussian to fitsig deviations.
  // Fill the fit info into the GaussFitInfo data structure and return it.

  VoigtFitInfo vfi;
  vfi.x = x;
  vfi.x_err = x_err;

  // Grab the results we are interested in.
  vfi.hmean = inhist -> GetMean();
  vfi.hmean_err = inhist -> GetMeanError();

  vfi.hrms = inhist -> GetRMS();
  vfi.hrms_err = inhist -> GetRMSError();

  // Display fit info on canvas.
  gStyle->SetOptFit(0011);
  // Fit things.
  TF1* fit = new TF1("fit", "[0]*TMath::Voigt(x - [1], [2], [3])", 87, 95);
  fit -> SetParNames("Constant", "Mean", "Sigma", "Gamma");

  fit -> SetParameters(inhist -> GetMaximum(), inhist -> GetMean(), inhist -> GetRMS(), vfi.gamma);
  //fit -> SetParLimits(0, 1, 1000);
  //fit -> SetParLimits(1, 1, 1);
  fit -> SetParLimits(2, 0.001, 2 * inhist -> GetRMS());
  //fit -> SetParLimits(3, 0.001, 7);
  fit -> FixParameter(3, 0.083985);
  bool converged = 0;
  int ntries = 0;

  // Make sure the fit converges.
  while(!converged) {
    if(ntries >= 50) break;
    std::cout << "==== " << inhist->GetTitle() << " ====" << std::endl;
    fit -> SetParameter(1, inhist->GetRMS());
    inhist -> Fit("fit");
    inhist -> Fit("fit","","", fit -> GetParameter(1) - fitsig * fit -> GetParameter(2), fit -> GetParameter(1) + fitsig * fit -> GetParameter(2));
    TString sconverge = gMinuit -> fCstatu.Data();
    converged = sconverge.Contains(TString("CONVERGED"));
    ntries++;
  }
  
  vfi.vmean = fit -> GetParameter(1);
  vfi.vmean_err = fit -> GetParError(1);
  vfi.vsigma = fit -> GetParameter(2);
  vfi.vsigma_err = fit -> GetParError(2);
  vfi.vgamma = fit -> GetParameter(3);
  vfi.vgamma_err = fit -> GetParError(3);
  
  return vfi;
}
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TH1F* DiMuPlottingSystem::hist1D(TString name, TString varToPlot, TString bininfo, TCut cuts)
{
// Plot a 2D histogram and return it.
    tree->Draw(varToPlot+">>"+name+bininfo, cuts, "");
    TH1F* hist = (TH1F*)gDirectory->Get(name);

    // We have to clone the histogram since it goes out of scope otherwise.
    //hist = (TH1F*)hist->Clone();
    hist->SetDirectory(0);
    return hist;
}

TH1F* DiMuPlottingSystem::hist1D(TString name, TString varToPlot, TString bininfo, TString cuts)
{
// Plot a 2D histogram and return it.
    tree->Draw(varToPlot+">>"+name+bininfo, cuts, "");
    TH1F* hist = (TH1F*)gDirectory->Get(name);

    // We have to clone the histogram since it goes out of scope otherwise.
    //hist = (TH1F*)hist->Clone();
    hist->SetDirectory(0);
    return hist;
}

TH1F* DiMuPlottingSystem::hist1D(TString name, TString varToPlot, TString bininfo)
{
// Plot a 2D histogram and return it.
    tree->Draw(varToPlot+">>"+name+bininfo, "", "");
    TH1F* hist = (TH1F*)gDirectory->Get(name);

    // We have to clone the histogram since it goes out of scope otherwise.
    //hist = (TH1F*)hist->Clone();
    hist->SetDirectory(0);
    return hist;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TH1F* DiMuPlottingSystem::hist1D(TString name, TString xtitle, TString title, Float_t& varToPlot, Int_t nbins, Int_t xmin, Int_t xmax)
{
// Return a 1D histogram of varToPlot where each entry is weighted so that the MC pilup distribution matches the pileup distribution in data.

    TH1F* hist = new TH1F(name, title, nbins, xmin, xmax);
    //hist->Sumw2();
    hist->GetXaxis()->SetTitle(title);
    Int_t sum_events = 0;

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

TH1F* DiMuPlottingSystem::PUWeightedHist1D(TString name, TString xtitle, TString title, Float_t& varToPlot, Int_t nbins, Int_t xmin, Int_t xmax)
{
// Return a 1D histogram of varToPlot where each entry is weighted so that the MC pilup distribution matches the pileup distribution in data.

    TH1F* hist = new TH1F(name, title, nbins, xmin, xmax);
    //hist->Sumw2();
    hist->GetXaxis()->SetTitle(title);
    Int_t sum_weights = 0;
    Int_t sum_events = 0;

    for(int i=0; i<getNumEvents()/reductionFactor; i++)
    {
        getEntry(i);

        sum_weights+=lumiWeights->weight(nPU);
        sum_events+=1;

        if(i < 20) std::cout << i << ": genWeight = " << genWeight << std::endl;
        if(i < 20) std::cout << i << ": varToPlot = " << varToPlot << std::endl;
        if(i < 20) std::cout << i << ": reco1.pt = " << reco1.pt << std::endl;
        if(i < 20) std::cout << i << ": reco1.eta = " << reco1.eta << std::endl;
        if(i < 20) std::cout << i << ": reco2.pt = " << reco2.pt << std::endl;
        if(i < 20) std::cout << i << ": reco2.eta = " << reco2.eta << std::endl;
        if(i < 20) std::cout << i << ": vertexInfo.nVertices = " << vertexInfo.nVertices << std::endl;
        if(i < 20) std::cout << i << ": nPU = " << nPU << std::endl;
        if(i < 20) std::cout << i << ": lumiWeight = " << lumiWeights->weight(nPU) << std::endl;

        hist->Fill(varToPlot, 1.0*genWeight*lumiWeights->weight(nPU));
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

TH1F* DiMuPlottingSystem::PUWeightedHist1D(TString name, TString xtitle, TString title, Int_t& varToPlot, Int_t nbins, Int_t xmin, Int_t xmax)
{
// Return a 1D histogram of varToPlot where each entry is weighted so that the MC pilup distribution matches the pileup distribution in data.

    TH1F* hist = new TH1F(name, title, nbins, xmin, xmax);
    //hist->Sumw2();
    hist->GetXaxis()->SetTitle(title);

    Int_t sum_weights = 0;
    Int_t sum_events = 0;

    for(int i=0; i<getNumEvents()/reductionFactor; i++)
    {
        getEntry(i);

        sum_weights+=lumiWeights->weight(nPU);
        sum_events+=1;

        if(i < 20) std::cout << i << ": genWeight = " << genWeight << std::endl;
        if(i < 20) std::cout << i << ": varToPlot = " << varToPlot << std::endl;
        if(i < 20) std::cout << i << ": reco1.pt = " << reco1.pt << std::endl;
        if(i < 20) std::cout << i << ": reco1.eta = " << reco1.eta << std::endl;
        if(i < 20) std::cout << i << ": reco2.pt = " << reco2.pt << std::endl;
        if(i < 20) std::cout << i << ": reco2.eta = " << reco2.eta << std::endl;
        if(i < 20) std::cout << i << ": vertexInfo.nVertices = " << vertexInfo.nVertices << std::endl;
        if(i < 20) std::cout << i << ": nPU = " << nPU << std::endl;
        if(i < 20) std::cout << i << ": lumiWeight = " << lumiWeights->weight(nPU) << std::endl;

        hist->Fill(varToPlot, 1.0*genWeight*lumiWeights->weight(nPU));
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

TH1F* DiMuPlottingSystem::weightedHist1D(TString name, TString xtitle, TString title, Float_t& varToPlot, Int_t nbins, Int_t xmin, Int_t xmax)
{
// Return a 1D histogram of varToPlot where each entry is weighted so that the MC pilup distribution matches the pileup distribution in data.

    TH1F* hist = new TH1F(name, title, nbins, xmin, xmax);
    //hist->Sumw2();
    hist->GetXaxis()->SetTitle(title);

    for(int i=0; i<getNumEvents()/reductionFactor; i++)
    {
        getEntry(i);
        hist->Fill(varToPlot, 1.0*genWeight);
    }
    hist->SetDirectory(0);
    return hist;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

TH1F* DiMuPlottingSystem::weightedHist1D(TString name, TString xtitle, TString title, Int_t& varToPlot, Int_t nbins, Int_t xmin, Int_t xmax)
{
// Return a 1D histogram of varToPlot where each entry is weighted so that the MC pilup distribution matches the pileup distribution in data.

    TH1F* hist = new TH1F(name, title, nbins, xmin, xmax);
    //hist->Sumw2();
    hist->GetXaxis()->SetTitle(title);

    for(int i=0; i<getNumEvents()/reductionFactor; i++)
    {
        getEntry(i);
        hist->Fill(varToPlot, 1.0*genWeight);
    }
    hist->SetDirectory(0);
    return hist;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////


Double_t DiMuPlottingSystem::scaleByXsec(Float_t luminosity)
{
// Calculate the scale factor according to the integrated luminosity of the data and the cross section of the sample.
    Double_t scalefactor = luminosity*xsec/(NoriginalWeighted/reductionFactor);
    return scalefactor;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::scaleHistByXsec(TH1F* hist, Float_t luminosity)
{
// Scale the histogram according to the cross section and luminosity
    hist->Scale(scaleByXsec(luminosity));
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////
TH2F* DiMuPlottingSystem::hist2D(TString name, TString varsToPlot, TString bininfo, TCut cuts) {
  // Plot a 2D histogram and return it.
  tree -> Draw(varsToPlot + ">>" + name + bininfo, cuts, "colz");
  TH2F* hist = (TH2F*) gDirectory -> Get(name);

  // We have to clone the histogram since it goes out of scope otherwise.
  //hist = (TH2F* )hist -> Clone();
  hist->SetDirectory(0);
  return hist;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

std::vector<GaussFitInfo> DiMuPlottingSystem::plotVarVsGaussFit(TH2F* hist, TList* hlist, Float_t fitsig) {
  // Given a 2D histogram look at the distribution of the y values within each xbin and fit these with gaussians.
  // Save the gauss fit 1D histograms into a TList and return a vector with the fit info for each of the histograms.

  std::vector<GaussFitInfo> gfis = std::vector<GaussFitInfo>();

  std::cout << "plotVarVsGaussFit nentries: " << hist->GetEntries() << std::endl; 
  // i=0 to i<=hist->GetNbinsX()+1 to include overflow bins
  for(Int_t i=1; i<=hist->GetNbinsX(); i++)
    {
      Float_t x = (hist->GetXaxis()->GetBinLowEdge(i) + hist->GetXaxis()->GetBinLowEdge(i) + hist->GetXaxis()->GetBinWidth(i))  / 2;
      Float_t x_err = hist->GetXaxis()->GetBinWidth(i)/2;
      
      TH1F* h = (TH1F*)hist->ProjectionY(TString("py")+Form("%d", i), i, i);
      h->SetTitle( TString(hist->GetXaxis()->GetTitle()) + Form("_%5.2f_to_%5.2f", hist->GetXaxis()->GetBinLowEdge(i), hist->GetXaxis()->GetBinLowEdge(i) + hist->GetXaxis()->GetBinWidth(i)) );
      GaussFitInfo gfi = fitGaussHisto(h, fitsig, x, x_err); 
      gfis.push_back(gfi);
      hlist->Add(h);
    }
  return gfis;
}

std::vector<VoigtFitInfo> DiMuPlottingSystem::plotVarVsVoigtFit(TH2F* hist, TList* hlist, Float_t fitsig) {
  std::vector<VoigtFitInfo> vfis = std::vector<VoigtFitInfo>();

  std::cout << "plotVarVsVoigtFit nentries: " << hist -> GetEntries() << std::endl;

  for (Int_t i = 1; i <= hist -> GetNbinsX(); ++i) {
    Float_t x = (hist -> GetXaxis() -> GetBinLowEdge(i) + hist -> GetXaxis() -> GetBinLowEdge(i) + hist -> GetXaxis() -> GetBinWidth(i)) / 2;
    Float_t x_err = hist -> GetXaxis() -> GetBinWidth(i) /2;
    
    TH1F* h = (TH1F*) hist -> ProjectionY(TString("py") + Form("%d", i), i, i);
    h -> SetTitle(TString(hist -> GetXaxis() -> GetTitle()) + Form("_%5.2f_to_%5.2f", hist -> GetXaxis() -> GetBinLowEdge(i), hist -> GetXaxis() -> GetBinLowEdge(i) + hist -> GetXaxis() -> GetBinWidth(i)));
    std::cout << "=========== Entries for " << h -> GetTitle() << ": " << h -> GetEntries() << std::endl << std::endl;
    VoigtFitInfo vfi = fitVoigtHisto(h, fitsig, x, x_err);
    vfis.push_back(vfi);
    hlist -> Add(h);
  }
  return vfis;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::plotGFI(std::vector<GaussFitInfo>& gfis, std::vector<TGraphErrors*>& plots, TString name, TString xtitle, TString ytitle)
{
  // Given a vector of gauss fit information make some plots of the fit mean and fit sigma.
  // Return the mean and sigma plots in the passed TGraphErrors vector where [0] is the mean and [1] is sigma. 

  TGraphErrors* gmean = new TGraphErrors();
  TGraphErrors* gsig = new TGraphErrors();

  for(unsigned int i=0; i<gfis.size(); i++)
    {
      gmean->SetPoint(i, gfis[i].x, gfis[i].gmean);
      gmean->SetPointError(i, gfis[i].x_err, gfis[i].gmean_err);

      gsig->SetPoint(i, gfis[i].x, gfis[i].gsigma);
      gsig->SetPointError(i, gfis[i].x_err, gfis[i].gsigma_err);
    }
  plots.push_back(gmean);
  gmean->GetXaxis()->SetTitle(xtitle);
  gmean->GetYaxis()->SetTitle(ytitle + " Gaussian Fit Mean");
  gmean->SetTitle(ytitle + " Gaussian Fit Mean");
  gmean->SetName(name + "_gmean");

  plots.push_back(gsig);
  gsig->GetXaxis()->SetTitle(xtitle);
  gsig->GetYaxis()->SetTitle(ytitle + " Gaussian Fit Sigma");
  gsig->SetTitle(ytitle + " Gaussian Fit Sigma");
  gsig->SetName(name +"_gsigma");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::plotVFI(std::vector<VoigtFitInfo>& vfis, std::vector<TGraphErrors*>& plots, TString name, TString xtitle, TString ytitle)
{
  TGraphErrors* vmean = new TGraphErrors();
  TGraphErrors* vsig = new TGraphErrors();

  for(unsigned int i=0; i < vfis.size(); i++) {
    vmean -> SetPoint(i, vfis[i].x, vfis[i].vmean);
    vmean -> SetPointError(i, vfis[i].x_err, vfis[i].vmean_err);
    
    vsig -> SetPoint(i, vfis[i].x, vfis[i].vsigma);
    vsig -> SetPointError(i, vfis[i].x_err, vfis[i].vsigma_err);
  }
  
  plots.push_back(vmean);
  vmean -> GetXaxis() -> SetTitle(xtitle);
  vmean -> GetYaxis() -> SetTitle(ytitle + " Voigt Fit Mean");
  vmean -> SetTitle(ytitle + " Voigt Fit Mean");
  vmean -> SetName(name + "_vmean");

  plots.push_back(vsig);
  vsig -> GetXaxis() -> SetTitle(xtitle);
  vsig -> GetYaxis() -> SetTitle(ytitle + " Voigt Fit Sigma");
  vsig -> SetTitle(ytitle + " Voigt Fit Sigma");
  vsig -> SetName(name +"_vsigma");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

double DiMuPlottingSystem::newError(double data, double data_err, double dy, double dy_err) {
  double data_new_err = data_err / dy;
  double dy_new_err = data * dy_err / pow(dy, 2.0);
  double quadrature = pow(data_new_err, 2.0) + pow(dy_new_err, 2.0);
  double error =  sqrt(quadrature);
  return error;
};

void DiMuPlottingSystem::plotVFIRatio(std::vector<VoigtFitInfo>& vfidata, std::vector<VoigtFitInfo>& vfidy, std::vector<TGraphErrors*>& plots, TString name, TString xtitle, TString ytitle) {
  TGraphErrors* vmean = new TGraphErrors();
  TGraphErrors* vsig = new TGraphErrors();

  for(unsigned int i=0; i < vfidata.size(); i++) {
    vmean -> SetPoint(i, vfidata[i].x, vfidata[i].vmean / vfidy[i].vmean);
    vmean -> SetPointError(i, vfidata[i].x_err, DiMuPlottingSystem::newError(vfidata[i].vmean, vfidata[i].vmean_err, vfidy[i].vmean, vfidy[i].vmean_err));

    vsig -> SetPoint(i, vfidata[i].x, vfidata[i].vsigma / vfidy[i].vsigma);
    vsig -> SetPointError(i, vfidata[i].x_err, DiMuPlottingSystem::newError(vfidata[i].vsigma, vfidata[i].vsigma_err, vfidy[i].vsigma, vfidy[i].vsigma_err));
  }

  plots.push_back(vmean);
  vmean -> GetXaxis() -> SetTitle(xtitle);
  vmean -> GetYaxis() -> SetTitle(ytitle + " Voigt Fit Mean");
  vmean -> SetTitle(ytitle + " Voigt Fit Mean");
  vmean -> SetName(name + "_vmean_ratio");

  plots.push_back(vsig);
  vsig -> GetXaxis() -> SetTitle(xtitle);
  vsig -> GetYaxis() -> SetTitle(ytitle + " Voigt Fit Sigma");
  vsig -> SetTitle(ytitle + " Voigt Fit Sigma");
  vsig -> SetName(name +"_vsigma_ratio");
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void DiMuPlottingSystem::arrangeStatBox(TCanvas* c)
{
  // Align the stat box better.

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

THStack* DiMuPlottingSystem::overlay(TList* list, TString title, TString xaxistitle, TString yaxistitle)
{
// Assumes data is in the last list location so that it appears on top of the other histograms.


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

  TIter next(list);
  TObject* object = 0;
  int i=0;

  std::vector<Int_t> colors;
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
      if(object == list->Last())
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

TH1F* DiMuPlottingSystem::addHists(TList* list, TString name)
{
// Add all of the histograms in the list into a single histogram

    TH1F* htotal = (TH1F*)list->First()->Clone(name);

    TIter next(list);
    TObject* object = 0;

    while ((object = next()))
    {
        TH1F* h = (TH1F*) object;
        if(object != list->First()) htotal->Add(h);
    }
    return htotal;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

TCanvas* DiMuPlottingSystem::stackedHistogramsAndRatio(TList* list, TString title, TString xaxistitle, TString yaxistitle)
{
    // Define two gaussian histograms. Note the X and Y title are defined
    // at booking time using the convention "Hist_title ; X_title ; Y_title"

    // Define the Canvas
    TCanvas *c = new TCanvas();
    c->SetCanvasSize(800,800);

    TH1F* first = (TH1F*) list->At(0);
    TH1F* last = (TH1F*) list->Last();

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
    THStack* stack = overlay(list, title, xaxistitle, yaxistitle);

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
