//////////////////////////////////////////////////////////////////////////
//                             Sample.cxx                                //
//=======================================================================//
//                                                                       //
//        Keeps track of the sample information. TTree, xsec, etc.       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "Sample.h"

///////////////////////////////////////////////////////////////////////////
// _______________________Constructor/Destructor_________________________//
///////////////////////////////////////////////////////////////////////////

Sample::Sample() 
{
    lumiWeights = 0;
    xsec = -999; 
    lumi = -999;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Sample::Sample(TString ifilename, TString iname)
{
    filename = ifilename;
    name = iname;
    treename = TString("dimuons/tree");
    file = new TFile(filename);
    tree = (TTree*)file->Get(treename);
    N = tree->GetEntries();

    lumiWeights = 0;
    xsec = -999; 
    lumi = -999;
    
    calculateNoriginal();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Sample::Sample(TString ifilename, TString iname, TString isampleType)
{
    filename = ifilename;
    name = iname;
    sampleType = isampleType;
    treename = TString("dimuons/tree");
    file = new TFile(filename);
    tree = (TTree*)file->Get(treename);
    N = tree->GetEntries();

    lumiWeights = 0;
    xsec = -999; 
    lumi = -999;
    
    calculateNoriginal();
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

Sample::~Sample() {
  // free pointed to memory!
  if (tree != 0) {
    delete tree;
  }
  if (file !=0) {
    delete file;
  }
  if (lumiWeights !=0) {
    delete lumiWeights;
  }
}

///////////////////////////////////////////////////////////////////////////
// _______________________Other Functions________________________________//
///////////////////////////////////////////////////////////////////////////

void Sample::setBranchAddresses(int whichCategories)
{
      // Link only to the branches we need to save a lot of time
      // run 1 category info 
      branches.recoDimuCands  = tree->GetBranch("pairs");
      branches.recoMuons      = tree->GetBranch("muons");
      branches.jets           = tree->GetBranch("jets");
      branches.mht            = tree->GetBranch("mht");
      branches.eventInfo      = tree->GetBranch("event");
      branches.nVertices      = tree->GetBranch("nVertices");

      branches.recoDimuCands->SetAddress(&vars.recoDimuCands);
      branches.recoMuons->SetAddress(&vars.recoMuons);
      branches.jets->SetAddress(&vars.jets);
      branches.mht->SetAddress(&vars.mht);
      branches.eventInfo->SetAddress(&vars.eventInfo);
      branches.nVertices->SetAddress(&vars.nVertices);

      // extra branches needed for run 2 categories
      if(whichCategories == 2)
      {    
          branches.recoElectrons = tree->GetBranch("eles");
          branches.recoElectrons->SetAddress(&vars.recoElectrons);
      }

      // extra branches needed for MC samples
      if(!sampleType.EqualTo("data"))
      {
          branches.nPU     = tree->GetBranch("nPU");
          branches.gen_wgt = tree->GetBranch("GEN_wgt");
          branches.pu_wgt  = tree->GetBranch("PU_wgt");
          branches.eff_wgt = tree->GetBranch("IsoMu_eff_3");

          branches.gen_wgt->SetAddress(&vars.gen_wgt);
          branches.nPU->SetAddress(&vars.nPU);
          branches.pu_wgt->SetAddress(&vars.pu_wgt);
          branches.eff_wgt->SetAddress(&vars.eff_wgt);
      }
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void Sample::calculateNoriginal()
{
// Calculate the number of original events using the meta data tree
// Calculate the weighted number of original events as well
    TTree* metadata = (TTree*)file->Get("dimuons/metadata");
    metadata->Draw("sumEventWeights>>eweights_"+name);
    TH1F* sumEventWeightsHist = (TH1F*) gDirectory->Get("eweights_"+name); 

    // There are many ttrees combined so the histogram has a numEvents entry for each
    // file. The total number is equal to the average number times the total number of entries.
    nOriginalWeighted = sumEventWeightsHist->GetEntries()*sumEventWeightsHist->GetMean();

    metadata->Draw("originalNumEvents>>nevents_"+name);
    TH1F* nEventsHist = (TH1F*) gDirectory->Get("nevents_"+name); 
    nOriginal = nEventsHist->GetEntries()*nEventsHist->GetMean();
    if(metadata !=0) delete metadata;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

int Sample::getEntry(int i)
{
    if(branches.nVertices != 0) branches.nVertices->GetEntry(i);
    if(branches.eventInfo != 0) branches.eventInfo->GetEntry(i);
    if(branches.mht != 0) branches.mht->GetEntry(i);

    if(branches.recoDimuCands != 0) branches.recoDimuCands->GetEntry(i);
    if(branches.recoMuons != 0) branches.recoMuons->GetEntry(i);
    if(branches.recoElectrons != 0) branches.recoElectrons->GetEntry(i);
    if(branches.jets != 0) branches.jets->GetEntry(i);

    if(branches.eff_wgt != 0) branches.eff_wgt->GetEntry(i);
    if(branches.pu_wgt != 0) branches.pu_wgt->GetEntry(i);
    if(branches.nPU != 0) branches.nPU->GetEntry(i);
    if(branches.gen_wgt != 0) branches.gen_wgt->GetEntry(i);

    return i;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

int Sample::getEntry(int i, TEntryList* list)
{
    int treenum = list->GetEntry(i);
    tree->GetEntry(treenum);
    return treenum;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

float Sample::getWeight()
{
// Assumes getEntry has already been called to load the appropriate values.
// Gets the weight for the histogram depending on the sample type 
    if(sampleType.EqualTo("data")) return 1.0;
    else if(lumiWeights == 0) return 1.0*vars.gen_wgt*vars.pu_wgt*vars.eff_wgt;
    else return 1.0*vars.gen_wgt*vars.pu_wgt*lumiWeights->weight(vars.nPU);
    //else return 1.0*vars.gen_wgt;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

float Sample::getScaleFactor(float luminosity)
{
// Scale the MC histograms based upon the data luminosity, the number of events
// that the CMSSW analyzer looked at, and the xsec for the process
    if(sampleType.EqualTo("data")) return 1.0;
    else return luminosity*xsec/nOriginalWeighted;
}
