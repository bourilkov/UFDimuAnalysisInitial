// Get the FEWZ dimu mass histos made by the python script (..python/fewz)
// Get the DY_MC dimu mass histograms made by fewzCompare.cxx
// Overlay and make a ratio plot to compare them. Do this for every category.

#include "Sample.h"
#include "DiMuPlottingSystem.h"
#include "CutSet.h"
#include "Cut.h"
#include "SelectionCuts.h"
#include "CategorySelection.h"
#include "JetSelectionTools.h"

#include "ParticleTools.h"
#include "EventTools.h"
#include "PUTools.h"
#include "SignificanceMetrics.cxx"

#include "TLorentzVector.h"

#include <sstream>
#include <map>
#include <vector>
#include <utility>

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

TH1F* getFewzTH1F(TString categoryName)
{
    std::cout << "getFEWZ looking at category, " << categoryName << std::endl;
    TString dir = "../python/fewz/fewz_predictions/";
    TH1F* fewzHisto = 0;
    TString fewzHistoName = "Q_llInvaria";

    if(categoryName.EqualTo("1Jet_Narrow"))
        fewzHisto = (TH1F*) TFile(dir+"NLO.out_mu13tev_higgs_PI_MRST2004qed_dimi_NLO_WNLO_110-160_1jet.root").Get(fewzHistoName);

    else if(categoryName.EqualTo("1Jet_Wide"))
        fewzHisto = (TH1F*) TFile(dir+"NLO.out_mu13tev_higgs_PI_MRST2004qed_dimi_NLO_WNLO_110-310_1jet.root").Get(fewzHistoName);

    else if(categoryName.EqualTo("Central_Central_Narrow"))
        fewzHisto = (TH1F*) TFile(dir+"NLO.out_mu13tev_higgs_PI_MRST2004qed_dimi_NLO_WNLO_110-160_cc.root").Get(fewzHistoName);

    else if(categoryName.EqualTo("Central_Central_Wide"))
        fewzHisto = (TH1F*) TFile(dir+"NLO.out_mu13tev_higgs_PI_MRST2004qed_dimi_NLO_WNLO_110-310_cc.root").Get(fewzHistoName);

    else if(categoryName.EqualTo("Narrow"))
        fewzHisto = (TH1F*) TFile(dir+"NLO.out_mu13tev_higgs_PI_MRST2004qed_dimi_NLO_WNLO_110-160_full.root").Get(fewzHistoName);

    else if(categoryName.EqualTo("Wide"))
        fewzHisto = (TH1F*) TFile(dir+"NLO.out_mu13tev_higgs_PI_MRST2004qed_dimi_NLO_WNLO_110-310_full.root").Get(fewzHistoName);

    else if(categoryName.EqualTo("Central_Not_Central_Narrow"))
    {
        TH1F* cc = (TH1F*) TFile(dir+"NLO.out_mu13tev_higgs_PI_MRST2004qed_dimi_NLO_WNLO_110-160_cc.root").Get(fewzHistoName);
        TH1F* ncnc = (TH1F*) TFile(dir+"NLO.out_mu13tev_higgs_PI_MRST2004qed_dimi_NLO_WNLO_110-160_ncnc.root").Get(fewzHistoName);
        TH1F* full = (TH1F*) TFile(dir+"NLO.out_mu13tev_higgs_PI_MRST2004qed_dimi_NLO_WNLO_110-160_full.root").Get(fewzHistoName);

        // central not central = full - not central not central - central central
        ncnc->Add(cc);
        full->Add(ncnc,-1);
        fewzHisto = full;
        delete cc;
        delete ncnc;
    }
    else if(categoryName.EqualTo("Central_Not_Central_Wide"))
    {
        TH1F* cc = (TH1F*) TFile(dir+"NLO.out_mu13tev_higgs_PI_MRST2004qed_dimi_NLO_WNLO_110-310_cc.root").Get(fewzHistoName);
        TH1F* ncnc = (TH1F*) TFile(dir+"NLO.out_mu13tev_higgs_PI_MRST2004qed_dimi_NLO_WNLO_110-310_ncnc.root").Get(fewzHistoName);
        TH1F* full = (TH1F*) TFile(dir+"NLO.out_mu13tev_higgs_PI_MRST2004qed_dimi_NLO_WNLO_110-310_full.root").Get(fewzHistoName);

        // central not central = full - not central not central - central central
        ncnc->Add(cc);
        full->Add(ncnc,-1);
        fewzHisto = full;
        fewzHisto = full;
        delete cc;
        delete ncnc;
    }

    std::cout << "fewzHistoName, integral" << std::endl;
    if(fewzHisto != 0)
    {
        std::cout << fewzHisto->GetName() << ", " << fewzHisto->Integral() << std::endl;
        fewzHisto->SetName("fewz_"+categoryName);
        fewzHisto->SetTitle("fewz_"+categoryName);
    }
    else std::cout << "isNULL" << std::endl;
    std::cout << std::endl;

    return fewzHisto;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

TH1F* getSampleTH1F(TString categoryName, TString filename, bool isDYMC)
{
    std::cout << "getSampleTH1F looking at category, " << categoryName << std::endl;
    TH1F* dyHisto = 0;
    TString dyHistoName = "histos/" + categoryName;
    if(isDYMC) dyHistoName += "_DYJetsToLL";
    else dyHistoName += "_Data";

    TFile* file = new TFile(filename);

    if(categoryName.EqualTo("1Jet_Narrow"))
        dyHisto = (TH1F*) file->Get(dyHistoName);

    else if(categoryName.EqualTo("1Jet_Wide"))
        dyHisto = (TH1F*) file->Get(dyHistoName);

    else if(categoryName.EqualTo("Central_Central_Narrow"))
        dyHisto = (TH1F*) file->Get(dyHistoName);

    else if(categoryName.EqualTo("Central_Central_Wide"))
        dyHisto = (TH1F*) file->Get(dyHistoName);

    else if(categoryName.EqualTo("Narrow"))
        dyHisto = (TH1F*) file->Get(dyHistoName);

    else if(categoryName.EqualTo("Wide"))
        dyHisto = (TH1F*) file->Get(dyHistoName);

    else if(categoryName.EqualTo("Central_Not_Central_Narrow"))
        dyHisto = (TH1F*) file->Get(dyHistoName);

    else if(categoryName.EqualTo("Central_Not_Central_Wide"))
        dyHisto = (TH1F*) file->Get(dyHistoName);

    std::cout << "dyHistoName, integral" << std::endl;
    if(dyHisto != 0) 
    {
        std::cout << dyHisto->GetName() << ", " << dyHisto->Integral() << std::endl;
        dyHisto->SetName("dy_mc_"+categoryName);
        dyHisto->SetTitle("dy_mc_"+categoryName);
    }
    else std::cout << "isNULL" << std::endl;
    std::cout << std::endl;

    delete file;

    return dyHisto;
}
//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    TH1::AddDirectory(kFALSE);

    // The filename with the DY_MC histos. Filenames for fewz histos are given in getFewzTH1F.
    TString filename = "rootfiles/00111_validate_dimu_mass_DY-FEWZ_MC_categories_3990.root";
    bool isDYMC = true;
    bool log = false;
    double luminosity = 3990;

    // get binary options for the cuts and filling from the string
    int slash = filename.First("/");
    int uscore = filename.First("_");
    TString cinfo = filename(slash+1, uscore-slash);

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) filename = TString(argv[i]);
        if(i==2) ss >> isDYMC;
        if(i==3) ss >> luminosity;
        if(i==4) ss >> log;
    }   

    TString varname = "dimu_mass";

    std::cout << "DY Histos Filename: " << filename << std::endl;
    std::cout << "luminosity        : " << luminosity << std::endl;
    std::cout << "Log?              : " << log  << std::endl;
    std::cout << std::endl;
    
    // Use this to plot some things if we wish
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();
    CategorySelectionFEWZ categorySelection;

    TList* varstacklist = new TList();  // list to save all of the stacks
    TList* histolist = new TList();     // list to save all of the histos

    TList* narrowlist = new TList();
    TList* widelist = new TList();

    // For each category we want to compare the dimu_mass for DY_MC with the dimu mass for FEWZ
    for(auto &c : categorySelection.categoryMap)
    {
        // Name the canvas 
        TString cname = c.first+"_stack";
        // Grab the fewz histogram for this category
        TH1F* fewzHisto = getFewzTH1F(c.first);
        TH1F* dyHisto = getSampleTH1F(c.first, filename, isDYMC);
        if(fewzHisto != 0 && dyHisto != 0)
        {
            c.second.histoList->Add(dyHisto);
            // if we don't care about this cateogry or it didn't work then move on
            // Get the dy histogram that we just made earlier
            //fewzHisto->Scale(dyHisto->Integral()/fewzHisto->Integral());
            fewzHisto->Scale(luminosity);

            // Make a clone with unit scaling for the overlay plots
            TH1F* fewzHistoClone = (TH1F*)fewzHisto->Clone();
            fewzHistoClone->Scale(1/fewzHistoClone->Integral());
            fewzHistoClone->SetName(c.first);

            // Add the fewzHisto to the list so that we can make a stack
            c.second.histoList->Add(fewzHisto);
            //                                              list of histos    ,name  , title, xtitle,  ytitle,      , rebin, fit, ratiotitle,  log,  stats, legend
            TCanvas* stack = dps->stackedHistogramsAndRatio(c.second.histoList, cname, cname, varname, "Num Entries", true, true, "FEWZ/DY_MC", log, false, 1);
            varstacklist->Add(stack);

            // Add the two histograms to the list
            histolist->Add(c.second.histoList);

            if(log) stack->SaveAs("../python/fewz/img/"+cinfo+cname+"_log.png");
            else stack->SaveAs("../python/fewz/img/"+cinfo+cname+".png");

            if(c.first.Contains("Wide")) widelist->Add(fewzHistoClone);
            if(c.first.Contains("Narrow")) narrowlist->Add(fewzHistoClone);
       }
    }

    TCanvas* cwideOverlay = dps->overlay(widelist, "FEWz_Overlay_Wide", "FEWz_Overlay_Wide", "dimu_mass", "xsec", log);
    std::cout << cwideOverlay->GetName() << std::endl;
    TCanvas* cnarrowOverlay = dps->overlay(narrowlist, "FEWz_Overlay_Narrow", "FEWz_Overlay_Narrow", "dimu_mass", "xsec", log);
    std::cout << cnarrowOverlay->GetName() << std::endl;

    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    filename.ReplaceAll("validate", Form("overlay_fewz_DYMC%d",(int)isDYMC));
    TFile* savefile = new TFile(filename, "RECREATE");
    TDirectory* stacks = savefile->mkdir("stacks");
    TDirectory* overlays = savefile->mkdir("overlays");
    TDirectory* histos = savefile->mkdir("histos");

    // save the different histos and stacks in the appropriate directories in the tfile
    stacks->cd();
    varstacklist->Write();

    overlays->cd();
    cwideOverlay->Write();
    cnarrowOverlay->Write();

    if(log) cwideOverlay->SaveAs("../python/fewz/img/"+cinfo+"FEWz_Compare_Overlay_Wide_log.png");
    else  cwideOverlay->SaveAs("../python/fewz/img/"+cinfo+"FEWz_Compare_Overlay_Wide.png");
    if (log) cnarrowOverlay->SaveAs("../python/fewz/img/"+cinfo+"FEWz_Compare_Overlay_Narrow_log.png");
    else cnarrowOverlay->SaveAs("../python/fewz/img/"+cinfo+"FEWz_Compare_Overlay_Narrow.png");

    histos->cd();
    histolist->Write();
    savefile->Close();

    return 0;
}
