// figure out the smearing in going from gen -> reco by making the reco/gen ratio plot for DY_MC 

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
    TString dir = "../python/fewz/";
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
        fewzHisto->GetXaxis()->SetTitle("dimu_mass");
    }
    else std::cout << "isNULL" << std::endl;
    std::cout << std::endl;

    return fewzHisto;
}

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

TH1F* getDYJetsTH1F(TString categoryName, TString filename)
{
    std::cout << "getDYJets looking at category, " << categoryName << std::endl;
    TH1F* dyHisto = 0;
    TString dyHistoName = "histos/" + categoryName + "_DYJetsToLL";
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
    // save the errors for the histogram correctly so they depend upon 
    // the number used to fill originally rather than the scaling
    TH1::SetDefaultSumw2();
    
    // Don't automatically associate the histogram with a tfile or directory
    // save it on your own
    TH1::AddDirectory(kFALSE);

    TString file_recocuts_recofill = "rootfiles/11110_validate_dimu_mass_DY-FEWZ_MC_categories_3990.root";
    TString file_gencuts_genfill = "rootfiles/00111_validate_dimu_mass_DY-FEWZ_MC_categories_3990.root";
    bool log = false;
    double luminosity = 3990;

    // get binary options for the cuts and filling from the string
    TString cinfo = "smearing_";

    for(int i=1; i<argc; i++)
    {   
        std::stringstream ss; 
        ss << argv[i];
        if(i==1) file_recocuts_recofill = TString(argv[i]);
        if(i==2) file_gencuts_genfill = TString(argv[i]);
        if(i==3) ss >> luminosity;
        if(i==4) ss >> log;
    }   

    TString varname = "dimu_mass";

    std::cout << "DY Histos reco-cut-reco-fill: " << file_recocuts_recofill << std::endl;
    std::cout << "DY Histos gen-cut-gen-fill  : " << file_gencuts_genfill << std::endl;
    std::cout << "luminosity                  : " << luminosity << std::endl;
    std::cout << "Log?                        : " << log  << std::endl;
    std::cout << std::endl;
    
    // Use this to plot some things if we wish
    DiMuPlottingSystem* dps = new DiMuPlottingSystem();
    CategorySelectionFEWZ categorySelection;

    TList* varstacklist = new TList();  // list to save all of the stacks
    TList* histolist = new TList();     // list to save all of the histos

    TList* narrowlist = new TList();
    TList* widelist = new TList();

    // For each category we want to compare the dimu_mass for DY with the dimu mass for FEWZ
    for(auto &c : categorySelection.categoryMap)
    {
        // Name the canvas 
        TString cname = c.first+"_"+varname+"_stack";
        // Grab the fewz histogram for this category
        TH1F* genHisto =  getDYJetsTH1F(c.first, file_gencuts_genfill);
        TH1F* recoHisto =  getDYJetsTH1F(c.first, file_recocuts_recofill);
        if(genHisto != 0 && recoHisto != 0)
        {
            genHisto->SetName(TString("gen_")+genHisto->GetName());
            recoHisto->SetName(TString("reco_")+genHisto->GetName());
            c.second.histoList->Add(genHisto);
            c.second.histoList->Add(recoHisto);

            //                                              list of histos    ,name  , title, xtitle,  ytitle,      , rebin, fit , ratiotitle, log, stats, legend
            TCanvas* stack = dps->stackedHistogramsAndRatio(c.second.histoList, cname, cname, varname, "Num Entries", true, true, "RECO/GEN", log, false, 1);
            varstacklist->Add(stack);

            // Add the two histograms to the list
            histolist->Add(c.second.histoList);

            if(log) stack->SaveAs("../python/fewz/img/"+cinfo+cname+"_log.png");
            else stack->SaveAs("../python/fewz/img/"+cinfo+cname+".png");
       }
    }

    std::cout << "  /// Saving plots..." << std::endl;
    std::cout << std::endl;
    TFile* savefile = new TFile("rootfiles/fewz_gen_smearing.root", "RECREATE");
    TDirectory* stacks = savefile->mkdir("stacks");
    TDirectory* histos = savefile->mkdir("histos");

    // save the different histos and stacks in the appropriate directories in the tfile
    stacks->cd();
    varstacklist->Write();

    histos->cd();
    histolist->Write();
    savefile->Close();

    return 0;
}
