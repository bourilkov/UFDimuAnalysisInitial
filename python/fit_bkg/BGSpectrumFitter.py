##############################################
# BGSpectrumFitter.py                      #
##############################################

#============================================
# import
#============================================

import PDFDatabase as pdfs
import prettytable
import string
import re
import argparse
from ROOT import *

import sys
sys.argv.append( '-b-' )

#============================================
# code
#============================================

class BGSpectrumFitter:
# object to make workspace, root files, and datacards needed for analytic shape or template
# limit setting via higgs combine.

    infilename = ''
    category = ''
    data_hist = 0
    bg_dy_hist = 0
    bg_ttbar_hist = 0
    bg_diboson_hist = 0
    bg_not_dy_hist = 0
    bg_all_hist = 0
    tfile = 0
    nuisance_params = []

    def __init__(self, infilename, category):
        self.infilename = infilename
        self.category = category
        self.tfile = TFile(infilename)
        self.setBGHists()
    
    def setBGHists(self):
        self.data_hist      = self.tfile.Get('net_histos/'+self.category+"_Net_Data")
        self.bg_dy_hist      = self.tfile.Get('net_histos/'+self.category+"_Drell_Yan")
        self.bg_ttbar_hist   = self.tfile.Get('net_histos/'+self.category+"_TTbar_Plus_SingleTop")
        self.bg_diboson_hist = self.tfile.Get('net_histos/'+self.category+"_Diboson_plus")
        self.bg_all_hist  = self.tfile.Get('net_histos/'+self.category+"_Net_Bkg")

        self.bg_not_dy_hist  = self.bg_ttbar_hist.Clone()
        self.bg_not_dy_hist.Add(self.bg_diboson_hist)
        self.bg_not_dy_hist.SetName(self.category+"_TTBar_Diboson_plus")
        self.bg_not_dy_hist.SetTitle(self.category+"_TTBar_Diboson_plus")
    
    
    def getX(self, histo):
        # suppress all messages except those that matter
        RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
        print "="*78
    
        # Get the dimu mass histogram to use for limit setting
        nbins   = histo.GetNbinsX()
        massmin = histo.GetBinLowEdge(1)
        massmax = massmin + nbins*histo.GetBinWidth(1)
    
        #----------------------------------------
        # create di-muon mass variable
        # syntax:
        # <name>[initial-val, min-val, max-val]
        #----------------------------------------
        x = RooRealVar('x','x',120, massmin, massmax)
        x.SetTitle('m_{#mu#mu}')
        x.setUnit('GeV')
  
        return x
    
    def fitAndSave(self, histo, pdfMmumu, x, xmin=110, xmax=310):
    # make workspace with signal model and background model for analytic shape fit.
    # save it to a root file.

        # suppress all messages except those that matter
        RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
        print "="*78

        # create binned dataset from histogram
        # needs to be named data_obs for higgs combine limit setting
        data = RooDataHist('data_obs', 'data_obs', RooArgList(x), histo)

        #----------------------------------------
        # fit and plot
        #----------------------------------------
        x.setRange("window",xmin-0.1, xmax+0.1)  # whole window for MC
        x.setRange("left",xmin-0.1, 120+0.1)    # exclude signal region for Data
        x.setRange("right",130-0.1, xmax+0.1)

      
        if "Data" in histo.GetName():
            pdfMmumu.fitTo(data, RooFit.Range("left,right"))
        else: 
            pdfMmumu.fitTo(data, RooFit.Range("window"))

        xframe = x.frame(RooFit.Name(histo.GetName()+"_Fit"), RooFit.Title(histo.GetName()+"_Fit"))
        xframe.GetXaxis().SetNdivisions(505)
        data.plotOn(xframe)
        pdfMmumu.plotOn(xframe)
        pdfMmumu.paramOn(xframe, RooFit.Format("NELU", RooFit.AutoPrecision(2)), RooFit.Layout(0.3, 0.95, 0.92) )
        chi2 = xframe.chiSquare(2)

        print "chi2    :     %7.3f"               % chi2
        print

        c1 = TCanvas(histo.GetName()+"_"+pdfMmumu.GetName()+"_c", histo.GetName()+"_"+pdfMmumu.GetName(), 10, 10, 600, 600)
        xframe.Draw()
        t = TLatex(.6,.6,"#chi^{2}/ndof = %7.3f" % chi2);  
        t.SetNDC(kTRUE);
        t.Draw();
        c1.SaveAs('img/'+c1.GetName()+'.png')
        c1.SaveAs(histo.GetName()+"_"+pdfMmumu.GetName()+'.root');

#----------------------------------------
# Let's fit some backgrounds using the
# object
#----------------------------------------

print('program is running ...')
#category = 'c_01_Jet_Tight_OE' 
#category = 'c_01_Jet_Tight_OO' 
#category = 'c_2_Jet_VBF_Loose' 
#category = 'c_2_Jet_GGF_Tight' 

categories = ['c_01_Jet_Tight_OE', 'c_01_Jet_Tight_BB', 'c_2_Jet_VBF_Loose', 'c_2_Jet_VBF_Tight', 'c_2_Jet_GGF_Tight']


for category in categories:
    wdm = BGSpectrumFitter('/home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/validate_blinded_dimu_mass_PF_110_310_nolow_run1categories_36814.root', category) 
    print wdm.infilename, wdm.category
    
    #----------------------------------------
    # Set up our x variable and get the histo
    # we want to fit
    #----------------------------------------
    
    #histo = wdm.bg_dy_hist
    #histo = wdm.bg_all_hist
    #histo = wdm.bg_not_dy_hist
    histo = wdm.data_hist
    
    x = wdm.getX(histo)
    
    lin_model, lin_params     = pdfs.linear(x)
    hgg_model, hgg_params     = pdfs.higgsGammaGamma(x)
    bwzr_model, bwzr_params   = pdfs.bwZredux(x)

    bwzg_model = 0
    bwzg_params = 0
    if "2_Jet" in category:
        bwzg_model, bwzg_params = pdfs.bwZGamma(x)
    else:
        bwzg_model, bwzg_params = pdfs.bwZGamma(x)

    #bwz_model, bwz_params      = pdfs.bwZ(x)
    #bwg_model, bwg_params     = pdfs.bwGamma(x)
    #bwzl_model, bwzl_params   = pdfs.bwZPlusLinear(x)
    #cheby_model, cheby_params = pdfs.chebychev(x)
    #bwzgl_model, bwzgl_params = pdfs.bwZGammaPlusLinear(x)
    
    #mix = RooRealVar("mix_exp_lin","mix",0.1,0,1)
    #hgg_lin_model = RooAddPdf("hggexp_plus_linear","hggexp_plus_linear", RooArgList(hgg_model,lin_model),RooArgList(mix))
    #fewz_lin_model = RooAddPdf("fewz_plus_linear","fewz_plus_linear", RooArgList(fewz_model,lin_model),RooArgList(mix))
    
    #model = lin_model
    #model = bwzg_model
    #model = bwzr_model
    model = hgg_model
    
    wdm.fitAndSave(histo, model, x, xmax=310)
