##############################################
# BGSpectrumFitter.py                      #
##############################################
# Makes roofit workspace and fits fewz       #
# dimu mass spectrum.                        #
##############################################

#============================================
# import
#============================================

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
    bg_dy_hist = 0
    bg_ttbar_hist = 0
    bg_diboson_hist = 0
    bg_not_dy_hist = 0
    tfile = 0
    nuisance_params = []

    def __init__(self, infilename, category):
        self.infilename = infilename
        self.category = category
        self.tfile = TFile(infilename)
        self.setBGHists()
    
    def setBGHists(self):
        self.bg_dy_hist      = self.tfile.Get('net_histos/'+self.category+"_Drell_Yan")
        self.bg_ttbar_hist   = self.tfile.Get('net_histos/'+self.category+"_TTbar_Plus_SingleTop")
        self.bg_diboson_hist = self.tfile.Get('net_histos/'+self.category+"_Diboson_plus")
    
    def fitAndSave(self, whichFunction = 'bw'):
    # make workspace with signal model and background model for analytic shape fit.
    # save it to a root file.
        # suppress all messages except those that matter
        RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
        print "="*78
    
        # Get the dimu mass histogram to use for limit setting
        nbins   = self.bg_dy_hist.GetNbinsX()
        massmin = self.bg_dy_hist.GetBinLowEdge(1)
        massmax = massmin + nbins*self.bg_dy_hist.GetBinWidth(1)
    
        #----------------------------------------
        # create di-muon mass variable
        # syntax:
        # <name>[initial-val, min-val, max-val]
        #----------------------------------------
        x = RooRealVar('x','x',120, massmin, massmax)
        x.SetTitle('m_{#mu#mu}')
        x.setUnit('GeV')
    
        # create binned dataset from histogram
        # needs to be named data_obs for higgs combine limit setting
        data = RooDataHist('data_obs', 'data_obs', RooArgList(x), self.bg_dy_hist)

        #----------------------------------------
        # breit weigner mixture (run1 bg)
        #----------------------------------------
        bwWidth =  RooRealVar("bwWidth","bwWidth",2.5,0,30)
        bwmZ =     RooRealVar("bwmZ","bwmZ",91.2,85,95)
        expParam = RooRealVar("expParam","expParam",-1e-03,-1e-01,1e-01)
        mixParam = RooRealVar("mixParam","mixParam",0.5,0,1)

        bwWidth.setConstant(True);
        bwmZ.setConstant(True);

        phoExpMmumu = RooGenericPdf("phoExpMmumu","exp(@0*@1)*pow(@0,-2)",RooArgList(x,expParam))
        bwExpMmumu  = RooGenericPdf("bwExpMmumu","exp(@0*@3)*(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",RooArgList(x,bwmZ,bwWidth,expParam))
        bwmodel     = RooAddPdf("bwmodel","bwmodel", RooArgList(bwExpMmumu,phoExpMmumu),RooArgList(mixParam))

        #----------------------------------------
        # falling exponential (hgammgamma bg)
        #----------------------------------------
        a1 = RooRealVar("a1", "a1", 5.0, -50, 50)          # nuisance parameter1 for the background fit
        a2 = RooRealVar("a2", "a2", -1.0, -50, 50)         # nuisance parameter2 for the background fit
        one = RooRealVar("one", "one", 1.0, -10, 10) 
        one.setConstant()
    
        f = RooFormulaVar("f", "-(@1*(@0/100)+@2*(@0/100)^2)", RooArgList(x, a1, a2))
        expmodel = RooExponential('expmodel', 'expmodel', f, one) # exp(1*f(x))

        #----------------------------------------
        # choose fit function
        #----------------------------------------
        pdfMmumu = 0
        if whichFunction is 'bw' : pdfMmumu = bwmodel
        if whichFunction is 'exp': pdfMmumu = expmodel

        #----------------------------------------
        # fit and plot
        #----------------------------------------
        x.setRange("window",110-0.1, 160+0.1) ;
        fr = pdfMmumu.fitTo(data, RooFit.Range("window"))

        xframe = x.frame(RooFit.Name(self.category+"_DY_Fit"), RooFit.Title(self.category+"_DY_Fit"))
        xframe.GetXaxis().SetNdivisions(505)
        data.plotOn(xframe)
        pdfMmumu.plotOn(xframe)
        pdfMmumu.paramOn(xframe, RooFit.Format("NELU", RooFit.AutoPrecision(2)), RooFit.Layout(0.4, 0.95, 0.92) )
        chi2 = xframe.chiSquare(2)

        print "="*80
        if whichFunction is 'bw':
            print "expParam:     %7.5f +\-%-7.5f GeV" % (expParam.getVal(), expParam.getError())
            print "mixParam:     %7.5f +\-%-7.5f GeV" % (mixParam.getVal(), mixParam.getError())
        if whichFunction is 'exp':
            print "a1      :     %7.5f +\-%-7.5f GeV" % (a1.getVal(), a1.getError())
            print "a2      :     %7.5f +\-%-7.5f GeV" % (a2.getVal(), a2.getError())
        print "chi2    :     %7.3f"               % chi2
        print

        c1 = TCanvas(self.category+'_bg_dy_fit_c', self.category+'_bg_dy_fit_c', 10, 10, 600, 600)
        xframe.Draw()
        t = TLatex(.6,.75,"#chi^{2}/ndof = %7.3f" % chi2);  
        t.SetNDC(kTRUE);
        t.Draw();
        c1.SaveAs('img/'+c1.GetName()+'.png')
        c1.SaveAs(c1.GetName()+'.root');
        c1.SetLogy(kTRUE)
        c1.Update()
        c1.SaveAs('img/'+c1.GetName()+'_log.png')

print('program is running ...')
# Needs the file with the dimu_mass plots created by categorize.cxx via running ./categorize 0 1
# also needs to know the category you want to make the root file and datacard for

category = 'c_01_Jet_Tight_OE' 

wdm = BGSpectrumFitter('/home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/validate_dimu_mass_PF_50_200_x69p2_8_0_X_MC_run1categories_36814.root', category) 
print wdm.infilename, wdm.category
wdm.fitAndSave()
