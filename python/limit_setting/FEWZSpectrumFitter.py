##############################################
# FEWZSpectrumFitter.py                      #
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

#============================================
# code
#============================================

class FEWZSpectrumFitter:
# object to make workspace, root files, and datacards needed for analytic shape or template
# limit setting via higgs combine.

    infilename = ''
    category = ''
    fewz_hist = 0
    tfile = 0
    nuisance_params = []

    def __init__(self, infilename, category):
        self.infilename = infilename
        self.category = category
        self.tfile = TFile(infilename)
        self.setFEWZHist()
    
    def setFEWZHist(self):
        self.fewz_hist = self.tfile.Get('histos/fewz_'+self.category)
    
    def fitAndSave(self):
    # make workspace with signal model and background model for analytic shape fit.
    # save it to a root file.
        # suppress all messages except those that matter
        RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
        print "="*78
    
        # Get the dimu mass histogram to use for limit setting
        nbins   = self.fewz_hist.GetNbinsX()
        massmin = self.fewz_hist.GetBinLowEdge(1)
        massmax = massmin + nbins*self.fewz_hist.GetBinWidth(1)
    
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
        data = RooDataHist('data_obs', 'data_obs', RooArgList(x), self.fewz_hist)

        bwWidth =  RooRealVar(self.category+"_bwWidth","bwWidth",2.5,0,30)
        bwmZ =     RooRealVar(self.category+"_bwmZ","bwmZ",91.2,85,95)
        expParam = RooRealVar(self.category+"_expParam","expParam",-1e-03,-1e-01,1e-01)
        mixParam = RooRealVar(self.category+"_mixParam","mixParam",0.5,0,1)

        bwWidth.setConstant(True);
        bwmZ.setConstant(True);

        phoExpMmumu = RooGenericPdf("phoExpMmumu","exp(@0*@1)*pow(@0,-2)",RooArgList(x,expParam))
        bwExpMmumu  = RooGenericPdf("bwExpMmumu","exp(@0*@3)*(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",RooArgList(x,bwmZ,bwWidth,expParam))
        pdfMmumu =    RooAddPdf("bak","bak", RooArgList(bwExpMmumu,phoExpMmumu),RooArgList(mixParam))

        fr = pdfMmumu.fitTo(data)

        print "="*80
        print "expParam:     %7.5f +\-%-7.5f GeV" % (expParam.getVal(), expParam.getError())
        print "mixParam:     %7.5f +\-%-7.5f GeV" % (mixParam.getVal(), mixParam.getError())
        print

        xframe = x.frame()
        xframe.GetXaxis().SetNdivisions(505)
        data.plotOn(xframe)
        pdfMmumu.plotOn(xframe)
        pdfMmumu.paramOn(xframe)
        c1 = TCanvas('fig_fewz_fit', 'fit', 10, 10, 500, 500)
        xframe.Draw()
        c1.SaveAs('.png')

print('program is running ...')
# Needs the file with the dimu_mass plots created by categorize.cxx via running ./categorize 0 1
# also needs to know the category you want to make the root file and datacard for
wdm = FEWZSpectrumFitter('/home/acarnes/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/00111_overlay_fewz_dimu_mass_DY-FEWZ_MC_categories_3990.root', 'Narrow') 
print wdm.infilename, wdm.category
wdm.fitAndSave()
