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
        self.setHists()
    
    def setHists(self):
        self.data_hist      = self.tfile.Get('net_histos/'+self.category+"_Net_Data")
        self.bg_dy_hist      = self.tfile.Get('net_histos/'+self.category+"_Drell_Yan_")
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
    
    def fit(self, histo, pdfMmumu, x, xmin=110, xmax=160, blinded=True, save=True):
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

      
        if "Data" in histo.GetName() and blinded:
            pdfMmumu.fitTo(data, RooFit.Range("left,right"))
        else: 
            pdfMmumu.fitTo(data, RooFit.Range("window"))

        xframe = x.frame(RooFit.Name(histo.GetName()+"_Fit"), RooFit.Title(histo.GetName()+"_Fit"))
        xframe.GetXaxis().SetNdivisions(505)
        data.plotOn(xframe)
        pdfMmumu.plotOn(xframe, RooFit.Name(pdfMmumu.GetName()))
        pdfMmumu.paramOn(xframe, RooFit.Format("NELU", RooFit.AutoPrecision(2)), RooFit.Layout(0.3, 0.95, 0.92) )
        chi2 = xframe.chiSquare()

        print "chi2    :     %7.3f"               % chi2
        print

        c1 = TCanvas(histo.GetName()+"_"+pdfMmumu.GetName()+"_c", histo.GetName()+"_"+pdfMmumu.GetName(), 10, 10, 600, 600)
        xframe.Draw()
        t = TLatex(.6,.6,"#chi^{2}/ndof = %7.3f" % chi2);  
        t.SetNDC(kTRUE);
        t.Draw();

        f = xframe.findObject(pdfMmumu.GetName());

        if(save):
            c1.SaveAs('img/'+c1.GetName()+'.png')
            c1.SaveAs(histo.GetName()+"_"+pdfMmumu.GetName()+'.root');

        return f;

#----------------------------------------
# Let's fit some backgrounds using the
# object
#----------------------------------------

print('program is running ...')

#categories = ['c12', 'c11', 'c10', 'c9', 'c8', 'c7', 'c6', 'c5', 'c4', 'c3', 'c2', 'c1', 'c0']
categories = ['c12']

filedir = '/home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/'
#filename = 'validate_blinded_dimu_mass_Roch_110_160_categories3_tree_categorization_final_36814_dyAMC_minpt10.root';
filename = 'validate_UNBLINDED_dimu_mass_Roch_110_160_categories3_tree_categorization_final_36814_dyAMC_minpt10.root';

for category in categories:
    wdm = BGSpectrumFitter(filedir+filename, category) 
    print wdm.infilename, wdm.category
    
    #----------------------------------------
    # Set up our x variable and get the histo
    # we want to fit
    #----------------------------------------
    
    #histo = wdm.bg_dy_hist
    #histo = wdm.bg_not_dy_hist
    #histo = wdm.bg_all_hist
    histo = wdm.data_hist
    
    x = wdm.getX(histo)
    
    bwzr_model, bwzr_params   = pdfs.bwZredux(x)
    bwzg_model, bwzg_params   = pdfs.bwZGamma(x)
    bernstein3_model, bernstein3_params   = pdfs.bernstein(x, order=3)
    bernstein4_model, bernstein4_params   = pdfs.bernstein(x, order=4)
    bernstein5_model, bernstein5_params   = pdfs.bernstein(x, order=5)
    bernstein6_model, bernstein6_params   = pdfs.bernstein(x, order=6)
    
    fr = wdm.fit(histo, bwzr_model, x, blinded=False,        save=False)
    fg = wdm.fit(histo, bwzg_model, x, blinded=False,        save=False)
    fb = wdm.fit(histo, bernstein5_model, x, order=6, blinded=False, save=False)

    print "%s at 125 GeV = %f " % (fr.GetName(), fr.Eval(125))
    print "%s at 125 GeV = %f " % (fg.GetName(), fg.Eval(125))
    print "%s at 125 GeV = %f " % (fb5.GetName(), fb.Eval(125))

    fr.SetLineColor(2)
    fg.SetLineColor(3)
    fb.SetLineColor(4)

    c = TCanvas()
    histo.Draw()
    fr.Draw("SAME")
    fg.Draw("SAME")
    fb.Draw("SAME")
    c.BuildLegend()
    c.SaveAs("fits.root")
