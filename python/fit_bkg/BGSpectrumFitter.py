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
import math as math
import numpy as np
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

categories = ['c12', 'c11', 'c10', 'c9', 'c8', 'c7', 'c6', 'c5', 'c4', 'c3', 'c2', 'c1', 'c0', 'root']
#categories = ['c12']

filedir = '/home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/'
#filename = 'validate_blinded_dimu_mass_Roch_110_160_categories3_tree_categorization_final_36814_dyAMC_minpt10.root';
filename = 'validate_UNBLINDED_dimu_mass_Roch_110_160_categories3_tree_categorization_final_36814_dyAMC_minpt10.root';

# fit values at 125 GeV for each category
bwzredux125 = [] 
bwzgamma125 = []
bernstein125 = []

max_diff_bwzr = []
max_diff_bwzg = []
max_diff_bern = []

order = 6 # order for bernstein poly

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
    bernstein_model, bernstein_params   = pdfs.bernstein(x, order=order)
    
    fr = wdm.fit(histo, bwzr_model, x, blinded=False,      save=False)
    fg = wdm.fit(histo, bwzg_model, x, blinded=False,      save=False)
    fb = wdm.fit(histo, bernstein_model, x, blinded=False, save=False)

    print "%s at 125 GeV = %f " % (fr.GetName(), fr.Eval(125))
    print "%s at 125 GeV = %f " % (fg.GetName(), fg.Eval(125))
    print "%s at 125 GeV = %f " % (fb.GetName(), fb.Eval(125))

    fr.SetLineColor(2)
    fg.SetLineColor(3)
    fb.SetLineColor(4)

    c = TCanvas("%s" % category, "%s" % category)
    l = TLegend(0.58, 0.67, 0.89, 0.89, "", "brNDC");

    # Draw data histogram
    histo.SetTitle(category)
    histo.Draw()

    # Draw bg fits
    fr.Draw("SAME")
    fg.Draw("SAME")
    fb.Draw("SAME")

    # Add bg fits to legend
    l.AddEntry(fr, fr.GetName(), "l")
    l.AddEntry(fg, fg.GetName(), "l")
    l.AddEntry(fb, fb.GetName(), "l")

    bwzredux125.append(fr.Eval(125))
    bwzgamma125.append(fg.Eval(125))
    bernstein125.append(fb.Eval(125))

    max_diff_bwzr.append(max(abs(fr.Eval(125) - fg.Eval(125)), abs(fr.Eval(125) - fb.Eval(125))))
    max_diff_bwzg.append(max(abs(fg.Eval(125) - fr.Eval(125)), abs(fg.Eval(125) - fb.Eval(125))))
    max_diff_bern.append(max(abs(fb.Eval(125) - fg.Eval(125)), abs(fb.Eval(125) - fr.Eval(125))))

    # Put fit values @ 125 GeV on canvas
    tr = TLatex(.4,.60,"%s(125 GeV) = %7.3f" % (fr.GetName(), fr.Eval(125))) 
    tg = TLatex(.4,.55,"%s(125 GeV) = %7.3f" % (fg.GetName(), fg.Eval(125)))
    tb = TLatex(.4,.50,"%s(125 GeV) = %7.3f" % (fb.GetName(), fb.Eval(125)))  
    tr.SetNDC(kTRUE);
    tg.SetNDC(kTRUE);
    tb.SetNDC(kTRUE);
    tr.Draw();
    tg.Draw();
    tb.Draw();

    l.Draw("SAME");
    c.SaveAs("fits_%s.root" % category)

print '\n=========== Fits at 125 GeV ==============\n'

print '{:<15}'.format('category') + '{:<15}'.format('bwz_redux_125') + '{:<15}'.format('bwzg_125') + '{:<15}'.format('bernstein%d_125' % order) 
for i in range(0,len(categories)):
    print ('{:<15}'.format(categories[i]) + '{:<15.2f}'.format(bwzredux125[i]) + '{:<15.2f}'.format(bwzgamma125[i]) + 
          '{:<15.2f}'.format(bernstein125[i]))

print '\n=========== max(bias)/sqrt(b) at 125 GeV ==============\n'

print '{:<15}'.format('category') + '{:<15}'.format('bwz_redux_125') + '{:<15}'.format('bwzg_125') + '{:<15}'.format('bernstein%d_125' % order) 
for i in range(0,len(categories)):
    b = np.mean([bwzredux125[i], bwzgamma125[i], bernstein125[i]])
    unc = math.sqrt(b)
    print ('{:<15}'.format(categories[i]) + '{:<15.2f}'.format(max_diff_bwzr[i]/unc) + '{:<15.2f}'.format(max_diff_bwzg[i]/unc) + 
          '{:<15.2f}'.format(max_diff_bern[i]/unc))

