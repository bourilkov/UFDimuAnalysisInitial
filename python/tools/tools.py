##############################################
# tools.py                                   #
##############################################

#============================================
# import
#============================================

import string
import re
import argparse
from ROOT import *

import sys
#sys.argv.append( '-b-' )

#============================================
# code
#============================================

def overlay(plot_list, title="overlay", savename="", name="", xtitle="", ytitle="", yrange=(-999,-999)):
    """ Overlay histograms or tgraphs onto the same canvas. plot_list = [('plotfile1.root', 'plotname1'), ('plotfile2.root', 'plotname2'), ... ]"""
    if savename == "": savename = title
    c = TCanvas()
    for i,plot_info in enumerate(plot_list):
        filename, plotname = plot_info
        print "%s, %s" % (filename, plotname)
        tfile = TFile(filename)
        plot = tfile.Get(plotname)
        plot.SetLineColor(i+1)
        if(i==0):
            plot.Draw()
            if yrange != (-999,-999): plot.GetYaxis().SetRangeUser(yrange[0], yrange[1])
            if xtitle!="": plot.GetXaxis().SetTitle(xtitle)
            if ytitle!="": plot.GetYaxis().SetTitle(ytitle)
            c.SetGridx()
            c.SetGridy()
            c.SetTitle(title)
            plot.SetTitle(title)
            c.SetName(title)
            if name!="": c.SetName(name)
        else:
            plot.Draw("SAME")
        c.SaveAs('imgs/%s.png' % savename)
        c.SaveAs('rootfiles/%s.root' % savename)

