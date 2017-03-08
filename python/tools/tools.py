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

colors = [58, 2, 8, 36, 91, 46, 50, 30, 9, 29, 3, 42, 98, 62, 74, 20, 29, 32, 49, 12, 3, 91]

def overlay(plot_list, title="overlay", savename="", name="", xtitle="", ytitle="", yrange=(-999,-999)):
    """ Overlay histograms or tgraphs onto the same canvas. plot_list = [('plotfile1.root', 'plotname1'), ('plotfile2.root', 'plotname2'), ... ]"""
    if savename == "": savename = title
    c = TCanvas()
    for i,plot in enumerate(plot_list):
        plot.SetLineColor(i)
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

##################################################################
#-----------------------------------------------------------------
##################################################################

def stack(stack_list, title="stack", log=True, name="stack", xtitle="Mass (GeV)", ytitle="", yrange=(-999,-999), ldim=[0.71,0.82,1.0,1.0]):
    """ Stack some histograms to compare to another histogram. Take the ratio chist/stack and put it in the pad in the bottom. """
    stack = THStack() 
    stack.SetTitle(title)
    stack.SetName(name)
    legend = TLegend(ldim[0], ldim[1], ldim[2], ldim[3], "", "brNDC")

    # keep track of these to auto zoom the y-axis
    minmax = 999999999

    for i,hist in enumerate(stack_list):
        print "%s: %f" % (hist.GetName(), hist.Integral())
        #  set the colors, etc
        hist.SetLineColor(colors[i])
        hist.SetFillColor(colors[i])
        hist.SetStats(0)

        if(i==0):
            minmax = hist.GetMaximum();
        if(hist.GetMaximum() < minmax):  minmax = hist.GetMaximum();   # the lowest maximum of the histograms

        stack.Add(hist) 
        legend.AddEntry(hist, hist.GetTitle(), "f")

    print  "stack integral: %f" % stack.GetStack().Last().Integral()
    return stack, legend, minmax;

##################################################################
#-----------------------------------------------------------------
##################################################################

def add(histo_list, name="Net_MC", title="Net MC"):
    """ Sum the histos given in the list of (file, name) pairs. """
    hsum = 0
    for i,hist in enumerate(histo_list):
        print "%s: %f" % (hist.GetName(), hist.Integral())
        if(i==0):
            hsum = hist.Clone(name) 
            hsum.SetTitle(title)
        
        else:
            hsum.Add(hist)

    print "%s: %f" % (hsum.GetName(), hsum.Integral())
    return hsum

##################################################################
#-----------------------------------------------------------------
##################################################################

def get(item_list):
    """ Return a list of ROOT objects given a list of (file, name) pairs. """
    rlist = []
    for i,item_info in enumerate(item_list):
        # first get the histograms from the rootfiles
        filename, itemname = item_info
        tfile = TFile(filename)
        gROOT.cd() # AHHHHH
        item = tfile.Get(itemname).Clone()
        print item.GetName()
        rlist.append(item)
    return rlist

##################################################################
#-----------------------------------------------------------------
##################################################################

def ratio(histo_list):
    """ Make a ratio histogram by dividing the Nth histogram and the sum of the 0 to N-1 histograms. """
    hsum = add(histo_list[0:-1])
    hlast = histo_list[-1].Clone()
    print "%s: %f, %s: %f, ratio: %f" % (hlast.GetName(), hlast.Integral(), hsum.GetName(), hsum.Integral(), hlast.Integral()/hsum.Integral())
    hlast.Divide(hsum)
    return hlast



##################################################################
#-----------------------------------------------------------------
##################################################################

def stackAndRatio(histo_list, title="stack", log=True, name="stack", xtitle="Mass (GeV)", ytitle="", yrange=(-999,-999), ytitleratio="Data/MC", ldim=[0.71,0.82,1.0,1.0]):

    c = TCanvas()
    pad1 = TPad("pad1", "pad1", 0,0.3,1,1.0)
    pad1.SetBottomMargin(0.025)
    pad1.Draw()
    pad1.cd()

    s, legend, minmax = stack(histo_list[0:-1], title=title, log=log, name=name, xtitle=xtitle, ytitle=ytitle, yrange=yrange, ldim=ldim)
    
    s.Draw("hist")
    s.GetXaxis().SetTitle("")
    s.GetYaxis().SetTitle(ytitle)
    s.GetYaxis().SetLabelFont(1);
    
    # overlay the mc error band
    sumhist = s.GetStack().Last().Clone()
    sumhist.SetFillStyle(3001)
    sumhist.SetLineColor(12)
    sumhist.SetFillColor(12)
    sumhist.Draw("E2SAME")
    legend.AddEntry(sumhist, "MC Error Band", "f")
    gPad.Modified()
    
    if(yrange==(-999,-999)): yrange = (minmax*0.1, s.GetMaximum()*10)
    if(log and yrange[0] > 0): 
        gPad.SetLogy(1)
    s.SetMaximum(yrange[1])
    s.SetMinimum(yrange[0])
    
    hlast = histo_list[-1]
    hlast.SetMarkerStyle(20)
    hlast.SetLineColor(1)
    hlast.SetFillColor(0)
    legend.AddEntry(hlast, hlast.GetTitle(), "p")
    hlast.Draw("epsames")

    # change to bottom pad and draw the ratio plot
    c.cd()
    pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetGridy()
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.2)
    pad2.Draw()
    pad2.cd()
    hratio = ratio(histo_list)
    hratio.SetMinimum(0.58)
    hratio.SetMaximum(1.42)
    hratio.SetMarkerStyle(20)
    hratio.Draw("ep")

    hratio.SetTitle("")
    hratio.GetYaxis().SetTitle(ytitleratio);
    hratio.GetYaxis().SetNdivisions(505);
    hratio.GetYaxis().SetTitleSize(20);
    hratio.GetYaxis().SetTitleFont(43);
    hratio.GetYaxis().SetTitleOffset(0.6*1.55);
    hratio.GetYaxis().SetLabelFont(43); 
    hratio.GetYaxis().SetLabelSize(15);

    hratio.GetXaxis().SetTitle(xtitle);
    hratio.GetXaxis().SetTitleSize(20);
    hratio.GetXaxis().SetTitleFont(43);
    hratio.GetXaxis().SetTitleOffset(0.8*4.);
    hratio.GetXaxis().SetLabelFont(43); 
    hratio.GetXaxis().SetLabelSize(15);


    c.SaveAs("stack.png")

