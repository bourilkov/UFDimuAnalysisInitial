##############################################
# tools.py                                   #
##############################################

#============================================
# import
#============================================

import string
import re
import argparse
import numpy as np
from ROOT import *
from array import *

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

def stackAndRatio(histo_list, rlist=0, title="stack", log=True, name="stack", xtitle="Mass (GeV)", ytitle="", yrange=(-999,-999), ytitleratio="Data/MC", ldim=[0.71,0.82,1.0,1.0]):

    if(rlist==0):
        rlist = histo_list

    for h in histo_list:
        h.SetStats(False)

    for h in rlist:
        h.SetStats(False)

    c = TCanvas()
    pad1 = TPad("pad1", "pad1", 0,0.3,1,1.0)
    pad1.SetBottomMargin(0.0125)
    pad1.Draw()
    pad1.cd()

    s, legend, minmax = stack(histo_list[0:-1], title=title, log=log, name=name, xtitle=xtitle, ytitle=ytitle, yrange=yrange, ldim=ldim)
    
    s.Draw("hist")
    s.GetXaxis().SetTitle("")
    s.GetYaxis().SetTitle(ytitle)
    s.GetYaxis().SetLabelFont(43); 
    s.GetYaxis().SetLabelSize(15);
    s.GetXaxis().SetLabelSize(0.000001);
    
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
    legend.Draw("SAME");

    # change to bottom pad and draw the ratio plot
    c.cd()
    pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetGridy()
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.2)
    pad2.Draw()
    pad2.cd()
    hratio = ratio(rlist)
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

##################################################################
#-----------------------------------------------------------------
##################################################################

def ratioError2(numerator, numeratorError2, denominator, denominatorError2):
    """ Propogation of errors for N_numerator/N_denominator
        Returns error squared, takes in numerator error squared and denom err squared as well as N_num and N_denom"""
    if(numerator   <= 0): return 0
    if(denominator <= 0): return 0
    ratio = numerator/denominator;

    return  ratio*ratio*( numeratorError2/numerator/numerator + denominatorError2/denominator/denominator );

##################################################################
#-----------------------------------------------------------------
##################################################################

def getRebinEdges(numerator, denominator, max_err = 0.10):
    """The ratio plots are a bit crazy with huge errors sometimes, so we want to rebin with variable binning
       such that the error is always low in each of the ratio plot bins"""

    newBins = []    # will return this as the new binning scheme based upon the ratio plot errors
    errVec = []     # error in each bin
    numName = numerator.GetName()
    isMassData = False
    if("dimu_mass" in numerator.GetXaxis().GetTitle() and "Data" in numName): isMassData = True
    print "numerator name: %s, numerator xtitle: %s, isMassData: %d" %(numName, numerator.GetXaxis().GetTitle(), isMassData)

    # check if the two histos are binned the same way
    compatible = numerator.GetNbinsX() == denominator.GetNbinsX() and numerator.GetBinLowEdge(1) == denominator.GetBinLowEdge(1) \
                 and numerator.GetBinLowEdge(numerator.GetNbinsX()) == denominator.GetBinLowEdge(denominator.GetNbinsX())
    if(not compatible):
        print "!!! ERROR: tools.getBinning >> numerator and denominator histograms do not have the same binning scheme!"
        return []

    if(numerator.Integral() <= 0 or denominator.Integral() <=0):
        print "!!! ERROR: tools.getBinning >>  numerator or denominator integral <= 0"
        return []

    # first low edge is the lowest edge by default
    newBins.append(numerator.GetBinLowEdge(1))     # edges for the new binning scheme
    sumNum = 0                                     
    sumErr2Num = 0
    sumDenom = 0
    sumErr2Denom = 0
    lastBinIsCutBin = 0

    # Collect bins together until the error is low enough for the corresponding ratio plot bin.
    # Once the error is low enough call that the new bin, move on and repeat.
    for i in range(1,numerator.GetNbinsX()+1):
        sumNum += numerator.GetBinContent(i)
        sumErr2Num += numerator.GetBinError(i)*numerator.GetBinError(i)

        sumDenom += denominator.GetBinContent(i)
        sumErr2Denom += denominator.GetBinError(i)*denominator.GetBinError(i)

        # squared error for the ratio, using propogation of errors
        # 0 if sumNum or sumDenom is 0
        ratioErr2 = ratioError2(sumNum, sumErr2Num, sumDenom, sumErr2Denom)

        percentError = 0
        if(ratioErr2>0): 
            percentError = TMath.Sqrt(ratioErr2)/(sumNum/sumDenom)

        # If the bins suddenly drop to zero or go up from zero then this is probably due to some cut that was applied
        # and we dont' want to combine the 0 bins on one side of the cut with the bins after the cut. It messes up the comparison plots.
        # it's better to see exactly where the cut caused everything to drop to zero.
        if(i!=numerator.GetNbinsX() and TMath.Abs(numerator.GetBinContent(i)-numerator.GetBinContent(i+1)) > 10
           and (numerator.GetBinContent(i) == 0 or numerator.GetBinContent(i+1) == 0)):
            newBins.append(numerator.GetBinLowEdge(i)+numerator.GetBinWidth(i))
            errVec.append(percentError)

            sumNum = 0
            sumErr2Num = 0

            sumDenom = 0
            sumErr2Denom = 0
            lastBinIsCutBin = True
        # we have the minimum error necessary create the bin in the vector
        # or just make a bin if we are in the blinded region of the mass spectrum
        elif(percentError !=0 and percentError <= max_err or (isMassData and numerator.GetBinLowEdge(i) >= 120 
             and numerator.GetBinLowEdge(i) < 130)):
            newBins.append(numerator.GetBinLowEdge(i)+numerator.GetBinWidth(i))
            errVec.append(percentError)

            sumNum = 0
            sumErr2Num = 0

            sumDenom = 0
            sumErr2Denom = 0

            lastBinIsCutBin = False
        # we got to the end of the histogram and the last bins can't add up
        # to get the error low enough, so we merge these last bins with the 
        # last bin set up in the new variable binning scheme.
        if(i==numerator.GetNbinsX() and (percentError > max_err or percentError == 0)):
            # get rid of the last bin added, but make sure there is at least one bin. Don't remove the low edge of the zero-th bin.
            # la la we need at least a value for the beginning and end, size 1 vector doesn't make sense and rebinning will fail.
            if(len(newBins) > 1 and not lastBinIsCutBin): 
                del newBins[-1]
            # replace with the end bin value in the numerator histo
            newBins.append(numerator.GetBinLowEdge(i)+numerator.GetBinWidth(i))
    return newBins

##################################################################
#-----------------------------------------------------------------
##################################################################

def rebin(hlist, edges):
    """Use the a variable binning scheme to rebin the histograms in the list."""
    rebin_hlist = []
    for h in hlist:
        hrebin = h.Rebin(len(edges)-1, h.GetName()+"_", array('d',edges))
        rebin_hlist.append(hrebin)
        
    return rebin_hlist

##################################################################
#-----------------------------------------------------------------
##################################################################

def scaleByBinWidth(hlist, edges):
    min_diff = 9999999999
    scale = []
    retlist = []

    # find minimum bin width
    for i in range(len(edges)-1):
        scale.append(edges[i+1]-edges[i])   # bin width for each bin
        if(edges[i+1] - edges[i] < min_diff): min_diff = edges[i+1]-edges[i]

    scale = np.array(scale)
    scale = scale/min_diff

    for h in hlist:
        hclone = h.Clone(h.GetName()+"_")
        # scale each bin in hrebin so that it is reduced by binwidth/minwidth
        for i,b in enumerate(scale):
            hclone.SetBinContent(i+1, 1/b*hclone.GetBinContent(i+1))
            hclone.SetBinError(i+1, 1/b*hclone.GetBinError(i+1))
        retlist.append(hclone)

    return retlist
