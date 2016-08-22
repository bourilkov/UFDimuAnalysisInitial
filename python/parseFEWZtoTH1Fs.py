##############################################
# parseFEWZtoTH1Fs.py                        #
##############################################
# read in a FEWZ data table and make         #
# TH1Fs from it.                             #
#                                            #
##############################################


#============================================
# import
#============================================

import string
import re
import argparse
from ROOT import *

#============================================
# code
#============================================

parser = argparse.ArgumentParser(description='Convert FEWZ data tables from txt file into TH1Fs and save to a root file.')
parser.add_argument("--infile", help="the name of the FEWZ txt file")
args = parser.parse_args()

infilename = ""
if args.infile:
    infilename = args.infile
else:
    infilename = 'fewz/NLO.out_mu13tev_higgs_PI_MRST2004qed_dimi_NLO_WNLO_110-160_cc.txt'

def parseFile(infilename):
    """ Parse a fewz output file into a dictionary where each key is the histo name and the value is the data. """

    file = open(infilename, 'r')
    
    athisto = False
    histtitle = ''

    # dictionary of histoname, histodata pairs
    histos = {}
    
    for line in file:
        # skip all lines until we get to a histogram
        if ' ----' not in line and athisto == False:
            continue
        # found a histogram, should be histos from here on out
        if '  ----' in line:
            athisto = True
            # get the name/title of the histogram
            histtitle = re.sub('[- ]', '', line).strip()
            print '\"'+histtitle+'\"'
            # create a new entry in the dict for this histogram
            histos[histtitle] = []
        if '\n' == line:
            continue
        # get the histogram data
        if len(line.split()) == 5 and ' ----' not in line:
            # x = [bin, weight, error]
            xstr = line.split()[0:3]
            if('undefined' in xstr[2]): xstr[2] = re.sub('undefined','-9999', xstr[2])
            x = map(float, xstr)
            print "    ",
            print x
            # add the data to the dictionary for the appropriate histogram
            histos[histtitle].append(x)

    return histos

def saveHistograms(histodict, savefilename):
    """ Take the dictionary of histo name, data pairs and save to TH1Fs in a root file. """

    # where to save the ROOT histograms
    savefile = TFile(savefilename, "RECREATE")

    for histoname, histodata in histodict.items():
        print 
        print histoname
        print "///////////////////////////////////////////"
        print histodata
        print 

        if len(histodata) <=1:
            continue

        w = (histodata[1][0] - histodata[0][0])/2
        xmin = histodata[0][0] - w                # xmin is the first xvalue minus half the bin width since FEWZ uses midpoints of bins as xvalue
        xmax = histodata[-1][0] + w               # xmax is the last xvalue plus half the bin width
        nbins = len(histodata)

        h = TH1F(histoname, histoname, nbins, xmin, xmax)

        # fill the histogram
        for datapoint in histodata:
            # AddBinContent(bin corresponding to xvalue, fillvalue)
            h.AddBinContent(h.FindBin(datapoint[0]), datapoint[1])

        # save the histogram
        h.SetStats(kFALSE)
        h.GetXaxis().SetTitle(histoname)
        h.Write()

    savefile.Write()
    savefile.Close()

savefilename = re.sub('.txt', '.root', infilename)
histodict = parseFile(infilename)
#saveHistograms(histodict, savefilename)
print 
print savefilename
print 

