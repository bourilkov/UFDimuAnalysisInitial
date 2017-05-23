##############################################
# adhoc_bias_test.py                         #
##############################################

#============================================
# import
#============================================

import PDFDatabase as pdfs
from BGSpectrumFitter import *
import prettytable
import string
import re
import argparse
import math as math
import numpy as np
from ROOT import *

import sys
sys.argv.append( '-b-' )

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
#filename = 'validate_UNBLINDED_dimu_mass_Roch_100_150_categories3_tree_categorization_final_36814_dyAMC-J_minpt10_b-4_sig-xlumi1.root'

# fit values at 125 GeV for each category
blinded = False
order = 8 # order for bernstein poly

vals_by_cat = []
diffs_by_cat = []
fitnames = []

for ic, category in enumerate(categories):
    wdm = BGSpectrumFitter(filedir+filename, category) 
    print wdm.infilename, wdm.category
    
    #----------------------------------------
    # Set up our x variable and get the histo
    # we want to fit
    #----------------------------------------
    
    histo = wdm.data_hist
    x = wdm.getX(histo)
    
    bwzr_model, bwzr_params   = pdfs.bwZreduxFixed(x)
    bwzg_model, bwzg_params   = pdfs.bwZGamma(x)
    bernstein_model, bernstein_params   = pdfs.bernstein(x, order=order)
    h2mupoly_model, h2mupoly_params   = pdfs.h2mupoly(x, order=order)

    models = [bwzr_model, bwzg_model, bernstein_model]
    #models = [bwzr_model, bwzg_model, h2mupoly_model]

    fits = []
    vals = []
    diffs = []

    for i,m in enumerate(models):
        f = wdm.fit(histo, m, x, blinded=False, save=False, xmin=110, xmax=160)
        v = f.Eval(125)
        vals.append(v)
        fits.append(f)
        if ic == 0: 
            fitnames.append(f.GetName())

    for vi in vals:
        idiff = []
        for vj in vals:
            idiff.append(abs(vi-vj))
        diffs.append(idiff)

    c = TCanvas("%s" % category, "%s" % category)
    l = TLegend(0.58, 0.67, 0.89, 0.89, "", "brNDC");

    # Draw data histogram
    histo.SetTitle(category)
    histo.Draw()

    for i,f in enumerate(fits):
        print "%s at 125 GeV = %f " % (f.GetName(), v)
        f.SetLineColor(i+1)
        f.Draw("SAME")
        l.AddEntry(f, f.GetName(), "l")
        t = TLatex(.4,.60-0.05*i,"%s(125 GeV) = %7.3f" % (f.GetName(), f.Eval(125))) 
        t.SetNDC(kTRUE);
        t.Draw();

    l.Draw("SAME");
    c.SaveAs("rootfiles/fits_%s.root" % category)

    vals_by_cat.append(vals)
    diffs_by_cat.append(diffs)

print '\n=========== Fits at 125 GeV ==============\n'

titles = '{:<25}'.format('category') 
fields = ''

for f in fitnames:
    titles+='{:<25}'.format(f)

for i,c in enumerate(categories):
    vals = vals_by_cat[i]
    fields+='{:<25}'.format(c)
    for v in vals:
        fields+='{:<25.2f}'.format(v)
    fields+='\n'


print titles
print fields

print '\n=========== max(bias)/sqrt(b) at 125 GeV ==============\n'

titles = '{:<25}'.format('category') 
fields = ''

for f in fitnames:
    titles+='{:<25}'.format(f)

for i,c in enumerate(categories):
    fields+='{:<25}'.format(c)
    vals = vals_by_cat[i]
    diffs = diffs_by_cat[i]
    b = np.mean(vals)
    unc = math.sqrt(b)
    for d in diffs:
        fields+='{:<25.2f}'.format(max(d)/unc)
    fields+='\n'


print titles
print fields

