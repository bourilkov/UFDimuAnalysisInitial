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

