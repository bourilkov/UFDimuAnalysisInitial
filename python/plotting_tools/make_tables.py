##############################################
# make_tables.py                             #
##############################################

import tools as tools
import sys
from ROOT import *
from categories_and_samples import *
import re
import argparse
import numpy as np

gROOT.SetBatch()

parser = argparse.ArgumentParser(description='Parse cmd line arguments.')
parser.add_argument('f')
opts = parser.parse_args()
filename = opts.f

#[function(number) for number in numbers if condition(number)]

#============================================
# code
#============================================

print "\n /// Making tables for %s ... \n" % filename

# Check for DY naming bug ------------------------------
dy_bug = False
tfile = TFile(filename)
tfile_dir = tfile.GetDirectory('net_histos')
for key in tfile_dir.GetListOfKeys():
    name = key.GetName()
    #print "(%d) Looping over %s ..." % (dy_bug, key)
    if 'Drell_Yan_' in name and name[-1] == '_':
        dy_bug = True
        break

# List of histograms to make the table from -----------------
in_list = []
for c in bdt_categories:
    for s in net_samples_all:
        if s == 'Drell_Yan' and dy_bug:
            s+='_'
        in_list.append((filename, 'net_histos/%s_%s' % (c, s))),

print "\n =========== GET ================== \n"
hlist = tools.get(in_list)

fwhm_widths = []
ffwhm_widths = []
fwhm_bounds = []
ffwhm_bounds = []
signal_net  = []
signal_fwhm = []
bkg_fwhm    = []

vbf_net  = []
vbf_fwhm = []
ggf_net  = []
ggf_fwhm = []
vh_net   = []
vh_fwhm  = []

bkg_ffwhm    = []
signal_ffwhm = []
signal_net_ffwhm = []
vbf_ffwhm = []
ggf_ffwhm = []
vh_ffwhm  = []

dy_fwhm           = []
ttbar_fwhm        = []
diboson_plus_fwhm = []

dy_ffwhm           = []
ttbar_ffwhm        = []
diboson_plus_ffwhm = []


for h in hlist:
    fwhm = 0
    #print "%s: %f" % (h.GetName(), h.Integral())
    if 'Net_Signal' in h.GetName():
        left,right,fwhm = tools.fwhm(h)
        fleft,fright,ffwhm = tools.fwhm_fit(h)
        print "FWHM %s: %f, %f" % (h.GetName(), fwhm, h.Integral(left, right))
        fwhm_widths.append(fwhm)
        ffwhm_widths.append(ffwhm)
        fwhm_bounds.append((left,right))
        ffwhm_bounds.append((fleft,fright))
        c = TCanvas()
        h.Draw()
        c.SaveAs(h.GetName()+".png")
        

for i,h in enumerate(hlist):
    left,right = fwhm_bounds[i/len(net_samples_all)]
    fleft,fright = ffwhm_bounds[i/len(net_samples_all)]
    if 'VBF' in h.GetName():
        vbf_net.append(h.Integral())
        vbf_fwhm.append(h.Integral(left,right))
        vbf_fwhm.append(h.Integral(fleft,fright))
    elif '_gg' in h.GetName():
        ggf_net.append(h.Integral())
        ggf_fwhm.append(h.Integral(left,right))
        ggf_ffwhm.append(h.Integral(fleft,fright))
    elif 'VH' in h.GetName():
        vh_net.append(h.Integral())
        vh_fwhm.append(h.Integral(left,right))
        vh_ffwhm.append(h.Integral(fleft,fright))
    elif 'Net_Signal' in h.GetName(): 
        signal_net.append(h.Integral())
        signal_fwhm.append(h.Integral(left,right))
        signal_ffwhm.append(h.Integral(fleft,fright))
    elif 'Drell_Yan' in h.GetName():
        dy_fwhm.append(h.Integral(left,right))
        dy_ffwhm.append(h.Integral(fleft,fright))
    elif 'TT' in h.GetName():
        ttbar_fwhm.append(h.Integral(left,right))
        ttbar_ffwhm.append(h.Integral(fleft,fright))
    elif 'Diboson' in h.GetName():
        diboson_plus_fwhm.append(h.Integral(left,right))
        diboson_plus_ffwhm.append(h.Integral(fleft,fright))
    elif 'Net_Bkg' in h.GetName():
        bkg_fwhm.append(h.Integral(left,right))
        #bkg_ffwhm.append(h.Integral(fleft,fright))
    elif 'Net_Data' in h.GetName():
        fit = tools.fit_bkg(h) 
        c = TCanvas()
        h.Draw()
        c.SaveAs(h.GetName()+".png")
        hfit = fit.GetHistogram()
        print "\nNBINS : %d vs %d \n" % (hfit.GetNbinsX(), h.GetNbinsX())
        #bkg_ffwhm.append(hfit.Integral(fleft,fright)*h.GetNbinsX()/hfit.GetNbinsX())
        bkg_ffwhm.append(fit.Integral(h.GetBinLowEdge(fleft), h.GetBinLowEdge(fright)+h.GetBinWidth(fright))/h.GetBinWidth(fright))

sensitivity_values = []
sensitivity_titles = ['category','signal_fwhm','signal_in_fwhm','bkg_in_fwhm', 's_over_b', 's_over_sqrt_b']
sensitivity_values.append(bdt_categories)
sensitivity_values.append(fwhm_widths)
sensitivity_values.append(signal_fwhm)
sensitivity_values.append(bkg_fwhm)
sensitivity_values.append(np.array(signal_fwhm)/np.array(bkg_fwhm))
sensitivity_values.append(np.array(signal_fwhm)/np.sqrt(bkg_fwhm))

sensitivity_values_fit = []
sensitivity_titles_fit = ['category','signal_fwhm','signal_in_fwhm','bkg_in_fwhm', 's_over_b', 's_over_sqrt_b']
sensitivity_values_fit.append(bdt_categories)
sensitivity_values_fit.append(ffwhm_widths)
sensitivity_values_fit.append(signal_ffwhm)
sensitivity_values_fit.append(bkg_ffwhm)
sensitivity_values_fit.append(np.array(signal_ffwhm)/np.array(bkg_ffwhm))
sensitivity_values_fit.append(np.array(signal_ffwhm)/np.sqrt(bkg_ffwhm))

signal_values = []
signal_titles = ['category','net_signal','ggf','vbf', 'vh']
signal_values.append(bdt_categories)
signal_values.append(signal_net)
signal_values.append(ggf_net)
signal_values.append(vbf_net)
signal_values.append(vh_net)


print bdt_categories
print signal_fwhm
print len(bdt_categories)
print len(signal_fwhm)

#root_VH: 10.753939
#root_H2Mu_VBF: 19.897926
#root_H2Mu_gg: 219.350811
#root_Diboson_plus: 6503.863871
#root_TTbar_Plus_SingleTop: 57263.086929
#root_Drell_Yan: 440643.825693
#root_Net_Signal: 748.973609
#root_Net_Bkg: 504410.776493
#root_Net_Data: 391503.000000

print "\n =========== SENSITIVITY TABLE ================ \n" 

width = 20
decimals = 2
titles = ''
values = ''
for i,t in enumerate(sensitivity_titles):
    titles+='{field:<{width}}'.format(field=t,width=width)

for i,c in enumerate(bdt_categories):
    for j in range(len(sensitivity_titles)): 
        if j == 0: values+='{field:<{width}}'.format(field=sensitivity_values[j][i],width=width)
        #else: values+='{field:<{width}.{decimals}}'.format(field=sensitivity_values[j][i],width=width, decimals=decimals)
        else: values+='{field:<{width}}'.format(field=sensitivity_values[j][i],width=width)
    values+='\n'

print titles
print values

print "\n =========== SENSITIVITY TABLE FITS ================ \n" 

width = 20
decimals = 2
titles = ''
values = ''
for i,t in enumerate(sensitivity_titles_fit):
    titles+='{field:<{width}}'.format(field=t,width=width)

for i,c in enumerate(bdt_categories):
    for j in range(len(sensitivity_titles_fit)): 
        if j == 0: values+='{field:<{width}}'.format(field=sensitivity_values_fit[j][i],width=width)
        #else: values+='{field:<{width}.{decimals}}'.format(field=sensitivity_values[j][i],width=width, decimals=decimals)
        else: values+='{field:<{width}}'.format(field=sensitivity_values_fit[j][i],width=width)
    values+='\n'

print titles
print values

print "\n =========== SIGNAL TABLE ================ \n" 

titles = ''
values = ''
for i,t in enumerate(signal_titles):
    titles+='{field:<{width}}'.format(field=t,width=width)

for i,c in enumerate(bdt_categories):
    for j in range(len(signal_titles)): 
        if j == 0: values+='{field:<{width}}'.format(field=signal_values[j][i],width=width)
        #else: values+='{field:<{width}.{decimals}}'.format(field=signal_values[j][i],width=width, decimals=decimals)
        else: values+='{field:<{width}}'.format(field=signal_values[j][i],width=width)
    values+='\n'

print titles
print values

#print '{:<15}'.format('category') + '{:<15}'.format('mean_roch') + '{:<15}'.format('mean_kamu') + '{:<15}'.format('sigma_roch') + '{:<15}'.format('sigma_kamu') 
#for i in range(0,len(ggf_ratios)):
#    print ('{:<15}'.format(categories[i]) + '{:<15.2f}'.format(ggf_mean_roch[i]) + '{:<15.2f}'.format(ggf_mean_kamu[i]) + 
#          '{:<15.2f}'.format(ggf_sigma_roch[i]) + '{:<15.2f}'.format(ggf_sigma_kamu[i]) )
#

