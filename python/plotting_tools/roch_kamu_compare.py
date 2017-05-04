##############################################
# roch_kamu_compare.py                       #
##############################################

import tools as tools
import numpy as np
from ROOT import *

#============================================
# code
#============================================

categories = [
              'root',
              'c0',
              'c1',
              'c2',
              'c3',
              'c4',
              'c5',
              'c6',
              'c7',
              'c8',
              'c9',
              'c10',
              'c11',
              'c12'
             ]

# BKG GluGlu VBF WPlusH WMinusH ttH ZH 
signals = [
           'H2Mu_gg',
           'H2Mu_VBF',
#           'H2Mu_WH_pos',
#           'H2Mu_WH_neg',
#           'H2Mu_ZH'
          ]

directory = '/home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/'
roch_filename = directory+'validate_UNBLINDED_dimu_mass_Roch_110_160_categories3_tree_categorization_final_36814_dyAMC_minpt10.root'
kamu_filename = directory+'validate_UNBLINDED_dimu_mass_KaMu_110_160_categories3_tree_categorization_final_36814_dyAMC_minpt10.root'

roch_list = []
kamu_list = []

for c in categories:
    for s in signals:
        roch_list.append((roch_filename, 'signal_histos/%s_%s' % (c, s))),
        kamu_list.append((kamu_filename, 'signal_histos/%s_%s' % (c, s))),

print "\n =========== GET ================== \n" 
roch_plots = tools.get(roch_list)
kamu_plots = tools.get(kamu_list)

#print "\n =========== FIT ================ \n" 
#
#ggf_ratios = [] 
#ggf_mean_roch = [] 
#ggf_mean_kamu = [] 
#ggf_sigma_roch = [] 
#ggf_sigma_kamu = [] 
#
#vbf_ratios = [] 
#vbf_mean_roch = [] 
#vbf_mean_kamu = [] 
#vbf_sigma_roch = [] 
#vbf_sigma_kamu = [] 
#
#for i in range(0,len(roch_plots)):
#    hroch = roch_plots[i]
#    hkamu = kamu_plots[i]
#    hroch.SetMarkerColor(1)
#    hroch.SetLineColor(1)
#    hroch.SetName(hroch.GetName()+"_Roch")
#    hroch.SetTitle(hroch.GetName()+"_vs_KaMu")
#    hkamu.SetMarkerColor(2)
#    hkamu.SetLineColor(2)
#    hkamu.SetName(hkamu.GetName()+"_KaMu")
#    hkamu.SetTitle(hkamu.GetName())
#    print "\n  /// %d /// Fits for for %s, %s" % (i, hroch.GetName(), hkamu.GetName())
#    c = TCanvas()
#    c.SetTitle(hroch.GetName()+"_vs_KaMu")
#    pad = TPad()
#    fit_roch = TF1("gaus_roch", "gaus", 123, 127);
#    fit_kamu = TF1("gaus_kamu", "gaus", 123, 127);
#    fit_roch.SetLineColor(1)
#    fit_kamu.SetLineColor(2)
#    hroch.Fit("gaus_roch", "", "", 125-1.2*hroch.GetRMS(), 125+1.2*hroch.GetRMS())
#    hkamu.Fit("gaus_kamu","" , "", 125-1.2*hroch.GetRMS(), 125+1.2*hroch.GetRMS())
#    mean_roch = fit_roch.GetParameter(1);
#    mean_kamu = fit_kamu.GetParameter(1);
#    sigma_roch = fit_roch.GetParameter(2);
#    sigma_kamu = fit_kamu.GetParameter(2);
#    sigma_ratio = sigma_roch/sigma_kamu;
#    t = TLatex(.6,.6,"#sigma_{Roch}/#sigma_{KaMu} = %4.3f" % sigma_ratio);
#    t.SetNDC(kTRUE);
#    hroch.Draw()
#    hkamu.Draw("SAME")
#    t.Draw("SAME")
#    c.SaveAs("fit_"+hroch.GetName()+"_vs_KaMu.png");
#    print "  +++ SIGMA kamu/roch: %f" % sigma_ratio
#    print ""
#    for c in categories:
#        if c+"_" in hroch.GetName() and "_gg_" in hroch.GetName():
#            ggf_ratios.append(sigma_ratio)
#            ggf_mean_roch.append(mean_roch)
#            ggf_mean_kamu.append(mean_kamu)
#            ggf_sigma_roch.append(sigma_roch)
#            ggf_sigma_kamu.append(sigma_kamu)
#        if c+"_" in hroch.GetName() and "_VBF_" in hroch.GetName():
#            vbf_ratios.append(sigma_ratio)
#            vbf_mean_roch.append(mean_roch)
#            vbf_mean_kamu.append(mean_kamu)
#            vbf_sigma_roch.append(sigma_roch)
#            vbf_sigma_kamu.append(sigma_kamu)
#        
#print "\n =========== LIST GGF RATIOS ================ \n" 
#for i in range(0,len(ggf_ratios)):
#    print "%s: %4.3f" % (categories[i], ggf_ratios[i]) 
#
#print "\n =========== LIST VBF RATIOS ================ \n" 
#for i in range(0,len(vbf_ratios)):
#    print "%s: %4.3f" % (categories[i], vbf_ratios[i]) 
#
#width = 10
#
#print "\n =========== LIST GGF INFO ================ \n" 
#
#print '{:<15}'.format('category') + '{:<15}'.format('mean_roch') + '{:<15}'.format('mean_kamu') + '{:<15}'.format('sigma_roch') + '{:<15}'.format('sigma_kamu') 
#for i in range(0,len(ggf_ratios)):
#    print ('{:<15}'.format(categories[i]) + '{:<15.2f}'.format(ggf_mean_roch[i]) + '{:<15.2f}'.format(ggf_mean_kamu[i]) + 
#          '{:<15.2f}'.format(ggf_sigma_roch[i]) + '{:<15.2f}'.format(ggf_sigma_kamu[i]) )
#
#print "\n =========== LIST VBF INFO ================ \n" 
#
#print '{:<15}'.format('category') + '{:<15}'.format('mean_roch') + '{:<15}'.format('mean_kamu') + '{:<15}'.format('sigma_roch') + '{:<15}'.format('sigma_kamu')
#for i in range(0,len(vbf_ratios)):
#    print ('{:<15}'.format(categories[i]) + '{:<15.2f}'.format(vbf_mean_roch[i]) + '{:<15.2f}'.format(vbf_mean_kamu[i]) + 
#          '{:<15.2f}'.format(vbf_sigma_roch[i]) + '{:<15.2f}'.format(vbf_sigma_kamu[i]) )

print "\n =========== REBIN and RATIO ================ \n" 

for i in range(0,len(roch_plots)):
    hroch = roch_plots[i]
    hkamu = kamu_plots[i]
    hroch.SetName(hroch.GetName()+"_Roch")
    hroch.SetTitle(hroch.GetTitle()+"_Roch")
    hkamu.SetName(hkamu.GetName()+"_KaMu")
    hkamu.SetTitle(hkamu.GetTitle()+"_KaMu")
    print "\n  /// %d /// Stack and Ratio for %s, %s" % (i, hroch.GetName(), hkamu.GetName())
    hlist = [hkamu, hroch]
    #newBins = tools.getRebinEdges(hroch, hkamu, max_err=0.1)
    newBins = np.linspace(110, 160, 50) 
    #print newBins
    rebin_hlist = tools.rebin(hlist, newBins)
    hroch = rebin_hlist[1]
    hkamu = rebin_hlist[0]
    tools.stackAndRatio(rebin_hlist, title=hroch.GetName()+"_vs_KaMu", ytitleratio='Roch/KaMu', log=False, yrange=(0, hroch.GetMaximum()*1.15), yrange_ratio=(0.85,1.15))
