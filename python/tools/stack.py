##############################################
# stack.py                                   #
##############################################

import tools as tools
import sys
sys.argv.append( '-b-' )
from ROOT import *

#============================================
# code
#============================================


dirname   = '/home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/'
filenames = ['validate_blinded_dimu_mass_PF_50_200_nolow_run1categories_36814.root', 
             'validate_blinded_dimu_mass_PF_50_200_low_run1categories_36814.root'
            ]

filename = dirname+filenames[0]

categories = ['c_ALL',
              'c_2_Jet',
              'c_2_Jet_VBF_Tight',
              'c_2_Jet_GGF_Tight',
              'c_2_Jet_VBF_Loose',
              'c_01_Jet',
              'c_01_Jet_Tight',
              'c_01_Jet_Tight_BB',
              'c_01_Jet_Tight_BO',
              'c_01_Jet_Tight_BE',
              'c_01_Jet_Tight_OO',
              'c_01_Jet_Tight_OE',
              'c_01_Jet_Tight_EE',
              'c_01_Jet_Loose',
              'c_01_Jet_Loose_BB',
              'c_01_Jet_Loose_BO',
              'c_01_Jet_Loose_BE',
              'c_01_Jet_Loose_OO',
              'c_01_Jet_Loose_OE',
              'c_01_Jet_Loose_EE'
              ]
category = 'c_01_Jet_Tight_BB'

samples = ['VH', 'H2Mu_VBF', 'H2Mu_gg', 'Diboson_plus', 'TTbar_Plus_SingleTop', 'Drell_Yan', 'Net_Data']
getlist = []

for sample in samples:
    getlist.append((filename, 'net_histos/%s_%s' % (category, sample)))

print "\n =========== GET ================== \n" 
hlist = tools.get(getlist)
print "\n =========== REBIN ================ \n" 
hmc = tools.add(hlist[0:-1])
hdata = hlist[-1]
newBins = tools.getRebinEdges(hdata, hmc, max_err=0.1)
print newBins
rebin_hlist = tools.rebin(hlist, newBins) #rebinned to var bin width
srebin_hlist = tools.scaleByBinWidth(rebin_hlist, newBins) # rebinned to var bin width and wide bins are scaled down
print "\n =========== STACK AND RATIO ====== \n" 
tools.stackAndRatio(srebin_hlist, title=category)
