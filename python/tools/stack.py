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


filename = '/home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/validate_blinded_dimu_mass_PF_50_200_low_run1categories_36814.root'
category = 'c_ALL'
samples = ['VH', 'H2Mu_VBF', 'H2Mu_gg', 'Diboson_plus', 'TTbar_Plus_SingleTop', 'Drell_Yan', 'Net_Data']
getlist = []

for sample in samples:
    getlist.append((filename, 'net_histos/%s_%s' % (category, sample)))

print "\n =========== GET ================== \n" 
hlist = tools.get(getlist)
print "\n =========== STACKANDRATIO ======== \n" 
tools.stackAndRatio(hlist, title=category)
