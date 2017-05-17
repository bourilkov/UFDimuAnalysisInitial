##############################################
# h2mu_poly.py                               #
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
#categories = ['c8']

filedir = '/home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/'
#filename = 'validate_blinded_dimu_mass_Roch_110_160_categories3_tree_categorization_final_36814_dyAMC_minpt10.root';
#filename = 'validate_UNBLINDED_dimu_mass_Roch_110_160_categories3_tree_categorization_final_36814_dyAMC_minpt10.root';
#filename = 'validate_blinded_dimu_mass_Roch_110_160_categories3_tree_categorization_final_36814_dyAMC-J_minpt10_b1_sig-xlumi1.root'
#filename = 'validate_UNBLINDED_dimu_mass_Roch_100_150_categories3_tree_categorization_final_36814_dyAMC-J_minpt10_b-4_sig-xlumi1.root'
filename = 'validate_blinded_dimu_mass_Roch_100_150_categories3_tree_categorization_final_36814_dyAMC-J_minpt10_b4_sig-xlumi1.root'

order = 4 # order for bernstein poly
blinded = True

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
    
    #h2mupoly_model, h2mupoly_params   = pdfs.h2mupoly(x, order=8)
    #h2mupolyf_model, h2mupolyf_params   = pdfs.h2mupolyf(x, order=15)
    #h2mupolypow_model, h2mupolypow_params   = pdfs.h2mupolypow(x, order=6)
    bwzr_model, bwzr_params = pdfs.bwZreduxFixed(x)
    bwzg_model, bwzg_params = pdfs.bwZGamma(x)
    
    f = wdm.fit(histo, bwzr_model, x, blinded=blinded, save=True)
