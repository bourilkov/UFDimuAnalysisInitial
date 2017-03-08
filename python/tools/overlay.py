##############################################
# overlay.py                                 #
##############################################

import tools as tools

#============================================
# code
#============================================

runs = ['Net_Data', 'RunB', 'RunC', 'RunD', 'RunE', 'RunF', 'RunG', 'RunH']


for run in runs:
    plot_list = [('/home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/zcalibration_eta_plus_muID0_data_MC.root', 'plots/resolution_%s_mass_PF_eta_plus' % run),
                 ('/home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/zcalibration_eta_plus_muID1_data_MC.root', 'plots/resolution_%s_mass_PF_eta_plus' % run),
                 ('/home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/zcalibration_eta_plus_muID2_data_MC.root', 'plots/resolution_%s_mass_PF_eta_plus' % run)
                ]   
    
    tools.overlay(plot_list, title="%s_Z_Resolution_and_Mu_ID" % run, savename="%s_resolution_muID_eta_plus" % run, 
            xtitle="#mu^{+} Eta", ytitle="Z Resolution", yrange=(0,5))

