##############################################
# systematics.py                             #
##############################################

import tools as tools

#============================================
# code
#============================================

categories = [
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
           'H2Mu_WH_pos',
           'H2Mu_WH_neg',
           'H2Mu_ZH'
          ]

directory = '/home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/'
filename = directory+'validate_UNBLINDED_dimu_mass_Roch_110_160_categories3_tree_categorization_final_36814_dyAMC_minpt10.root'

in_list = []
in_list_JES_up = []
in_list_JES_down = []
in_list_PU_up = []
in_list_PU_down = []

for c in categories:
    for s in signals:
        in_list.append((filename, 'signal_histos/%s_%s' % (c, s))),
        in_list_JES_up.append((filename, 'signal_histos/%s_%s_JES_up' % (c, s))),
        in_list_JES_down.append((filename, 'signal_histos/%s_%s_JES_down' % (c, s))),
        in_list_PU_up.append((filename, 'signal_histos/%s_%s_PU_up' % (c, s))),
        in_list_PU_down.append((filename, 'signal_histos/%s_%s_PU_down' % (c, s))),

std_plots = tools.get(in_list)
JES_up_plots = tools.get(in_list_JES_up)
JES_down_plots = tools.get(in_list_JES_down)
PU_up_plots = tools.get(in_list_PU_up)
PU_down_plots = tools.get(in_list_PU_down)

print ""
print "=========================================="
print "JES Down/Up"
print "=========================================="

for i in range(0,len(std_plots)):
    if i%5 == 0: print ""
    JES_up_ratio = JES_up_plots[i].Integral()/std_plots[i].Integral()
    JES_down_ratio = JES_down_plots[i].Integral()/std_plots[i].Integral()
    print "JES %s: %5.3f/%5.3f" % (std_plots[i].GetName(), JES_down_ratio, JES_up_ratio)

print ""
print "=========================================="
print "PU Down/Up"
print "=========================================="

for i in range(0,len(std_plots)):
    if i%5 == 0: print ""
    PU_up_ratio = PU_up_plots[i].Integral()/std_plots[i].Integral()
    PU_down_ratio = PU_down_plots[i].Integral()/std_plots[i].Integral()
    print "PU %s: %5.3f/%5.3f" % (std_plots[i].GetName(), PU_down_ratio, PU_up_ratio)
