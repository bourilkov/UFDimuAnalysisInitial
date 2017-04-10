##############################################
# overlay.py                                 #
##############################################

import tools as tools

#============================================
# code
#============================================

categories = [
              'T_0p1230_lt_bdt_score_0p400_lt_bdt_score_0p050_lt_bdt_score_n0p400',
              'T_0p1423_lt_bdt_score_0p400_gt_bdt_score_0p050_gt_dimu_max_abs_eta_0p900_gt_bdt_score_0p250_gt_dimu_max_abs_eta_1p900',
              'T_0p1580_lt_bdt_score_0p400_lt_bdt_score_0p050_gt_bdt_score_n0p400_gt_dimu_max_abs_eta_1p900',
              'T_0p2123_gt_bdt_score_0p400_lt_bdt_score_0p730_lt_bdt_score_0p650_gt_dimu_max_abs_eta_0p900_gt_dimu_max_abs_eta_1p900',
              'T_0p2948_lt_bdt_score_0p400_gt_bdt_score_0p050_lt_dimu_max_abs_eta_0p900_lt_bdt_score_0p250',
              'T_0p3127_lt_bdt_score_0p400_gt_bdt_score_0p050_lt_dimu_max_abs_eta_0p900_gt_bdt_score_0p250',
              'T_0p3528_lt_bdt_score_0p400_gt_bdt_score_0p050_gt_dimu_max_abs_eta_0p900_gt_bdt_score_0p250_lt_dimu_max_abs_eta_1p900',
              'T_0p3757_lt_bdt_score_0p400_gt_bdt_score_0p050_gt_dimu_max_abs_eta_0p900_lt_bdt_score_0p250',
              'T_0p3778_lt_bdt_score_0p400_lt_bdt_score_0p050_gt_bdt_score_n0p400_lt_dimu_max_abs_eta_1p900',
              'T_0p4006_gt_bdt_score_0p400_lt_bdt_score_0p730_gt_bdt_score_0p650',
              'T_0p4267_gt_bdt_score_0p400_lt_bdt_score_0p730_lt_bdt_score_0p650_gt_dimu_max_abs_eta_0p900_lt_dimu_max_abs_eta_1p900',
              'T_0p4621_gt_bdt_score_0p400_lt_bdt_score_0p730_lt_bdt_score_0p650_lt_dimu_max_abs_eta_0p900',
              'T_0p6890_gt_bdt_score_0p400_gt_bdt_score_0p730'
             ]

h = 'H2Mu_gg'

i=0
for c in categories:
    in_list = [
               ('../../bin/rootfiles/validate_M120_blinded_dimu_mass_PF_110_160_nolow_categories3_tree_categorization_final_36814_dyMG.root', 'net_histos/%s_%s_120' % (c, h)),
               ('../../bin/rootfiles/validate_M130_blinded_dimu_mass_PF_110_160_nolow_categories3_tree_categorization_final_36814_dyMG.root', 'net_histos/%s_%s_130' % (c, h)),
               ('../../bin/rootfiles/validate_blinded_dimu_mass_PF_110_160_nolow_categories3_tree_categorization_final_36814_dyMG.root'     , 'net_histos/%s_%s' % (c, h))
              ]

    plot_list = tools.get(in_list)
    plot_list[0].SetName('M120')
    plot_list[0].SetTitle('M120')
    plot_list[1].SetName('M130')
    plot_list[1].SetTitle('M130')
    plot_list[2].SetName('M125')
    plot_list[2].SetTitle('M125')
    
    tools.overlay(plot_list, title="c%d_120_125_130_%s" % (i, h), savename="c%d_120_125_130_%s" % (i,h), 
            xtitle="dimu_mass", ytitle="Events / 1 GeV", ldim=[0.5, 0.7, 0.88, 0.88])
    i+=1

