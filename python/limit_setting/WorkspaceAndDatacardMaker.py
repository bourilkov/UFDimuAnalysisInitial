##############################################
# WorkspaceAndDatacardMaker.py               #
##############################################
# make a workspace for a category to use     #
# with higgs combine                         #
#                                            #
##############################################

#============================================
# import
#============================================

import prettytable
import string
import re
import argparse
from ROOT import *

#============================================
# code
#============================================

class WorkspaceAndDatacardMaker:

    infilename = ''
    category = ''
    signal_hist = 0
    bkg_hist = 0
    net_hist= 0
    tfile = 0
    nuisance_params = []

    def __init__(self, infilename, category):
        self.infilename = infilename
        self.category = category
        self.tfile = TFile(infilename)
        self.setBackgroundHist()
        self.setSignalHist()
        self.setNetMCHist()
    
    def setBackgroundHist(self):
        self.bkg_hist = self.tfile.Get('net_histos/'+self.category+"_Net_Bkg")
    
    def setSignalHist(self):
        self.signal_hist = self.tfile.Get('net_histos/'+self.category+"_Net_Signal")
    
    def setNetMCHist(self):
        # add up the signal and background
        self.net_hist = self.signal_hist.Clone()
        self.net_hist.Add(self.bkg_hist)
        self.net_hist.SetName(self.category+'_Net_MC')
        self.net_hist.SetTitle(self.category+'_Net_MC')
    
    def makeWorkspace(self):
        # suppress all messages except those that matter
        RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
        print "="*78
        wspace = RooWorkspace(self.category)
    
        # Get the dimu mass histogram to use for limit setting
        nbins   = self.net_hist.GetNbinsX()
        massmin = self.net_hist.GetBinLowEdge(1)
        massmax = massmin + nbins*self.net_hist.GetBinWidth(1)
    
        #----------------------------------------
        # create di-muon mass variable
        # syntax:
        # <name>[initial-val, min-val, max-val]
        #----------------------------------------
        wspace.factory('x[120.0, %f, %f]' % (massmin, massmax))
        wspace.var('x').SetTitle('m_{#mu#mu}')
        wspace.var('x').setUnit('GeV')
    
        # define the set obs = (x)
        wspace.defineSet('obs', 'x')
    
        # make the set obs known to Python
        obs  = wspace.set('obs')
    
        # create binned dataset from histogram
        # needs to be named data_obs for higgs combine limit setting
        data = RooDataHist('data_obs', 'data_obs', RooArgList(obs), self.net_hist)
        # add data to workspace.
        # the RooCmdArg() is a workaround a PyROOT "feature"
        getattr(wspace,'import')(data, RooCmdArg())
    
        #----------------------------------------
        # create background model
        #----------------------------------------
        #wspace.factory('bmodel_norm[1.0, 0.1, 10]') # higgs combine fits the background and sets the norm on its own
        #wspace.var('bmodel_norm').setConstant()     # so we don't need this when we use the wkspace in conjuction with higgs combine
        wspace.factory('a1[ 5.0, -50, 50]')          # nuisance parameter1 for the background fit
        wspace.factory('a2[-1.0, -50, 50]')          # nuisance parameter2 for the background fit
        self.nuisance_params.append('a1')
        self.nuisance_params.append('a2')
    
        # define the background model 
        # should take this as an input instead of hard coding it
        # so that we can automatically get limits for different models
        wspace.factory('expr::f("-(a1*(x/100)+a2*(x/100)^2)",a1,a2,x)')
        # exp(c*x), c = 1 and x = f
        wspace.factory('Exponential::bmodel_'+self.category+'(f, 1)')
        bmodel  = wspace.pdf('bmodel')
    
        #----------------------------------------
        # create signal model
        #----------------------------------------
        #wspace.factory('smodel_norm[1.0, 0.001, 1000.0]') # higgs combine automatically runs through different normalizations
        #wspace.var('smodel_norm').setConstant()           # to get the limits so we don't need to figure it out here
        wspace.factory('mass[125, %f, %f]' % (massmin, massmax))
        wspace.factory('w[1.0, 0.1, 10]')
        # define the signal model 
        # should take this as an input instead of hard coding it
        # so that we can automatically get limits for different models
        wspace.factory('Gaussian::smodel_'+self.category+'(x, mass, w)')
        wspace.var('mass').setConstant()                   # just set the mass to 125 for now
        self.nuisance_params.append('w')
        self.nuisance_params.append('mass')
    
        #----------------------------------------
        # save data and signal & bg models for use
        # with higgs combine
        #----------------------------------------
        wspace.SaveAs(self.category+'_s.root')

    def makeShapeDatacard(self):
        # make whitespace the size of category
        width = max(len('smodel_'+self.category), len('process'))
        width+=4

        f = open(self.category+'_s.txt', 'w') 
        f.write('imax *\n')
        f.write('jmax *\n')
        f.write('kmax *\n')
        f.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        f.write('shapes * * '+self.category+'_s.root '+self.category+':$PROCESS\n')
        f.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        f.write('bin            '+self.category+'\n')
        f.write('observation    -1.0\n')
        f.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        f.write('bin'.ljust(width)+self.category.ljust(width)+self.category.ljust(width)+'\n')
        f.write('process'.ljust(width)+('smodel_'+self.category).ljust(width)+('bmodel_'+self.category).ljust(width)+'\n')
        f.write('process'.ljust(width)+'0'.ljust(width)+'1'.ljust(width)+'\n')
        f.write('rate'.ljust(width)+'{0:.2f}'.format(self.signal_hist.Integral()).ljust(width)+'{0:.2f}'.format(self.bkg_hist.Integral()).ljust(width)+'\n')
        f.write('----------------------------------------------------------------------------------------------------------------------------------\n')
      
        # get maximum length of the strings to figure out the width of the columns for the systematics section 
        pwidth = len('param')
        for n in self.nuisance_params:
            if len(n) > pwidth: pwdith = len(n)
        pwidth+=4
      
        # write the systematics section
        for n in self.nuisance_params:
            f.write(n.ljust(pwidth)+'param'.ljust(pwidth)+'0.0'.ljust(pwidth)+'0.1'.ljust(pwidth)+'\n')

    def makeTemplateRootFile(self):
        s = self.signal_hist.Clone()
        b = self.bkg_hist.Clone()
        data = self.net_hist.Clone()

        s.SetName('signal_'+self.category)
        b.SetName('bkg_'+self.category)
        data.SetName('data_obs')

        f = TFile(self.category+'_t.root', 'RECREATE');
        f.cd()

        s.Write()
        b.Write()
        data.Write()

        f.Write()
        f.Close()
        self.tfile.cd()

    def makeTemplateDatacard(self):
        # make whitespace the size of category
        width = max(len('signal_'+self.category), len('process'))
        width+=4

        f = open(self.category+'_t.txt', 'w') 
        f.write('imax *\n')
        f.write('jmax *\n')
        f.write('kmax *\n')
        f.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        f.write('shapes * * '+self.category+'_t.root $PROCESS\n')
        f.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        f.write('bin            '+self.category+'\n')
        f.write('observation    -1.0\n')
        f.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        f.write('bin'.ljust(width)+self.category.ljust(width)+self.category.ljust(width)+'\n')
        f.write('process'.ljust(width)+('signal_'+self.category).ljust(width)+('bkg_'+self.category).ljust(width)+'\n')
        f.write('process'.ljust(width)+'0'.ljust(width)+'1'.ljust(width)+'\n')
        f.write('rate'.ljust(width)+'{0:.2f}'.format(self.signal_hist.Integral()).ljust(width)+'{0:.2f}'.format(self.bkg_hist.Integral()).ljust(width)+'\n')
        f.write('----------------------------------------------------------------------------------------------------------------------------------\n')

print('program is running ...')
wdm = WorkspaceAndDatacardMaker('/home/acarnes/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/validate_dimu_mass_110_160_x69p2_8_0_X_MC_categories_27217.root', 'GGF_Tight') 
print wdm.infilename, wdm.category
wdm.makeWorkspace()
wdm.makeShapeDatacard()
wdm.makeTemplateRootFile()
wdm.makeTemplateDatacard()
