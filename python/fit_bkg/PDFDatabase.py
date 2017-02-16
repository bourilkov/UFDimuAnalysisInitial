##############################################
# PDFDatabase.py                             #
##############################################

#============================================
# import
#============================================

import prettytable
import string
import re
import argparse
from ROOT import *

import sys
sys.argv.append( '-b-' )

#----------------------------------------
# linear
#----------------------------------------
def linear(x):
    m = RooRealVar("m", "m", -0.2, -50, 50)          # nuisance parameter1 for the background fit
    b = RooRealVar("b", "b", 39, 0, 1000)            # nuisance parameter2 for the background fit
    
    linear_model = RooGenericPdf("linear_model", "@1*@0+@2", RooArgList(x, m, b))
    return linear_model, [m, b]

#--------------------------------------------------------
# breit weigner mixture scaled by falling exp (run1 bg)
#--------------------------------------------------------
def bwZGamma(x):
    bwWidth =  RooRealVar("bwWidth","widthZ",2.5,0,30)
    bwmZ =     RooRealVar("bwmZ","mZ",91.2,90,92)
    expParam = RooRealVar("expParam","expParam",-1e-03,-1e-01,1e-01)
    mixParam = RooRealVar("mixParam","mix",0.5,0,1)
    
    bwWidth.setConstant(True);
    bwmZ.setConstant(True);
    
    phoExpMmumu = RooGenericPdf("phoExpMmumu","exp(@0*@1)*pow(@0,-2)",RooArgList(x,expParam))
    bwExpMmumu  = RooGenericPdf("bwExpMmumu","exp(@0*@3)*(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",RooArgList(x,bwmZ,bwWidth,expParam))
    bwmodel     = RooAddPdf("bw_model","bw_model", RooArgList(bwExpMmumu,phoExpMmumu),RooArgList(mixParam))

    return bwmodel, [bwWidth, bwmZ, expParam, mixParam, phoExpMmumu, bwExpMmumu]
    
#--------------------------------------------------------
# breit weigner mixture scaled by falling exp (run1 bg)
#--------------------------------------------------------
def bwZGammaPlusLinear(x):
    bwWidth =  RooRealVar("bwl_widthZ","widthZ",2.5,0,30)
    bwmZ =     RooRealVar("bwl_mZ","mZ",91.2,85,95)
    expParam = RooRealVar("bwl_expParam","expParam",-1e-03,-1e-01,1e-01)

    bwWidth.setConstant(True);
    bwmZ.setConstant(True);

    slopeParam = RooRealVar("bwl_slope", "slope", -0.2, -50, 0)          
    offsetParam = RooRealVar("bwl_offset", "offset", 39, 0, 1000)            
    
    mix1 = RooRealVar("bwl_mix1","mix1",0.95,0.5,1)
    mix2 = RooRealVar("bwl_mix2","mix2",0.05,0,0.5)

    linMmumu = RooGenericPdf("bwl_linMmumu", "@1*@0+@2", RooArgList(x, slopeParam, offsetParam))
    phoExpMmumu = RooGenericPdf("bwl_phoExpMmumu","exp(@0*@1)*pow(@0,-2)",RooArgList(x,expParam))
    bwExpMmumu  = RooGenericPdf("bwl_bwExpMmumu","exp(@0*@3)*(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",RooArgList(x,bwmZ,bwWidth,expParam))
    model     = RooAddPdf("bwl_model","bwl_model", RooArgList(bwExpMmumu,phoExpMmumu,linMmumu),RooArgList(mix1, mix2))

    return model, [bwWidth, bwmZ, expParam, mix1, mix2, slopeParam, offsetParam, phoExpMmumu, bwExpMmumu, linMmumu]
    
#----------------------------------------
# falling exponential (hgammgamma bg)
#----------------------------------------
def higgsGammaGamma(x):
    a1 = RooRealVar("a1", "a1", 5.0, -50, 50)          # nuisance parameter1 for the background fit
    a2 = RooRealVar("a2", "a2", -1.0, -50, 50)         # nuisance parameter2 for the background fit
    one = RooRealVar("one", "one", 1.0, -10, 10) 
    one.setConstant()
    
    f = RooFormulaVar("f", "-(@1*(@0/100)+@2*(@0/100)^2)", RooArgList(x, a1, a2))
    expmodel = RooExponential('hggexp_model', 'hggexp_model', f, one) # exp(1*f(x))

    return expmodel, [a1, a2, one, f]

def chebychev(x, order=7): 
    #c0 = RooRealVar("c0","c0", 1.0,-1.0,1.0)
    #c1 = RooRealVar("c1","c1", 1.0,-1.0,1.0)
    #c2 = RooRealVar("c2","c2", 1.0,-1.0,1.0)

    args = RooArgList()
    params = []
    for i in range(0,order):
        c = RooRealVar("c"+str(i),"c"+str(i), 1.0/2**i,-1.0,1.0)
        args.add(c)
        params.append(c)

    chebychev = RooChebychev("chebychev"+str(order),"chebychev"+str(order), x,args) 
    return chebychev, params
