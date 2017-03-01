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
def bwGamma(x):
    expParam = RooRealVar("bwg_expParam","expParam",-1e-03,-1e-01,1e-01)
    bwmodel = RooGenericPdf("bwg_model","exp(@0*@1)*pow(@0,-2)",RooArgList(x,expParam))

    return bwmodel, [expParam]
#--------------------------------------------------------
# breit weigner mixture scaled by falling exp (run1 bg)
#--------------------------------------------------------
def bwZ(x):
    bwWidth =  RooRealVar("bwz_Width","widthZ",2.5,0,30)
    bwmZ =     RooRealVar("bwz_mZ","mZ",91.2,90,92)
    expParam = RooRealVar("bwz_expParam","expParam",-1e-03,-1e-02,1e-02)
    
    bwWidth.setConstant(True);
    bwmZ.setConstant(True);
    
    bwmodel  = RooGenericPdf("bwz_model","exp(@0*@3)*(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",RooArgList(x,bwmZ,bwWidth,expParam))
    return bwmodel, [bwWidth, bwmZ, expParam]

#--------------------------------------------------------
# breit weigner mixture scaled by falling exp (run1 bg)
#--------------------------------------------------------
def bwZGamma(x, mix_min=0.001):
    bwWidth =  RooRealVar("bwzg_Width","widthZ",2.5,0,30)
    bwmZ =     RooRealVar("bwzg_mZ","mZ",91.2,90,92)
    expParam = RooRealVar("bwzg_expParam","expParam",-0.0053,-0.0073,-0.0033)
    mixParam = RooRealVar("bwzg_mixParam","mix",0.379,0.359,0.399)
    
    bwWidth.setConstant(True);
    bwmZ.setConstant(True);
    
    phoExpMmumu = RooGenericPdf("phoExpMmumu","exp(@0*@1)*pow(@0,-2)",RooArgList(x,expParam))
    bwExpMmumu  = RooGenericPdf("bwExpMmumu","exp(@0*@3)*(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",RooArgList(x,bwmZ,bwWidth,expParam))
    bwmodel     = RooAddPdf("bwzg_model","bw_model", RooArgList(bwExpMmumu,phoExpMmumu),RooArgList(mixParam))

    return bwmodel, [bwWidth, bwmZ, expParam, mixParam, phoExpMmumu, bwExpMmumu]
    
#--------------------------------------------------------
# breit weigner mixture scaled by falling exp (run1 bg)
#--------------------------------------------------------
def bwZGammaPlusLinear(x):
    bwWidth =  RooRealVar("bwzgl_widthZ","widthZ",2.5,0,30)
    bwmZ =     RooRealVar("bwzgl_mZ","mZ",91.2,85,95)
    expParam = RooRealVar("bwzgl_expParam","expParam",-1e-03,-1e-01,1e-01)

    bwWidth.setConstant(True);
    bwmZ.setConstant(True);

    slopeParam = RooRealVar("bwl_slope", "slope", -0.2, -50, 0)          
    offsetParam = RooRealVar("bwl_offset", "offset", 39, 0, 1000)            
    
    mix1 = RooRealVar("bwzgl_mix1","mix1",0.95,0,1)
    mix2 = RooRealVar("bwzgl_mix2","mix2",0.05,0,1)

    linMmumu = RooGenericPdf("bwzgl_linMmumu", "@1*@0+@2", RooArgList(x, slopeParam, offsetParam))
    phoExpMmumu = RooGenericPdf("bwzgl_phoExpMmumu","exp(@0*@1)*pow(@0,-2)",RooArgList(x,expParam))
    bwExpMmumu  = RooGenericPdf("bwzgl_bwExpMmumu","exp(@0*@3)*(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",RooArgList(x,bwmZ,bwWidth,expParam))
    model     = RooAddPdf("bwzgl_model","bwl_model", RooArgList(bwExpMmumu,phoExpMmumu,linMmumu),RooArgList(mix1, mix2))

    return model, [bwWidth, bwmZ, expParam, mix1, mix2, slopeParam, offsetParam, phoExpMmumu, bwExpMmumu, linMmumu]
    
#--------------------------------------------------------
# breit weigner scaled by falling exp + line
#--------------------------------------------------------
def bwZPlusLinear(x):
    bwWidth =  RooRealVar("bwzl_widthZ","widthZ",2.5,0,30)
    bwmZ =     RooRealVar("bwzl_mZ","mZ",91.2,85,95)
    expParam = RooRealVar("bwzl_expParam","expParam",-1e-03,-1e-01,1e-01)

    bwWidth.setConstant(True);
    bwmZ.setConstant(True);

    slopeParam = RooRealVar("bwzl_slope", "slope", -0.2, -50, 0)          
    offsetParam = RooRealVar("bwzl_offset", "offset", 39, 0, 1000)            
    
    mix1 = RooRealVar("bwzl_mix1","mix1",0.95,0,1)

    linMmumu = RooGenericPdf("bwzl_linMmumu", "@1*@0+@2", RooArgList(x, slopeParam, offsetParam))
    bwExpMmumu  = RooGenericPdf("bwzl_bwExpMmumu","exp(@0*@3)*(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",RooArgList(x,bwmZ,bwWidth,expParam))
    model     = RooAddPdf("bwzl_model","bwzl_model", RooArgList(bwExpMmumu,linMmumu),RooArgList(mix1))

    return model, [bwWidth, bwmZ, expParam, mix1, slopeParam, offsetParam, bwExpMmumu, linMmumu]

#----------------------------------------
# fewz falling exponential
#----------------------------------------
def bwZredux(x):
    a1 = RooRealVar("bwz_redux_a1", "a1", 1.39, 0.7, 2.1)
    a2 = RooRealVar("bwz_redux_a2", "a2", 0.46, 0.30, 0.62)
    a3 = RooRealVar("bwz_redux_a3", "a3", -0.26, -0.40, -0.12)

    #a1.setConstant()
    #a2.setConstant()
    #a3.setConstant()
    
    f = RooFormulaVar("bwz_redux_f", "(@1*(@0/100)+@2*(@0/100)^2)", RooArgList(x, a2, a3))
    expmodel = RooGenericPdf("bwz_redux_model", "exp(@2)*(2.5)/(pow(@0-91.2,@1)+0.25*pow(2.5,@1))", RooArgList(x, a1, f))
    return expmodel, [a1, a2, a3, f]

#----------------------------------------
# hgg falling exponential
#----------------------------------------
def higgsGammaGamma(x):
    a1 = RooRealVar("hgg_a1", "a1", 5.0, -50, 50)          # nuisance parameter1 for the background fit
    a2 = RooRealVar("hgg_a2", "a2", -1.0, -50, 50)         # nuisance parameter2 for the background fit
    one = RooRealVar("hgg_one", "one", 1.0, -10, 10) 
    one.setConstant()
    
    f = RooFormulaVar("hgg_f", "-(@1*(@0/100)+@2*(@0/100)^2)", RooArgList(x, a1, a2))
    expmodel = RooExponential('hggexp_model', 'hggexp_model', f, one) # exp(1*f(x))

    return expmodel, [a1, a2, one, f]

#----------------------------------------
# chebychev
#----------------------------------------
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
