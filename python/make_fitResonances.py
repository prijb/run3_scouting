#!/bin/env python3                                                                                                                                               
import ROOT, array, random, copy                                                                                                    
from ROOT import TCanvas, TFile, TH1, TH1F, TF1, gSystem                                                                                        
import ROOT, array, random, copy                                                                                        
from ROOT import RooCmdArg, RooArgSet, kFALSE, RooLinkedList, kBlue, kRed, kBlack, kOpenStar, kWhite, kGray                                    
from ROOT import gStyle, TStyle, TGraph, TGraphErrors, TMath, TMultiGraph, TLine, gPad, TGaxis, TLegend, TText, TLatex, TColor, TPaveText      
from ROOT import TAttFill, TLegend, TRatioPlot, TPad, THStack, TFileCollection                                                      
from ROOT import kBlue, kRed, kBlack, kWhite, kAzure, kOrange, kPink, kGreen, kYellow, kCyan, kMagenta                                         
from ROOT import RooRealVar, RooDataHist, RooCBShape, RooDoubleCBFast, RooGaussian, RooExponential, RooAddPdf, RooFit, RooArgList, RooPolynomial
from ROOT import RooCmdArg, RooArgSet, RooLinkedList, kFALSE, kBlue, kRed, kBlack, kGray
from ROOT import gStyle, TCanvas, TLatex, TLegend, TFile
import math                                                                                                                             
import os                                                                                                            
import argparse                                                                                                                         
import sys
from datetime import date
import numpy as np
import mplhep as hep
import matplotlib.pyplot as plt
                                                                                                       

ROOT.gROOT.SetBatch()                                                                                          
ROOT.gStyle.SetOptStat(0)                                                                                    
ROOT.gStyle.SetOptTitle(0)

def getValues(histo):
    values = []
    bins = []
    for n in range(1, histo.GetNbinsX()+1):
        values.append(histo.GetBinContent(n))
        bins.append(histo.GetBinLowEdge(n))
    bins.append(histo.GetBinLowEdge(n) + histo.GetBinWidth(n))
    return np.array(values), np.array(bins)

def fitDataset(name, dataset, mass, fit_range=(2.6, 3.6), outDir="output"):
    # Signal parameters initialization
    stddev = 0.02*mass
    minstddev = 0.001
    maxstddev = 5.0*stddev
    #
    mmean = mass
    #
    alphaR = 1.0;
    minalphaR = 0.1;
    maxalphaR = 10.0;
    #
    alphaL = 4.0;
    minalphaL = 0.1;
    maxalphaL = 10.0;
    #
    nR = 6.0;
    minnR = 1.0;
    maxnR = 15.0;
    #
    nL = 2.0;
    minnL = 1.0;
    maxnL = 5.0;
    #
    #
    ## Frame:
    mfit = RooRealVar("mfit", "mfit", fit_range[0], fit_range[1])
    data = dataset.reduce(mfit, "%f < mfit && mfit < %f"%(fit_range[0], fit_range[1]))
    data.Print()
    x = data.get().find("mfit")
    x.setRange("fitRange", fit_range[0], fit_range[1])
    nBins = 100
    histogram = data.createHistogram("mfit", 100)  # 50 bins
    data_hist = RooDataHist("data_hist", "Binned data", RooArgSet(x), histogram)
    #
    ## Signal pdfs
    mean_ = RooRealVar("mean", "mean of crystal ball", 3.1, 2.85, 3.25)
    sigma_ = RooRealVar("sigma", "width of crystal ball", stddev, minstddev, maxstddev)
    gauss = RooGaussian("gauss", "gauss", x, mean_, sigma_)
    #
    alphaL_ = RooRealVar("alphaL", "alphaL", alphaL, minalphaL, maxalphaL)
    alphaR_ = RooRealVar("alphaR", "alphaR", alphaR, minalphaR, maxalphaR)
    nL_ = RooRealVar("nL", "nL", nL, minnL, maxnL)
    nR_ = RooRealVar("nR", "nR", nR, minnR, maxnR)
    dcb = RooDoubleCBFast("dcb", "dcb", x, mean_, sigma_, alphaL_, nL_, alphaR_, nR_)
    #
    #
    ## Background pdfs
    # Background exponential
    expo_slope  = RooRealVar("expo_slope","expo_slope",-0.1,-20.0,20.0);
    background = RooExponential("background_exponential", "background_exponential",x,expo_slope);
    #b0 = RooRealVar("b0", "constant term", 1.0, 0.0, 10.0)
    #b1 = RooRealVar("b1", "linear term", 0.0, -1.0, 1.0)
    #background = RooPolynomial("background", "background", x, RooArgList(b0, b1))

    # Combine the models
    frac_gauss = RooRealVar("frac_gauss", "fraction of Gaussian", 0.5, 0.0, 1.0)
    frac_dcb = RooRealVar("frac_dcb", "fraction of crystal ball", 0.3, 0.0, 1.0)
    model = RooAddPdf("model", "DCB + Gaussian + Background", RooArgList(dcb, gauss, background), RooArgList(frac_dcb, frac_gauss))

    # Fit
    nMaxFitAttempts = 1
    nFits = 0
    while (nFits < nMaxFitAttempts):
        fit_result = model.fitTo(data_hist, RooFit.Range(fit_range[0], fit_range[1]), RooFit.Save())
        nFits+=1
        if fit_result.status()==0:
            break

    # Create a frame to plot the fit result
    frame = x.frame(RooFit.Title("Mass"))
    data_hist.plotOn(frame)
    model.plotOn(frame)
    #
    frame.SetXTitle("Mass [GeV]")
    frame.GetXaxis().SetTitleOffset(1.3)
    frame.SetYTitle("Events / 0.01 GeV")
    #
    #c1 = TCanvas("c1", "")
    #frame.Draw()
    #c1.SaveAs("pruebaFit.png")

    # Get Chi2 and normalize by ndof and number of bins
    chi2 = frame.chiSquare(fit_result.floatParsFinal().getSize())  # chi^2
    ndof = fit_result.floatParsFinal().getSize()  # Number of degrees of freedom
    chi2_over_ndof = chi2 / (nBins - ndof)
    print("----------------------------------------")
    print("nBins = ", nBins)
    print("chi2_over_ndof = ", chi2_over_ndof)
    print("----------------------------------------")    

    ## Plotting
    histo, bins = getValues(histogram)
    y_err = histo**0.5
    xval = np.linspace(fit_range[0], fit_range[1], 200)
    vbkg = []  
    vdcb = []  
    vgaus = []  
    vmodel = []  
    model_integral_value = model.createIntegral(RooArgSet(x)).getVal()
    bkg_integral_value = background.createIntegral(RooArgSet(x)).getVal()
    gauss_integral_value = gauss.createIntegral(RooArgSet(x)).getVal()
    dcb_integral_value = dcb.createIntegral(RooArgSet(x)).getVal()
    binw = histogram.GetBinLowEdge(2)-histogram.GetBinLowEdge(1)
    #    
    for val in xval:
        x.setVal(val)
        frac_bkg = 1.0 - frac_dcb.getVal()-frac_gauss.getVal()
        vbkg.append(background.getVal()*data.sumEntries()*binw*(frac_bkg)/bkg_integral_value)
        vdcb.append(dcb.getVal()*data.sumEntries()*binw*frac_dcb.getVal()/dcb_integral_value)
        vgaus.append(gauss.getVal()*data.sumEntries()*binw*frac_gauss.getVal()/gauss_integral_value)
        vmodel.append(vbkg[-1]+vdcb[-1]+vgaus[-1])
    #
    plt.style.use(hep.style.CMS)
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_ylabel(r'Events / 0.01 GeV', fontsize=24)
    ax.set_xlabel(r'Dimuon mass (GeV)', fontsize=24)
    ax.set_xlim(2.6,3.6)
    hep.cms.label("Preliminary", data=True, lumi=1, year=2024, com='13.6')
    ax.set_ylim(0.0,1.2*max(histo))
    ax.errorbar(bins[:-1], histo, yerr=y_err, fmt='o', capsize=5, label='Data', color='k', markersize=8)
    ax.plot(xval, vbkg, label='Background', color='firebrick', lw = 3)
    ax.plot(xval, vgaus, label='Gaussian', color='deepskyblue', lw = 3)
    ax.plot(xval, vdcb, label='Double Crystal Ball', color='violet', lw = 3)
    ax.plot(xval, vmodel, label='Signal + Background', color='slateblue', lw = 3)
    if 'L1' in name:
        ax.text(0.45, 0.95, 'L1' + name.split('L1')[1], fontsize=15, color='black', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
    ax.text(0.65, 0.8, 'Fit parameters:', fontsize=18, color='black', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
    ax.text(0.65, 0.76, r'Mean = %.3f GeV'%(mean_.getVal()), fontsize=18, color='black', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
    ax.text(0.65, 0.72, r'$\sigma =$ %.3f GeV'%(sigma_.getVal()), fontsize=18, color='black', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
    ax.legend(loc='upper left', fontsize = 18, frameon = True, ncol=1)
    fig.savefig('%s/%s_fit.png'%(outDir,name), dpi=140)  

    return fit_result, model

if __name__=="__main__":

    dNames = []
    dNames.append("d_Dimuon_lxy0p0to0p2_inclusive")
    dNames.append("d_Dimuon_lxy0p2to1p0_inclusive")
    dNames.append("d_Dimuon_lxy1p0to2p4_inclusive")
    dNames.append("d_Dimuon_lxy2p4to3p1_inclusive")
    dNames.append("d_Dimuon_lxy3p1to7p0_inclusive")
    dNames.append("d_Dimuon_lxy7p0to11p0_inclusive")
    dNames.append("d_Dimuon_lxy11p0to16p0_inclusive")
    dNames.append("d_Dimuon_lxy16p0to70p0_inclusive")

    # Input
    inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Sep-25-2024_RooDatasets_unblind"
    files = os.listdir(inDir)
    #
    eras = []
    eras.append("DataC")
    eras.append("DataD")
    eras.append("DataE")
    eras.append("DataF")
    eras.append("DataG")
    #
    today = date.today().strftime("%b-%d-%Y")
    outDir = "plotsMasking_%s"%(today)
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    #
    # Loop:
    for d_,d in enumerate(dNames):
        for e,era in enumerate(eras):
            if e==0:
                fin = ROOT.TFile.Open("%s/%s"%(inDir,files[0]))
                dataset = fin.Get(d).Clone()
            for f,file in enumerate(files):
                if era in file:
                    if e==0 and f==0:
                        continue
                    print("Reading %s"%(file))
                    fin = ROOT.TFile.Open("%s/%s"%(inDir,file))
                    tds = fin.Get(d).Clone()
                    dataset.append(tds)
        print(dataset.numEntries())
        # Fitting:
        fitDataset("JPsi_%s"%(d), dataset, mass = 3.1, fit_range=(2.6, 3.6), outDir=outDir)




