import ROOT
import numpy
import copy
import os,sys
from datetime import date
#import plotUtils

ROOT.gROOT.SetBatch(1)

user = os.environ.get("USER")
today= date.today().strftime("%b-%d-%Y")

thisDir = os.environ.get("PWD")
inDir  = "%s/datacards_all_May-06-2024_2023/"%thisDir

outDir = ("%s/combineFits_"%(thisDir))+today
if not os.path.exists(outDir):
    os.makedirs(outDir)

mass = "2.5"
ctau = "1"

dNames = []
"""
dNames.append(["ch1", "d_FourMu_sep"])
dNames.append(["ch2", "d_FourMu_osv"])
dNames.append(["ch3", "d_Dimuon_lxy0p0to0p2_iso0_ptlow"])
dNames.append(["ch4", "d_Dimuon_lxy0p0to0p2_iso0_pthigh"])
dNames.append(["ch5", "d_Dimuon_lxy0p0to0p2_iso1_ptlow"])
dNames.append(["ch6", "d_Dimuon_lxy0p0to0p2_iso1_pthigh"])
dNames.append(["ch7", "d_Dimuon_lxy0p2to1p0_iso0_ptlow"])
dNames.append(["ch8", "d_Dimuon_lxy0p2to1p0_iso0_pthigh"])
dNames.append(["ch9", "d_Dimuon_lxy0p2to1p0_iso1_ptlow"])
dNames.append(["ch10", "d_Dimuon_lxy0p2to1p0_iso1_pthigh"])
dNames.append(["ch11", "d_Dimuon_lxy1p0to2p4_iso0_ptlow"])
dNames.append(["ch12", "d_Dimuon_lxy1p0to2p4_iso0_pthigh"])
dNames.append(["ch13", "d_Dimuon_lxy1p0to2p4_iso1_ptlow"])
dNames.append(["ch14", "d_Dimuon_lxy1p0to2p4_iso1_pthigh"])
dNames.append(["ch15", "d_Dimuon_lxy2p4to3p1_iso0_ptlow"])
dNames.append(["ch16", "d_Dimuon_lxy2p4to3p1_iso0_pthigh"])
dNames.append(["ch17", "d_Dimuon_lxy2p4to3p1_iso1_ptlow"])
dNames.append(["ch18", "d_Dimuon_lxy2p4to3p1_iso1_pthigh"])
dNames.append(["ch19", "d_Dimuon_lxy3p1to7p0_iso0_ptlow"])
dNames.append(["ch20", "d_Dimuon_lxy3p1to7p0_iso0_pthigh"])
dNames.append(["ch21", "d_Dimuon_lxy3p1to7p0_iso1_ptlow"])
dNames.append(["ch22", "d_Dimuon_lxy3p1to7p0_iso1_pthigh"])
"""
dNames.append(["ch23", "d_Dimuon_lxy7p0to11p0_iso0_ptlow"])
dNames.append(["ch24", "d_Dimuon_lxy7p0to11p0_iso0_pthigh"])
dNames.append(["ch25", "d_Dimuon_lxy7p0to11p0_iso1_ptlow"])
dNames.append(["ch26", "d_Dimuon_lxy7p0to11p0_iso1_pthigh"])
dNames.append(["ch27", "d_Dimuon_lxy11p0to16p0_iso0_ptlow"])
dNames.append(["ch28", "d_Dimuon_lxy11p0to16p0_iso0_pthigh"])
dNames.append(["ch29", "d_Dimuon_lxy11p0to16p0_iso1_ptlow"])
dNames.append(["ch30", "d_Dimuon_lxy11p0to16p0_iso1_pthigh"])
dNames.append(["ch31", "d_Dimuon_lxy16p0to70p0_iso0_ptlow"])
dNames.append(["ch32", "d_Dimuon_lxy16p0to70p0_iso0_pthigh"])
dNames.append(["ch33", "d_Dimuon_lxy16p0to70p0_iso1_ptlow"])
dNames.append(["ch34", "d_Dimuon_lxy16p0to70p0_iso1_pthigh"])
dNames.append(["ch35", "d_Dimuon_lxy0p0to0p2_non-pointing"])
dNames.append(["ch36", "d_Dimuon_lxy0p2to1p0_non-pointing"])
dNames.append(["ch37", "d_Dimuon_lxy1p0to2p4_non-pointing"])
dNames.append(["ch38", "d_Dimuon_lxy2p4to3p1_non-pointing"])
dNames.append(["ch39", "d_Dimuon_lxy3p1to7p0_non-pointing"])
dNames.append(["ch40", "d_Dimuon_lxy7p0to11p0_non-pointing"])
dNames.append(["ch41", "d_Dimuon_lxy11p0to16p0_non-pointing"])
dNames.append(["ch42", "d_Dimuon_lxy16p0to70p0_non-pointing"])

def drawLabels(year="all",lumi=59.83+41.48+19.5+16.8,plotData=False):
    # Labels
    latex = ROOT.TLatex()
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(0.04)
    latex.SetNDC(True)

    latexCMS = ROOT.TLatex()
    latexCMS.SetTextFont(61)
    latexCMS.SetTextSize(0.055)
    latexCMS.SetNDC(True)

    latexCMSExtra = ROOT.TLatex()
    latexCMSExtra.SetTextFont(52)
    latexCMSExtra.SetTextSize(0.04)
    latexCMSExtra.SetNDC(True)

    legoffset = 0.0
    latexSel = ROOT. TLatex()
    latexSel.SetTextAlign(31)
    latexSel.SetTextFont(42)
    latexSel.SetTextSize(0.04)
    latexSel.SetNDC(True)

    yearenergy=""
    if year!="all" or lumi<100.0:
        if year!="all":
            yearenergy="%.1f fb^{-1} (%s, 13 TeV)"%(lumi,year)
        else:
            yearenergy="%.1f fb^{-1} (2022, 13.6 TeV)"%(lumi)
    else:
        yearenergy="%.0f fb^{-1} (13 TeV)"%(lumi)
    if plotData:
        cmsExtra="Preliminary"
    else:
        cmsExtra="Simulation"

    # Draw CMS headers
    expoffset=0
    if doRatio:
        latex.DrawLatex(0.95, 0.93+expoffset, yearenergy);
        latexCMS.DrawLatex(0.15,0.93+expoffset,"CMS");
        latexCMSExtra.DrawLatex(0.24,0.93+expoffset, cmsExtra);
    else:
        latex.DrawLatex(0.95, 0.93+expoffset, yearenergy);
        latexCMS.DrawLatex(0.11,0.93+expoffset,"CMS");
        latexCMSExtra.DrawLatex(0.21,0.93+expoffset, cmsExtra);


def getLegend(ch,gd,p,hp,hs,smodel,smass,sigsf=-1.0,plotSignal=True,plotData=True,plotOnlySignal=False,hmc=None,hgauss=None,hdcb=None):
    legend = ROOT.TLegend(0.5,0.65,0.91,0.91)
    legend.SetLineColor(0)
    legend.SetTextSize(0.04)
    legend.SetLineWidth(0)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    if plotOnlySignal:
        legend.SetHeader(smodel)
        legend.AddEntry(hmc,"Signal Monte Carlo","EP")
        legend.AddEntry(hs,"Total signal fit", "L")
        legend.AddEntry(hgauss,"Gaussian", "L")
        legend.AddEntry(hdcb,"Double Crystal Ball", "L")
    else:
        if plotData:
            legend.AddEntry(gd,"Data","EP")
        else:
            legend.AddEntry(gd,"Data (toy)","EP")
        for pn,pp in enumerate(p):
            legend.AddEntry(hp[pn],pp,"L")
        if plotSignal:
            if sigsf>0:
                legend.AddEntry(hs,"%s, (x%.4f)"%(smodel,sigsf), "L")
            else:
                legend.AddEntry(hs,"%s"%(smodel), "L")
    return legend

ROOT.gStyle.SetOptStat(0)

colors = [ROOT.kAzure, ROOT.kRed+1, ROOT.kGreen+2, ROOT.kOrange]

temp = "tempdir"
if os.path.exists(f"{thisDir}/{temp}"):
    os.system(f"rm {thisDir}/{temp}/*")
    os.system(f"rmdir {thisDir}/{temp}/")
os.makedirs(f"{thisDir}/{temp}") 
for ch,sr in dNames:

    ## Combine to get fits
    os.system(f"cp {inDir}/card_{ch}_HTo2ZdTo2mu2x_M{mass}_ctau{ctau}_2023.txt {thisDir}/{temp}/")
    os.chdir(f"{thisDir}/{temp}/")
    os.system(f"combine -M FitDiagnostics card_{ch}_HTo2ZdTo2mu2x_M2.5_ctau1_2023.txt --saveShapes --saveWithUncertainties --saveNormalizations --cminDefaultMinimizerStrategy 0 -v 0")
    os.chdir(f"{thisDir}")

    outfile = ROOT.TFile(f"{thisDir}/{temp}/fitDiagnosticsTest.root")
    gdata = outfile.Get(f"shapes_fit_s/{ch}/data")
    htotal = outfile.Get(f"shapes_fit_s/{ch}/total")
    hsig = outfile.Get(f"shapes_fit_s/{ch}/total_signal")
    hbkg = outfile.Get(f"shapes_fit_s/{ch}/total_background")
    for n in range(0, gdata.GetN()+1):
        gdata.SetPointY(n, htotal.GetBinWidth(1)*gdata.GetPointY(n))
        gdata.SetPointEYhigh(n, htotal.GetBinWidth(1)*gdata.GetErrorYhigh(n))
        gdata.SetPointEYlow(n, htotal.GetBinWidth(1)*gdata.GetErrorYlow(n))
    htotal.Scale(htotal.GetBinWidth(1))
    hsig.Scale(hsig.GetBinWidth(1))
    hbkg.Scale(hbkg.GetBinWidth(1))
    haxis = htotal.Clone("haxis")
    haxis.SetMinimum(0)
    haxis.SetMaximum(2.0*max(gdata.GetY()))

    ## Set style on objects
    gdata.SetMarkerStyle(20)
    gdata.SetMarkerColor(ROOT.kBlack)
    gdata.SetLineWidth(2)
    htotal.SetLineWidth(3)
    htotal.SetLineColor(ROOT.kBlue+2)
    hsig.SetLineWidth(3)
    hsig.SetLineColor(ROOT.kOrange)
    hbkg.SetLineWidth(3)
    hbkg.SetLineColor(ROOT.kCyan)

    xlabel = "m_{#mu#mu} (GeV)" if "Dimuon" in sr else "m_{4#mu} (GeV)"
    haxis.GetXaxis().SetTitle(xlabel)

    ## Create canvas
    can = ROOT.TCanvas("can","",800, 600)
    can.cd()

    pads = []
    pads.append(ROOT.TPad("1","1",0,0,1,1))
    pads[0].SetTopMargin(0.08)
    pads[0].SetLeftMargin(0.15)
    pads[0].SetRightMargin(0.05)
    pads[0].Draw()

    pads[0].cd()
    haxis.Draw("AXIS")
    gdata.Draw("P,SAME")
    htotal.Draw("LINE, SAME")
    hsig.Draw("LINE, SAME")
    hbkg.Draw("LINE, SAME")

    legend = ROOT.TLegend(0.15,0.85,0.95,0.92)
    legend.SetTextSize(0.035)
    legend.SetLineWidth(1)
    legend.SetNColumns(4)
    legend.AddEntry(gdata, "Data")
    legend.AddEntry(htotal, "S+B post-fit")
    legend.AddEntry(hsig, "Signal post-fit")
    legend.AddEntry(hbkg, "Background post-fit")

    pads[0].Update()
    pads[0].RedrawAxis()
    legend.Draw()
    
    ## Search region labels
    #
    latexExtra = ROOT.TLatex()
    latexExtra.SetTextFont(42)
    latexExtra.SetTextSize(0.04)
    latexExtra.SetNDC(True)
    #
    latexExtraBold = ROOT.TLatex()
    latexExtraBold.SetTextFont(62)
    latexExtraBold.SetTextSize(0.04)
    latexExtraBold.SetNDC(True)

    ## Save canvas
    can.SaveAs("%s/%s.png"%(outDir,sr))
    can.SaveAs("%s/%s.pdf"%(outDir,sr))

os.system(f"rm {temp}/*")
os.system(f"rmdir {temp}")

