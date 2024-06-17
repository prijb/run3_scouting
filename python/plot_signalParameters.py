import ROOT
import numpy
import copy
import os,sys
from datetime import date
#import plotUtils

ROOT.gROOT.SetBatch(1)

user = os.environ.get("USER")
today= date.today().strftime("%b-%d-%Y")

doRel = True

wsname = "wfit"
thisDir = os.environ.get("PWD")
inDir  = "%s/fitResults_2022_freeze/"%thisDir

useCategorizedSignal = True
useCategorizedBackground = True

outDir = ("%s/signalParameters_"%(thisDir))+today
if not os.path.exists(outDir):
    os.makedirs(outDir)

dNames = []
dNames.append("d_FourMu_sep")
dNames.append("d_FourMu_osv")
dNames.append("d_Dimuon_lxy0p0to0p2_iso0_ptlow")
dNames.append("d_Dimuon_lxy0p0to0p2_iso0_pthigh")
dNames.append("d_Dimuon_lxy0p0to0p2_iso1_ptlow")
dNames.append("d_Dimuon_lxy0p0to0p2_iso1_pthigh")
dNames.append("d_Dimuon_lxy0p2to1p0_iso0_ptlow")
dNames.append("d_Dimuon_lxy0p2to1p0_iso0_pthigh")
dNames.append("d_Dimuon_lxy0p2to1p0_iso1_ptlow")
dNames.append("d_Dimuon_lxy0p2to1p0_iso1_pthigh")
dNames.append("d_Dimuon_lxy1p0to2p4_iso0_ptlow")
dNames.append("d_Dimuon_lxy1p0to2p4_iso0_pthigh")
dNames.append("d_Dimuon_lxy1p0to2p4_iso1_ptlow")
dNames.append("d_Dimuon_lxy1p0to2p4_iso1_pthigh")
dNames.append("d_Dimuon_lxy2p4to3p1_iso0_ptlow")
dNames.append("d_Dimuon_lxy2p4to3p1_iso0_pthigh")
dNames.append("d_Dimuon_lxy2p4to3p1_iso1_ptlow")
dNames.append("d_Dimuon_lxy2p4to3p1_iso1_pthigh")
dNames.append("d_Dimuon_lxy3p1to7p0_iso0_ptlow")
dNames.append("d_Dimuon_lxy3p1to7p0_iso0_pthigh")
dNames.append("d_Dimuon_lxy3p1to7p0_iso1_ptlow")
dNames.append("d_Dimuon_lxy3p1to7p0_iso1_pthigh")
dNames.append("d_Dimuon_lxy7p0to11p0_iso0_ptlow")
dNames.append("d_Dimuon_lxy7p0to11p0_iso0_pthigh")
dNames.append("d_Dimuon_lxy7p0to11p0_iso1_ptlow")
dNames.append("d_Dimuon_lxy7p0to11p0_iso1_pthigh")
dNames.append("d_Dimuon_lxy11p0to16p0_iso0_ptlow")
dNames.append("d_Dimuon_lxy11p0to16p0_iso0_pthigh")
dNames.append("d_Dimuon_lxy11p0to16p0_iso1_ptlow")
dNames.append("d_Dimuon_lxy11p0to16p0_iso1_pthigh")
dNames.append("d_Dimuon_lxy16p0to70p0_iso0_ptlow")
dNames.append("d_Dimuon_lxy16p0to70p0_iso0_pthigh")
dNames.append("d_Dimuon_lxy16p0to70p0_iso1_ptlow")
dNames.append("d_Dimuon_lxy16p0to70p0_iso1_pthigh")
dNames.append("d_Dimuon_lxy0p0to0p2_non-pointing")
dNames.append("d_Dimuon_lxy0p2to1p0_non-pointing")
dNames.append("d_Dimuon_lxy1p0to2p4_non-pointing")
dNames.append("d_Dimuon_lxy2p4to3p1_non-pointing")
dNames.append("d_Dimuon_lxy3p1to7p0_non-pointing")
dNames.append("d_Dimuon_lxy7p0to11p0_non-pointing")
dNames.append("d_Dimuon_lxy11p0to16p0_non-pointing")
dNames.append("d_Dimuon_lxy16p0to70p0_non-pointing")
#dNames.append("")

years = []
#years.append("2022")
#years.append("2017")
#years.append("2016APV")
#years.append("2016nonAPV")
###
years.append("allEras")

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

f = ROOT.TFile.Open("utils/signalFitParameters.root")

params = []
params.append(['gsigma', 'splines', '#sigma'])
params.append(['gmean', 'splinem', '#mu'])
params.append(['gnL', 'splinenL', 'n_{L}'])
params.append(['gnR', 'splinenR', 'n_{R}'])
params.append(['gaL', 'splineaL', 'a_{L}'])
params.append(['gaR', 'splineaR', 'a_{R}'])

colors = [ROOT.kAzure, ROOT.kRed+1, ROOT.kGreen+2, ROOT.kOrange]

for p,par in enumerate(params):

    can = ROOT.TCanvas("can","",800, 600)
    can.cd()

    pads = []
    pads.append(ROOT.TPad("1","1",0,0,1,1))
    pads[0].SetTopMargin(0.08)
    pads[0].SetLeftMargin(0.15)
    pads[0].SetRightMargin(0.05)
    pads[0].Draw()

    pads[0].cd()

    ## Create axis plot
    ymax = 0.01
    h_axis = ROOT.TH1D("h_axis","", 100, 0.0, 50.0)
    h_axis.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    h_axis.GetYaxis().SetTitle("Parameter")
    h_axis.SetMinimum(0.0)
    h_axis.SetMaximum(1.0)
    h_axis.GetYaxis().SetRangeUser(0.0,1.0)
    h_axis.GetXaxis().SetTitleSize(0.045)
    h_axis.GetXaxis().SetLabelSize(0.04)
    h_axis.GetYaxis().SetTitleSize(0.045)
    h_axis.GetYaxis().SetLabelSize(0.04)
    h_axis.GetXaxis().SetTitleOffset(1.1)
    h_axis.Draw("")

    legend = ROOT.TLegend(0.15,0.85,0.95,0.92)
    #legend.SetLineColor(0)
    legend.SetTextSize(0.04)
    legend.SetLineWidth(1)
    legend.SetNColumns(3)

    relgraph = []
    for t,ctau in enumerate([1, 10, 100]):

        gname = "%s_%smm"%(par[0], str(ctau))
        sname = "%s_%smm"%(par[1], str(ctau))
        print(gname, sname)

        graph = f.Get(gname)
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(colors[t])
        graph.SetMarkerSize(1.2)
        graph.SetLineWidth(2)

        spline = f.Get(gname)
        spline.SetLineWidth(2)
        spline.SetLineColor(colors[t])
        spline.SetLineStyle(2)

        relgraph.append(ROOT.TGraph())
        relgraph[-1].SetName(gname + "_rel")
        relgraph[-1].SetMarkerStyle(20)
        relgraph[-1].SetMarkerColor(colors[t])
        relgraph[-1].SetMarkerSize(1.2)
        relgraph[-1].SetLineWidth(2)
        for i in range(0, graph.GetN()):
            relgraph[-1].AddPoint(graph.GetX()[i], graph.GetY()[i]/graph.GetX()[i])

        #minR=0.5*min(graph.GetY())
        #maxR=1.5*max(graph.GetY())
        
        #graph.Draw("AP")
        if not t:
            ylabel = par[2]
            if doRel: 
                ylabel += "/m_{#mu#mu}"
            h_axis.GetYaxis().SetTitle(ylabel)

        if doRel:
            relgraph[-1].Draw("P, SAME")
            ymax = max(relgraph[-1].GetY()) if max(relgraph[-1].GetY()) > ymax else ymax
            legend.AddEntry(relgraph[-1], f"c#tau = {ctau} mm", "p")
        else:
            ymax = max(graph.GetY()) if max(graph.GetY()) > ymax else ymax
            graph.Draw("P, SAME")
            spline.Draw("L, SAME")
            legend.AddEntry(graph, f"c#tau = {ctau} mm", "p")

    #llabel = "M_{Z_{D}} = %s, c#tau = %s mm"%(m, t)
    #legend = getLegend(binidx,g_data,pn,hp,hs,llabel,float(m), -1, plotSignal)
    #year="all"
    #drawLabels(year,lumi,useData)
    #legend.Draw("same")
    h_axis.GetYaxis().SetRangeUser(0.0,1.4*ymax)
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
    can.SaveAs("%s/%s.png"%(outDir,par[0]))
    can.SaveAs("%s/%s.pdf"%(outDir,par[0]))

f.Close()

