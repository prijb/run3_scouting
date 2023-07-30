import ROOT
import os,sys,json
import argparse
from datetime import date    
import numpy as np
import copy
sys.path.append('utils')
import plotUtils

ROOT.gROOT.SetBatch(1)

user = os.environ.get("USER")
today= date.today().strftime("%b-%d-%Y")

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--inDir", default=os.environ.get("PWD")+"/outputHistograms_"+today, help="Choose output directory. Default: '"+os.environ.get("PWD")+"/outputHistograms_"+today+"'")
parser.add_argument("--inSamples", default=[], nargs="+", help="Choose sample(s); for all data samples in input directory, choose 'Data'")
parser.add_argument("--outDir", default=os.environ.get("PWD")+"/plots_"+today, help="Choose output directory. Default: '"+os.environ.get("PWD")+"/plots_"+today+"'")
parser.add_argument("--outSuffix", default="", help="Choose output directory. Default: ''")
parser.add_argument("--data", default=False, action="store_true", help="Plot data")
parser.add_argument("--signal", default=False, action="store_true", help="Plot signal")
parser.add_argument("--generator", default=False, action="store_true", help="Plot GEN-level histograms")
parser.add_argument("--year", default="2022", help="Year to be processes. Default: 2022")
parser.add_argument("--partialUnblinding", default=False, action="store_true", help="Have processed x% (default: x=50) of available data")
parser.add_argument("--partialUnblindingFraction", default="0.5", help="Fraction of available data that have been processed")
parser.add_argument("--dimuonMassSel", default=[], nargs="+", help="Selection on dimuon mass: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonPtSel", default=[], nargs="+", help="Selection on dimuon pT: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--fourmuonMassSel", default=[], nargs="+", help="Selection on four-muon mass: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--fourmuonPtSel", default=[], nargs="+", help="Selection on four-muon pT: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassSelForFourMuon", default=[], nargs="+", help="Selection on dimuon mass in four-muon system: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassDiffSelForFourMuon", default=["0.05"], nargs="+", help="Selection on dimuon mass difference / mean in four-muon system: first (or only) value is *upper* cut, second (optional) value is *lower* cut")
parser.add_argument("--dimuonPtSelForFourMuon", default=[], nargs="+", help="Selection on dimuon pT in four-muon system: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassSelForFourMuonOSV", default=[], nargs="+", help="Selection on dimuon mass in four-muon system from overlapping SV: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassDiffSelForFourMuonOSV", default=["0.05"], nargs="+", help="Selection on dimuon mass difference / mean in four-muon system from overlapping SV: first (or only) value is *upper* cut, second (optional) value is *lower* cut")
parser.add_argument("--dimuonPtSelForFourMuonOSV", default=[], nargs="+", help="Selection on dimuon pT in four-muon system from overlapping SV: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--lxySel", default=[], nargs="+", help="Selection on lxy: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--relaxedSVSel", default=False, action="store_true", help="Extend range for 1D histograms of SV selection variables")
parser.add_argument("--doRatio", default=False, action="store_true", help="Plot ratios for 1D histograms")
parser.add_argument("--shape", default=False, action="store_true", help="Shape normalization")
parser.add_argument("--logY", default=False, action="store_true", help="Log-scale for Y axis")
parser.add_argument("--zoomMass", default=[], nargs="+", help="Zoom mass [1] < mass < [0] GeV")
parser.add_argument("--zoomLxy", default=[], nargs="+", help="Zoom lxy [1] < lxy < [0] cm")
parser.add_argument("--zoomPixel2D", default=False, action="store_true", help="Zoom -20 <(x,y)< 20 cm")
parser.add_argument("--pdf", default=False, action="store_true", help="Output format: .pdf. Default: .png")
args = parser.parse_args()

isData   = args.data
isSignal = args.signal
samples = args.inSamples
if len(samples)<1:
    print("Please, specify list of samples to be plotted")
    exit()

for s in samples:
    if "Data" in s:
        isData = True
        break

doRatio = args.doRatio
if len(samples)<=1:
    doRatio = False

roff = 0.0
if not doRatio:
    roff=0.01

skimFraction = float(args.partialUnblindingFraction)
skimEvents = args.partialUnblinding
if not isData:
    skimEvents = False
luminosity = 35.144834218
luminosity2022B = 0.085456340
luminosity2022C = 5.070714851
luminosity2022D = 3.006332096
luminosity2022E = 5.866801170
luminosity2022F = 18.006671456
luminosity2022G = 3.108858306
if samples[0]!="Data" and len(samples)==1:
    ts = samples[0]
    if ts=="DataB":
        luminosity = luminosity2022B
    elif ts=="DataC":
        luminosity = luminosity2022C
    elif ts=="DataD":
        luminosity = luminosity2022D
    elif ts=="DataE":
        luminosity = luminosity2022E
    elif ts=="DataF":
        luminosity = luminosity2022F
    elif ts=="DataG":
        luminosity = luminosity2022G
luminosity      = 0.1*luminosity
luminosity2022B = 0.1*luminosity2022B
luminosity2022C = 0.1*luminosity2022C
luminosity2022D = 0.1*luminosity2022D
luminosity2022E = 0.1*luminosity2022E
luminosity2022F = 0.1*luminosity2022F
luminosity2022G = 0.1*luminosity2022G
if skimEvents and skimFraction>0.0:
    luminosity      = skimFraction*luminosity
    luminosity2022B = skimFraction*luminosity2022B
    luminosity2022C = skimFraction*luminosity2022C
    luminosity2022D = skimFraction*luminosity2022D
    luminosity2022E = skimFraction*luminosity2022E
    luminosity2022F = skimFraction*luminosity2022F
    luminosity2022G = skimFraction*luminosity2022G

colors = dict()
colors["Data"]  = ROOT.kBlack
colors["DataB"] = ROOT.kYellow+1
colors["DataC"] = ROOT.kOrange+1
colors["DataD"] = ROOT.kRed+1
colors["DataE"] = ROOT.kPink+1
colors["DataF"] = ROOT.kMagenta+1
colors["DataG"] = ROOT.kViolet+1
colors["signal"] = [ROOT.kBlue+1, ROOT.kAzure+1, ROOT.kCyan+1, ROOT.kTeal+1, ROOT.kGreen+1]

legnames = dict()
legnames["Data"]  = "Data"
legnames["DataB"] = "Run2022B (%.2f/fb)"%luminosity2022B
legnames["DataC"] = "Run2022C (%.2f/fb)"%luminosity2022C
legnames["DataD"] = "Run2022D (%.2f/fb)"%luminosity2022D
legnames["DataE"] = "Run2022E (%.2f/fb)"%luminosity2022E
legnames["DataF"] = "Run2022F (%.2f/fb)"%luminosity2022F
legnames["DataG"] = "Run2022G (%.2f/fb)"%luminosity2022G

for s in samples:
    if "Signal" in s:
        legnames[s] = s.replace("_"," ")

samplecol   = []
for s in samples:
    if "Signal" in s:
        samplecol.append("signal")
    else:
        samplecol.append(s)

indir  = args.inDir
outdir = args.outDir
if args.relaxedSVSel:
    outdir = outdir+"_relaxedSVselection"
if args.outSuffix!="":
    outdir = outdir+"_"+args.outSuffix
if not os.path.exists(outdir):
    os.makedirs(outdir)
os.system('cp '+os.environ.get("PWD")+'/utils/index.php '+outdir)

infiles = []
hname = "histograms"
if args.generator:
    hname = "histograms_GEN"
for s in samples:
    if not os.path.isfile("%s/%s_%s_%s_all.root"%(indir,hname,s,args.year)):
        os.system('hadd '+indir+'/'+hname+'_'+s+'_'+args.year+'_all.root $(find '+indir+' -name "'+hname+'_'+s+'*_'+args.year+'_*.root")')
    infiles.append(indir+'/'+hname+'_'+s+'_'+args.year+'_all.root')

if len(infiles)<1:
    print("No matching input file was found in %s."%indir)
    exit()

fin = ROOT.TFile.Open(infiles[0],"r")
listkeys = fin.GetListOfKeys()
size = listkeys.GetSize()
h1dn = []
h2dn = []
for i in range(0,size):
    if "TH1" in listkeys.At(i).GetClassName():
        h1dn.append(listkeys.At(i).GetName())
    elif "TH2" in listkeys.At(i).GetClassName():
        h2dn.append(listkeys.At(i).GetName())
h1d = []
h2d = []
inf = []

ncl = 1
xol = 0.0
if len(samples)>1:
    ncl=2
    xol=0.25
if len(samples)>2:
    ncl=3
    xol=0.5
leg = ROOT.TLegend(0.69-xol, 0.89-0.06*len(samples)/ncl, 0.89, 0.89)
leg.SetNColumns(ncl)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetLineWidth(0)

nSigSamples = 0
for fn,f in enumerate(infiles):
    inf.append(ROOT.TFile(f))
    if len(h1dn)>0:
        h1d.append([])
    if len(h2dn)>0:
        h2d.append([])
    for hn in h1dn:
        ht = inf[fn].Get(hn).Clone("%s_%d"%(hn,fn))
        #ht.SetDirectory(0)
        h1d[fn].append(copy.deepcopy(ht))
        if "signal" not in samplecol[fn]:
            h1d[fn][len(h1d[fn])-1].SetLineColor(colors[samplecol[fn]])
            h1d[fn][len(h1d[fn])-1].SetMarkerColor(colors[samplecol[fn]])
        else:
            h1d[fn][len(h1d[fn])-1].SetLineColor(colors[samplecol[fn]][nSigSamples])
            h1d[fn][len(h1d[fn])-1].SetMarkerColor(colors[samplecol[fn]][nSigSamples])
            nSigSamples = nSigSamples+1
    for hn in h2dn:
        ht = inf[fn].Get(hn).Clone("%s_%d"%(hn,fn))
        #ht.SetDirectory(0)
        h2d[fn].append(copy.deepcopy(ht))

    if len(h1d[fn])>0:
        if "Data" in samples[fn]:
            leg.AddEntry(h1d[fn][0], legnames[samples[fn]], "PEL")
        else:
            leg.AddEntry(h1d[fn][0], legnames[samples[fn]], "L")

# Draw histograms
unityArea = args.shape
doLogy = args.logY

# Labels
yearenergy = "%.2f fb^{-1} (%s, 13.6 TeV)"%(luminosity,args.year)
cmsExtra = "Preliminary"
if not isData:
    yearenergy = "(13.6 TeV)"
    cmsExtra = "Simulation"

drawCMSOnTop = True
#
latex = ROOT.TLatex()
latex.SetTextFont(42)
latex.SetTextAlign(31)
latex.SetTextSize(0.04)
latex.SetNDC(True)
#
latexCMS = ROOT.TLatex()
latexCMS.SetTextFont(61)
latexCMS.SetTextSize(0.04)
latexCMS.SetNDC(True)
#
latexCMSExtra = ROOT.TLatex()
latexCMSExtra.SetTextFont(52)
latexCMSExtra.SetTextSize(0.04)
latexCMSExtra.SetNDC(True)

h1dr     = []
h1dr_den = []
h_axis       = []
h_axis_ratio = []
minX         = []
maxX         = []
minY         = []
maxY         = []
line         = []
for fn in range(len(infiles)):
    h1dr.append([])
for hn,hnn in enumerate(h1dn):
    isZoom = False
    tminY = 1e100
    tmaxY = 0.0
    for fn in range(len(infiles)):
        sfSVrange = 1.0
        if args.relaxedSVSel:
            sfSVrange = 5.0
        if "_type" not in hnn:
            xmin=None
            xmax=None
            if "lxy" in hnn:
                if len(args.zoomLxy)>0:
                    xmax=float(args.zoomLxy[0])
                    isZoom = True
                if len(args.zoomLxy)>1:
                    xmin=float(args.zoomLxy[1])
                    isZoom = True
                if "selass" in hnn or "dimuon" in hnn:
                    if len(args.lxySel)>0:
                        xmin=float(args.lxySel[0])
                        isZoom = True
                    if len(args.lxySel)>1:
                        xmax=float(args.lxySel[1])
                        isZoom = True
            if "mass" in hnn:
                if len(args.zoomMass)>0:
                    xmax=float(args.zoomMass[0])
                    isZoom = True
                if len(args.zoomMass)>1:
                    xmin=float(args.zoomMass[1])
                    isZoom = True
                if "dimuon" in hnn:
                    if len(args.dimuonMassSel)>0:
                        xmin=float(args.dimuonMassSel[0])
                        isZoom = True
                    if len(args.dimuonMassSel)>1:
                        xmax=float(args.dimuonMassSel[1])
                        isZoom = True
            if "sv" in hnn and ("xerr" in hnn or "yerr" in hnn):
                xmax=0.05*sfSVrange
            if "sv" in hnn and "zerr" in hnn:
                xmax=0.10*sfSVrange
            if "sv" in hnn and "chi2ndof" in hnn and not "muon" in hnn:
                xmax=5.0*sfSVrange
            if not isZoom:
                plotUtils.PutUnderflowInFirstBin(h1d[fn][hn],xmin)
                plotUtils.PutOverflowInLastBin(h1d[fn][hn],xmax)
        h1dr[fn].append(h1d[fn][hn].Clone("%s_ratio"%hnn))
        if fn==0:
            h1dr_den.append(h1d[fn][hn].Clone("%s_denominator"%hnn))
            if samples[fn]=="Data":
                for b in range(1,h1dr[fn][hn].GetNbinsX()+1):
                    h1dr_den[hn].SetBinError(b,0.0)
            if xmin!=None:
                tb = h1d[fn][hn].GetXaxis().FindBin(xmin)
                h1d [fn][hn].GetXaxis().SetRangeUser(h1d[fn][hn].GetXaxis().GetBinLowEdge(tb),h1d[fn][hn].GetXaxis().GetBinUpEdge(h1d[fn][hn].GetNbinsX()))
                h1dr[fn][hn].GetXaxis().SetRangeUser(h1d[fn][hn].GetXaxis().GetBinLowEdge(tb),h1d[fn][hn].GetXaxis().GetBinUpEdge(h1d[fn][hn].GetNbinsX()))
            if xmax!=None:
                tb = h1d[fn][hn].GetXaxis().FindBin(xmax)
                if h1d[fn][hn].Integral(tb, -1)<=0.0 or (isZoom and not xmax>h1d[fn][hn].GetXaxis().GetBinLowEdge(tb)):
                    tb = tb-1
                h1d [fn][hn].GetXaxis().SetRangeUser(h1d[fn][hn].GetXaxis().GetBinLowEdge(1),h1d[fn][hn].GetXaxis().GetBinUpEdge(tb))
                h1dr[fn][hn].GetXaxis().SetRangeUser(h1d[fn][hn].GetXaxis().GetBinLowEdge(1),h1d[fn][hn].GetXaxis().GetBinUpEdge(tb))
            if xmin!=None and xmax!=None:
                tbm = h1d[fn][hn].GetXaxis().FindBin(xmin)
                tbM = h1d[fn][hn].GetXaxis().FindBin(xmax)
                if h1d[fn][hn].Integral(tbM, -1)<=0.0 or (isZoom and not xmax>h1d[fn][hn].GetXaxis().GetBinLowEdge(tbM)):
                    tbM = tbM-1
                h1d [fn][hn].GetXaxis().SetRangeUser(h1d[fn][hn].GetXaxis().GetBinLowEdge(tbm),h1d[fn][hn].GetXaxis().GetBinUpEdge(tbM))
                h1dr[fn][hn].GetXaxis().SetRangeUser(h1d[fn][hn].GetXaxis().GetBinLowEdge(tbm),h1d[fn][hn].GetXaxis().GetBinUpEdge(tbM))
            if "Data" in samples[fn]:
                h1d [fn][hn].SetBinErrorOption(ROOT.TH1.kPoisson)
                h1dr[fn][hn].SetBinErrorOption(ROOT.TH1.kPoisson)
            h1d[fn][hn].GetYaxis().SetLabelSize(0.025)
            h1d[fn][hn].GetYaxis().SetMaxDigits(3)
            h1d[fn][hn].SetLineWidth(2)
            ytitle = h1d[fn][hn].GetYaxis().GetTitle()
            if not unityArea:
                tmaxY = max(tmaxY, h1d[fn][hn].GetMaximum())
                tminY = 0.0
                if doLogy:
                    tminY = 0.9
        if unityArea:
            ytitle = ytitle.replace("Events", "Fraction of events")
            integral = h1d[fn][hn].Integral(0,-1)
            if integral>0.0:
                h1d[fn][hn].Scale(1.0/integral)
                h1dr[fn][hn].Scale(1.0/integral)
                if fn==0:
                    h1dr_den[hn].Scale(1.0/integral)
            tmaxY = max(tmaxY, h1d[fn][hn].GetMaximum())
            if doLogy:
                minb = 1e100
                for b in range(1,h1d[fn][hn].GetNbinsX()+1):
                    if h1d[fn][hn].GetBinContent(b)<minb and h1d[fn][hn].GetBinContent(b)>0.0:
                        minb = h1d[fn][hn].GetBinContent(b)
                if integral>0.0:
                    tminY = min(tminY, max(0.1*minb,0.1/integral))
                else:
                    tminY = min(tminY, 0.1)
                    tmaxY = max(tmaxY, 1.0)
        h1dr[fn][hn].Divide(h1dr_den[hn])
        if "hsv" in hnn and not "mind" in hnn and not "maxd" in hnn:
            ytitle = ytitle.replace("Events", "Number of SVs")
            ytitle = ytitle.replace("events", "SVs")
        if "hmuon" in hnn and hnn!="hmuon_mindr" and hnn!="hmuon_maxdr":
            ytitle = ytitle.replace("Events", "Number of muons")
            ytitle = ytitle.replace("events", "muons")
        h1d[fn][hn].GetYaxis().SetTitle(ytitle)
        if fn==0:
            h_axis.append(h1d[fn][hn].Clone("axis_%s"%hnn))
            h_axis[hn].Reset("ICE")
            h_axis[hn].GetYaxis().SetTitleSize(0.04)
            h_axis[hn].GetYaxis().SetTitleOffset(1.35)
            h_axis[hn].GetXaxis().SetTitleSize(0.04)
            h_axis[hn].GetXaxis().SetTitleOffset(1.25-roff*10.0)
            h_axis[hn].GetXaxis().SetLabelSize(0.035)
            h_axis[hn].GetYaxis().SetLabelSize(0.035)
            h_axis_ratio.append(h_axis[hn].Clone("axis_ratio_%s"%hnn))
            if xmax==None:
                maxX.append(h_axis[hn].GetXaxis().GetBinUpEdge(plotUtils.GetLastBin(h_axis[hn])))
            else:
                maxX.append(min(xmax, h_axis[hn].GetXaxis().GetBinUpEdge(plotUtils.GetLastBin(h_axis[hn]))))
            if xmin==None:
                minX.append(h_axis[hn].GetXaxis().GetBinLowEdge(plotUtils.GetFirstBin(h_axis[hn])))
            else:
                minX.append(max(xmin, h_axis[hn].GetXaxis().GetBinLowEdge(plotUtils.GetFirstBin(h_axis[hn]))))
            line.append(ROOT.TLine(minX[hn], 1.0, maxX[hn], 1.0))
            minY.append(tminY)
            maxY.append(tmaxY)
        else:
            if tminY<minY[hn]:
                minY[hn]=tminY
            if tmaxY>maxY[hn]:
                maxY[hn]=tmaxY

ROOT.gStyle.SetOptStat(0)
can = ROOT.TCanvas("can","",600, 600)
minR = []
maxR = []
for hn,hnn in enumerate(h1dn):
    pads = []
    if doRatio:
        pads.append(ROOT.TPad("1","1",0.0,0.18,1.0,1.0))
        pads.append(ROOT.TPad("2","2",0.0,0.0,1.0,0.19))
        pads[0].SetTopMargin(0.08)
        pads[0].SetBottomMargin(0.13)
        pads[0].SetLeftMargin(0.10)
        pads[0].SetRightMargin(0.05)
        pads[1].SetLeftMargin(0.10)
        pads[1].SetRightMargin(0.05)
    else:
        pads.append(ROOT.TPad("1","1",0,0,1,1))
        pads[0].SetTopMargin(0.08)
        pads[0].SetLeftMargin(0.10)
        pads[0].SetRightMargin(0.05)
    if doRatio:
        minR.append(0.0)
        maxR.append(2.0)
        ty = np.array([])
        tmax=maxR[hn]
        for fn in range(len(infiles)):
            tg = ROOT.TGraph(h1dr[fn][hn])
            ty = tg.GetY()
            if len(ty)>0:
                tmax = np.amax(ty)
            if tmax>maxR[hn]:
                maxR[hn]=tmax*1.05
            if maxR[hn]>5.0:
                minR[hn]=0.1
        h_axis_ratio[hn].GetYaxis().SetRangeUser(minR[hn],maxR[hn])
        h_axis_ratio[hn].SetMinimum(minR[hn])
        h_axis_ratio[hn].SetMaximum(maxR[hn])
        h_axis_ratio[hn].SetTitle(";;Ratio")
        h_axis_ratio[hn].GetYaxis().SetTitleSize(0.16)
        h_axis_ratio[hn].GetYaxis().SetTitleOffset(0.25)
        if doLogy:
            h_axis_ratio[hn].GetYaxis().SetTitleOffset(0.3)
        h_axis_ratio[hn].GetYaxis().SetNdivisions(505)
        h_axis_ratio[hn].GetYaxis().SetLabelSize(0.12)
        h_axis_ratio[hn].GetYaxis().CenterTitle()
        h_axis_ratio[hn].GetYaxis().SetTickLength(0.02)
        h_axis_ratio[hn].GetXaxis().SetLabelSize(0.0)
        h_axis_ratio[hn].GetXaxis().SetTitle("")
        h_axis_ratio[hn].GetXaxis().SetTickSize(0.06)

        pads[0].Draw()
        pads[1].Draw()
        pads[1].cd()
        if maxR[hn]>5.0:
            pads[1].SetLogy()
        pads[1].SetTickx()
        pads[1].SetTicky()
        h_axis_ratio[hn].Draw("")
        for fn in range(1,len(infiles)):
            if "Data" in samples[fn]:
                h1dr[fn][hn].Draw("SAME,P,E")
            else:
                h1dr[fn][hn].Draw("SAME,HIST,E")
        line[hn].SetLineStyle(2)
        line[hn].SetLineColor(colors[samplecol[0]])
        line[hn].SetLineWidth(1)
        line[hn].Draw("SAME")
        pads[1].Modified()
        pads[1].Update()

    else:
        pads[0].Draw()

    pads[0].cd()
    pads[0].SetTicky()
    if doLogy:
        pads[0].SetLogy()
    #h_axis[hn].GetYaxis().SetMoreLogLabels()
    if ("hits" in hnn) or ("layers" in hnn) or ("nmu" in hnn) or h1d[0][hn].GetXaxis().GetBinLowEdge(1)<0.0:
        maxY[hn] = 1.25*maxY[hn]
        if doLogy:
            maxY[hn] = 5.0*maxY[hn]
    if doLogy:
        maxY[hn] = 10.0*maxY[hn]
    else:
        maxY[hn] = 2.0*maxY[hn]
    print(minY[hn], maxY[hn])
    h_axis[hn].GetYaxis().SetRangeUser(minY[hn], maxY[hn])
    h_axis[hn].SetMinimum(minY[hn])
    h_axis[hn].SetMaximum(maxY[hn])
    h_axis[hn].Draw("")
    for fn in range(len(infiles)):
        if "Data" in samples[fn]:
            h1d[fn][hn].Draw("PE,same")
        else:
            h1d[fn][hn].Draw("hist,same")
    leg.Draw("same")
    pads[0].RedrawAxis()
    pads[0].Update()
    can.cd()
    latex.DrawLatex(0.95, 0.945-roff, yearenergy);
    if drawCMSOnTop:
        latexCMS.DrawLatex(0.11,0.945-roff,"CMS");
        latexCMSExtra.DrawLatex(0.21,0.945-roff, cmsExtra);
    else:
        latexCMS.DrawLatex(0.14,0.875,"CMS");
        latexCMSExtra.DrawLatex(0.14,0.825, cmsExtra);
    outname = "%s"%hnn
    if "h_" in outname:
        outname = outname.replace("h_","")
    elif "hsv" in outname or "hmuon" in outname or "hsel" in outname:
        outname = outname.replace("hsv","sv")
        outname = outname.replace("hmuon","muon")
        outname = outname.replace("hsel","sel")
    elif "hdi" in outname or "hfour" in outname:
        outname = outname.replace("hdi","di")
        outname = outname.replace("hfour","four")
    if unityArea:
        outname = outname+"_unitArea"
    if doLogy:
        outname = outname+"_logY"
    can.SaveAs("%s/%s.png"%(outdir,outname))
    if args.pdf:
        can.SaveAs("%s/%s.pdf"%(outdir,outname))
    can.Clear()
    del pads

for hn,hnn in enumerate(h2dn):
    for fn in range(len(infiles)):
        h = h2d[fn][hn].Clone()
        maxc = 100.0
        if args.zoomPixel2D:
            maxc =  20.0
        if len(args.lxySel)>1 and float(args.lxySel[1])<maxc:
            maxc = float(args.lxySel[1])
        if len(args.zoomLxy)>0 and float(args.zoomLxy[0])<maxc:
            maxc = float(args.zoomLxy[0])
        minc = -1.0*maxc
        h.GetXaxis().SetRangeUser(minc, maxc)
        h.GetYaxis().SetRangeUser(minc, maxc)
        can.cd()
        ROOT.gPad.SetLogy(0)
        if doLogy:
            ROOT.gPad.SetLogz()
        h.GetYaxis().SetLabelSize(0.025)
        h.GetYaxis().SetMaxDigits(3)
        ztitle = h.GetZaxis().GetTitle()
        integral = h.Integral(0,-1,0,-1)
        if unityArea:
            ztitle = ztitle.replace("Events", "Fraction of events")
            ztitle = ztitle.replace("Number", "Fraction")
            if integral>0:
                h.Scale(1.0/integral)
            if not doLogy:
                h.GetZaxis().SetRangeUser(0.0,1.1*h.GetMaximum())
            else:
                if integral>0:
                    h.GetZaxis().SetRangeUser(1.0/integral,1.1*h.GetMaximum())
                else:
                    h.GetZaxis().SetRangeUser(1e-1,1.1*h.GetMaximum())
        if "nmuons" in h.GetName():
            ztitle = ztitle.replace("Events", "Number of muons")
            ztitle = ztitle.replace("events", "muons")
        h.GetZaxis().SetTitle(ztitle)
        can.cd()
        if unityArea:
            ROOT.gStyle.SetPaintTextFormat(".2f");
        h.Draw("colz")
        can.Update()
        palette = h.GetListOfFunctions().FindObject("palette")
        if palette:
            palette.SetX2NDC(0.925)
        can.Update()
        h.GetZaxis().SetTitle(ztitle)
        h.GetZaxis().SetLabelSize(0.025)
        h.GetZaxis().SetMaxDigits(2)
        if not doLogy:
            h.GetZaxis().SetRangeUser(0.0,1.1*h.GetMaximum())
        else:
            if integral>0.0:
                h.GetZaxis().SetRangeUser(1.0/integral,2.0*h.GetMaximum())
            else:
                h.GetZaxis().SetRangeUser(1e-1,1.1*h.GetMaximum())
        latex.DrawLatex(0.90, 0.91, yearenergy);
        if drawCMSOnTop:
            latexCMS.DrawLatex(0.1,0.91,"CMS");
            latexCMSExtra.DrawLatex(0.19,0.91, cmsExtra);
        else:
            latexCMS.DrawLatex(0.13,0.86,"CMS");
            latexCMSExtra.DrawLatex(0.13,0.815, cmsExtra);
        ROOT.gPad.RedrawAxis()
        outname = "%s_%s"%(hnn,samples[fn])
        if "h_" in outname:
            outname = outname.replace("h_","")
        elif "hsv" in outname:
            outname = outname.replace("hsv","sv")
        if unityArea:
            outname = outname+"_unitArea"
        if doLogy:
            outname = outname+"_logZ"
        can.SaveAs("%s/%s.png"%(outdir,outname))
        if args.pdf:
            can.SaveAs("%s/%s.pdf"%(outdir,outname))
        can.Clear()

del can
