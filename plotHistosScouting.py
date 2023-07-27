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
parser.add_argument("--lxySelForFourMuon", default=[], nargs="+", help="Selection on lxy: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--shape", default=False, action="store_true", help="Shape normalization")
parser.add_argument("--logY", default=False, action="store_true", help="Log-scale for Y axis")
parser.add_argument("--zoomMass", default=False, action="store_true", help="Zoom mass < 35 GeV")
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
legnames["DataB"] = "Run2022B"
legnames["DataC"] = "Run2022C"
legnames["DataD"] = "Run2022D"
legnames["DataE"] = "Run2022E"
legnames["DataF"] = "Run2022F"
legnames["DataG"] = "Run2022G"

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

leg = ROOT.TLegend(0.6, 0.875-0.06*len(samples), 0.875, 0.875)
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
            h1d[fn][len(h1d[fn])-1].SetLineColor(colors[samplecol[fn][nSigSamples]])
            h1d[fn][len(h1d[fn])-1].SetMarkerColor(colors[samplecol[fn][nSigSamples]])
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

skimFraction = float(args.partialUnblindingFraction)
skimEvents = args.partialUnblinding
if not isData:
    skimEvents = False
luminosity = 35.144834218
if samples[0]!="Data" and len(samples)==1:
    ts = samples[0]
    if ts=="DataB":
        luminosity = 0.085456340
    elif ts=="DataC":
        luminosity = 5.070714851
    elif ts=="DataD":
        luminosity = 3.006332096
    elif ts=="DataE":
        luminosity = 5.866801170
    elif ts=="DataF":
        luminosity = 18.006671456
    elif ts=="DataG":
        luminosity = 3.108858306
luminosity = 0.1*luminosity
if skimEvents and skimFraction>0.0:
    luminosity = skimFraction*luminosity

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

ROOT.gStyle.SetOptStat(0)
can = ROOT.TCanvas("can","",600, 600)
for hn,hnn in enumerate(h1dn):
    for fn in range(len(infiles)):
        h = h1d[fn][hn].Clone()
        #print(fn, h1d[fn][hn].GetName())
        if "_type" not in hnn:
            xmin=None
            xmax=None
            if "fourmuon_lxy" in hnn:
                xmax=5.0
            if "fourmuon_mass" in hnn or "dimuon_mass" in hnn and args.zoomMass:
                xmax=35.0
            plotUtils.PutUnderflowInFirstBin(h1d[fn][hn],xmin)
            plotUtils.PutOverflowInLastBin(h1d[fn][hn],xmax)
        if fn==0:
            if xmin!=None:
                h1d[fn][hn].GetXaxis().SetRangeUser(h1d[fn][hn].GetXaxis().GetBinLowEdge(h1d[fn][hn].GetXaxis().FindBin(xmin)),h1d[fn][hn].GetXaxis().GetBinUpEdge(h1d[fn][hn].GetNbinsX()))
            if xmax!=None:
                h1d[fn][hn].GetXaxis().SetRangeUser(h1d[fn][hn].GetXaxis().GetBinLowEdge(1),h1d[fn][hn].GetXaxis().GetBinUpEdge(h1d[fn][hn].GetXaxis().FindBin(xmax)))
            if isData:
                h1d[fn][hn].SetBinErrorOption(ROOT.TH1.kPoisson)
            h1d[fn][hn].GetYaxis().SetLabelSize(0.025)
            h1d[fn][hn].GetYaxis().SetMaxDigits(3)
            h1d[fn][hn].GetYaxis().SetRangeUser(0.0, 1.1*h1d[fn][hn].GetMaximum())
            h1d[fn][hn].SetLineWidth(2)
            ytitle = h1d[fn][hn].GetYaxis().GetTitle()
            if doLogy:
                h1d[fn][hn].GetYaxis().SetRangeUser(0.9, 2.0*h1d[fn][hn].GetMaximum())
        if unityArea:
            ytitle = ytitle.replace("Events", "Fraction of events")
            integral = h1d[fn][hn].Integral(0,-1)
            if h1d[fn][hn].Integral(0,-1)>0:
                h1d[fn][hn].Scale(1.0/h1d[fn][hn].Integral(0,-1))
            h1d[fn][hn].GetYaxis().SetRangeUser(0.0,1.1*h1d[fn][hn].GetMaximum())
            if doLogy:
                minb = 1e9
                for b in range(1,h1d[fn][hn].GetNbinsX()+1):
                    if h1d[fn][hn].GetBinContent(b)<minb:
                        minb = h1d[fn][hn].GetBinContent(b)
                if integral>0.0:
                    h1d[fn][hn].GetYaxis().SetRangeUser(max(0.9*minb,0.9/integral), 2.0*h1d[fn][hn].GetMaximum())
                else:
                    h1d[fn][hn].GetYaxis().SetRangeUser(1e-3, 1.0)
        if "hsv" in hnn and not "mind" in hnn and not "maxd" in hnn:
            ytitle = ytitle.replace("Events", "Number of SVs")
            ytitle = ytitle.replace("events", "SVs")
        if "hmuon" in hnn and hnn!="hmuon_mindr" and hnn!="hmuon_maxdr":
            ytitle = ytitle.replace("Events", "Number of muons")
            ytitle = ytitle.replace("events", "muons")
        h1d[fn][hn].GetYaxis().SetTitle(ytitle)
        can.cd()
        if fn==0:
            if doLogy:
                ROOT.gPad.SetLogy()
            if "Data" in samples[fn]:
                #print("Drawing %s, %d"%(hnn,fn))
                h1d[fn][hn].Draw("PE")
            else:
                #print("Drawing %s, %d"%(hnn,fn))
                h1d[fn][hn].Draw("hist")
        else:
            #print("Drawing %s, %d"%(hnn,fn))
            if "Data" in samples[fn]:
                h1d[fn][hn].Draw("PE,same")
            else:
                h1d[fn][hn].Draw("hist,same")                
    leg.Draw("same")
    can.cd()
    can.Update()
    latex.DrawLatex(0.90, 0.91, yearenergy);
    if drawCMSOnTop:
        latexCMS.DrawLatex(0.1,0.91,"CMS");
        latexCMSExtra.DrawLatex(0.19,0.91, cmsExtra);
    else:
        latexCMS.DrawLatex(0.13,0.86,"CMS");
        latexCMSExtra.DrawLatex(0.13,0.815, cmsExtra);
    ROOT.gPad.RedrawAxis()
    can.Update()
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

for hn,hnn in enumerate(h2dn):
    for fn in range(len(infiles)):
        h = h2d[fn][hn].Clone()
        can.cd()
        ROOT.gPad.SetLogy(0)
        h.GetYaxis().SetLabelSize(0.025)
        h.GetYaxis().SetMaxDigits(3)
        ztitle = h.GetZaxis().GetTitle()
        if unityArea:
            ztitle = ztitle.replace("Events", "Fraction of events [%]")
            if h.Integral(0,-1,0,-1)>0:
                h.Scale(100.0/h.Integral(0,-1,0,-1))
            h.GetZaxis().SetRangeUser(0.0,1.1*h.GetMaximum())
        if "nmuons" in h.GetName():
            ztitle = ztitle.replace("Events", "Number of muons")
            ztitle = ztitle.replace("events", "muons")
        h.GetZaxis().SetTitle(ztitle)
        can.cd()
        if unityArea:
            ROOT.gStyle.SetPaintTextFormat(".2f");
        h.Draw("colz,text")
        can.Update()
        palette = h.GetListOfFunctions().FindObject("palette")
        if palette:
            palette.SetX2NDC(0.925)
        can.Update()
        h.GetZaxis().SetLabelSize(0.025)
        h.GetZaxis().SetMaxDigits(2)
        h.GetZaxis().SetRangeUser(0.0, 1.1*h.GetMaximum())
        latex.DrawLatex(0.90, 0.91, yearenergy);
        if drawCMSOnTop:
            latexCMS.DrawLatex(0.1,0.91,"CMS");
            latexCMSExtra.DrawLatex(0.19,0.91, cmsExtra);
        else:
            latexCMS.DrawLatex(0.13,0.86,"CMS");
            latexCMSExtra.DrawLatex(0.13,0.815, cmsExtra);
        ROOT.gPad.RedrawAxis()
        outname = "%s_%s"%(h.GetName(),samples[fn])
        if "h_" in outname:
            outname = outname.replace("h_","")
        if unityArea:
            outname = outname+"_unitArea"
        can.SaveAs("%s/%s.png"%(outdir,outname))
        if args.pdf:
            can.SaveAs("%s/%s.pdf"%(outdir,outname))
        can.Clear()

del can
