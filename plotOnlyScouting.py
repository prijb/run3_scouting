import ROOT
import os,sys,json
from datetime import date    
import numpy as np
import copy
sys.path.append('utils')
import plotUtils

ROOT.gROOT.SetBatch(1)

isData = False
indir = sys.argv[1]
infiles = []
if len(sys.argv)==2:
    if not os.path.isfile("%s/histograms_data_all.root"%indir):
        os.system('hadd '+indir+'/histograms_data_all.root $(find '+indir+' -name "histograms_*Data*.root")')
    infiles.append("%s/histograms_data_all.root"%indir)
    isData = True
else:
    for fn in range(2,len(sys.argv)):
        if os.path.isfile("%s/histograms_%s.root"%(indir,sys.argv[fn])):
            infiles.append("%s/histograms_%s.root"%(indir,sys.argv[fn]))
        if "Data" in sys.argv[fn] or "data" in sys.argv[fn]:
            isData=True

if len(infiles)<1:
    exit()

user = os.environ.get("USER")
today= date.today().strftime("%b-%d-%Y")
outdir = '/home/users/%s/public_html/Run3Scouting/plotsOnly_%s'%(user,today)
#outdir = indir
if not os.path.exists(outdir):
    os.makedirs(outdir)
os.system('cp '+os.environ.get("PWD")+'/utils/index.php '+outdir)

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
colors = [1,2,4,6]
if not isData:
    colors = [2,4,6]
print(len(infiles))
for fn,f in enumerate(infiles):
    inf.append(ROOT.TFile(f))
    h1d.append([])
    h2d.append([])
    for hn in h1dn:
        print(f, hn)
        ht = inf[fn].Get(hn).Clone("%s_%d"%(hn,fn))
        #ht.SetDirectory(0)
        h1d[fn].append(copy.deepcopy(ht))
        h1d[fn][len(h1d[fn])-1].SetLineColor(colors[fn])
        h1d[fn][len(h1d[fn])-1].SetMarkerColor(colors[fn])

# Draw histograms
unityArea = True
doLogy = True

# Labels
yearenergy = "1.46 fb^{-1} (Run2022D, 13.6 TeV)"
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
        print(fn, h1d[fn][hn].GetName())
        if "_type" not in hnn:
            xmin=None
            xmax=None
            if "fourmuon_lxy" in hnn:
                xmax=5.0
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
                h1d[fn][hn].GetYaxis().SetRangeUser(max(0.9*minb,0.9/integral), 2.0*h1d[fn][hn].GetMaximum())
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
            if isData:
                print("Drawing %s, %d"%(hnn,fn))
                h1d[fn][hn].Draw("PE")
            else:
                print("Drawing %s, %d"%(hnn,fn))
                h1d[fn][hn].Draw("hist")
        else:
            print("Drawing %s, %d"%(hnn,fn))
            h1d[fn][hn].Draw("hist,same")
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
    elif "hsv" in outname or "hmuon" in outname:
        outname = outname.replace("hsv","sv")
        outname = outname.replace("hmuon","muon")
    elif "hdi" in outname or "hfour" in outname:
        outname = outname.replace("hdi","di")
        outname = outname.replace("hfour","four")
    can.SaveAs("%s/%s.png"%(outdir,outname))
    can.Clear()

#for h in h2d:
#    can.cd()
#    ROOT.gPad.SetLogy(0)
#    h.GetYaxis().SetLabelSize(0.025)
#    h.GetYaxis().SetMaxDigits(3)
#    ztitle = h.GetZaxis().GetTitle()
#    if unityArea:
#        ztitle = ztitle.replace("Events", "Fraction of events [%]")
#        if h.Integral(0,-1,0,-1)>0:
#            h.Scale(100.0/h.Integral(0,-1,0,-1))
#        h.GetZaxis().SetRangeUser(0.0,1.1*h.GetMaximum())
#    if "nmuons" in h.GetName():
#        ztitle = ztitle.replace("Events", "Number of muons")
#        ztitle = ztitle.replace("events", "muons")
#    h.GetZaxis().SetTitle(ztitle)
#    can.cd()
#    if unityArea:
#        ROOT.gStyle.SetPaintTextFormat(".2f");
#        h.Draw("colz,text")
#    can.Update()
#    palette = h.GetListOfFunctions().FindObject("palette")
#    if palette:
#        palette.SetX2NDC(0.925)
#    can.Update()
#    h.GetZaxis().SetLabelSize(0.025)
#    h.GetZaxis().SetMaxDigits(2)
#    h.GetZaxis().SetRangeUser(0.0, 1.1*h.GetMaximum())
#    latex.DrawLatex(0.90, 0.91, yearenergy);
#    if drawCMSOnTop:
#        latexCMS.DrawLatex(0.1,0.91,"CMS");
#        latexCMSExtra.DrawLatex(0.19,0.91, cmsExtra);
#    else:
#        latexCMS.DrawLatex(0.13,0.86,"CMS");
#        latexCMSExtra.DrawLatex(0.13,0.815, cmsExtra);
#    ROOT.gPad.RedrawAxis()
#    outname = "%s"%h.GetName()
#    if "h_" in outname:
#        outname = outname.replace("h_","")
#    can.SaveAs("%s/%s.png"%(outdir,outname))
#    can.Clear()

del can
