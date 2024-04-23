import os,sys,math
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(1)
drawObserved = True
drawPoints = True
maskSMResonances = True
typeOfLimit = "BRH" # "r", "xsec", "xsecBR" "BRH" 
xsec = 1.0 # in pb, used to normalize the MC
xsec_h = 59.8 # higgs cross section in pb at 13.6 GeV, used to normalize the MC
scaleToFullLumi = True
doRatio = True
compare = True

model = sys.argv[1]
ctau = sys.argv[2]


if typeOfLimit=="r":
    ylabel = "95% CL upper limit on #sigma/#sigma_{theory}"
elif typeOfLimit=="xsec":
    ylabel = "95% CL upper limit on #sigma"
    if model=="HTo2ZdTo2mu2x":
        ylabel = "95% CL upper limit on #sigma(h#rightarrowZ_{D}Z_{D}) [pb]"
elif typeOfLimit=="BRH":
    ylabel = "95% CL upper limit on Br(h#rightarrowXX)"
    if model=="HTo2ZdTo2mu2x":
        ylabel = "95% CL upper limit on Br(h#rightarrowZ_{D}Z_{D})"
    if model=="DarkShowers":
        ylabel = "95% CL upper limit on Br(h#rightarrow#Psi#Psi)"
elif typeOfLimit=="xsecBR":
    ylabel = "95% CL upper limit on #sigmaxBR"
    if model=="HTo2ZdTo2mu2x":
        ylabel = "95% CL upper limit on #sigma(h#rightarrowZ_{D}Z_{D})xB(Z_{D}#rightarrow#mu#mu) [pb]"

fins = []
for i in range(4, len(sys.argv)):
    fins.append(sys.argv[i])

legtext = sys.argv[3]
legnames = []
for entry in legtext.split(':'):
    legnames.append(entry)

sets = []
colors = [ROOT.kRed, ROOT.kGreen+2]
for fname in fins:
    fin = open(fname, 'r')
    massl = []
    expl  = []
    for l in fin.readlines():
        if l.startswith("#"):
            continue
        ls = l.split(",")
        if ls[0]!=model:
            continue
        if int(ls[2])!=int(ctau):
            continue
        scale = 1.0
        if typeOfLimit=="xsec":
            scale = xsec
        if typeOfLimit=="BRH":
            scale = xsec/xsec_h
        elif typeOfLimit=="xsecBR":
            if model=="HTo2ZdTo2mu2x":
                BR = 1.0
                with open('data/hahm-mass_brs.txt', 'r') as f:
                    for line in f.readlines():
                        if "mZd" in line: continue
                        brs = line.split('\t')
                        while '' in brs:
                            brs.remove('')
                        brs[-1] = brs[-1].strip('\n')
                        if (abs(float(brs[0])-float(ls[1])) < 0.0001):
                            BR = float(brs[2])
                            break
            scale = xsec*BR
            print(float(ls[1]), float(ls[4]), BR)
        if scaleToFullLumi:
            scale = scale / math.sqrt(10.0)
        if float(ctau)>10 and float(ls[1]) < 1.5: continue
        massl.append(float(ls[1]))
        expl .append(scale*float(ls[4]))
    fin.close()
    massv = np.array(massl,"d")
    expv  = np.array(expl ,"d")

    print(massv, expv)
    sets.append(ROOT.TGraph(len(massv),massv,expv))
    sets[-1].SetMarkerStyle(20)
    sets[-1].SetLineWidth(2)

## Compare with Run 2 limits

cfile = None
if compare:
    print("> Request to open files to compare")
    if model=="HTo2ZdTo2mu2x":
        if typeOfLimit=="BRH":
            if ctau=="1":
                cfile = ROOT.TFile("data/run-2/scouting/HEPData-ins1997201-v2-Figure_8a.root")
                cgexp = cfile.Get("Figure 8a/Graph1D_y2")
            if ctau=="100":
                cfile = ROOT.TFile("data/run-2/scouting/HEPData-ins1997201-v2-Figure_8b.root")
                cgexp = cfile.Get("Figure 8b/Graph1D_y2")
    cgexp.SetLineColor(ROOT.kBlue)
    cgexp.SetMarkerColor(ROOT.kBlue)
    cgexp.SetMarkerStyle(20)
    cgexp.SetLineStyle(2)
    cgexp.SetLineWidth(2)
    print("> Files accessed")


ratio = ROOT.TGraph()
ratio0 = ROOT.TGraph()
ratio1 = ROOT.TGraph()
for i in range(0, len(sets[0].GetX())):
    y = sets[1].GetY()[i]/sets[0].GetY()[i]
    ratio.AddPoint(sets[0].GetX()[i], y)
    if sets[0].GetX()[i] >= 0.7:
        ratio0.AddPoint(sets[0].GetX()[i], sets[0].GetY()[i]/cgexp.Eval(sets[0].GetX()[i]))
        ratio1.AddPoint(sets[1].GetX()[i], sets[1].GetY()[i]/cgexp.Eval(sets[0].GetX()[i]))
ratio.SetLineWidth(2)
ratio0.SetLineWidth(2)
ratio0.SetLineColor(colors[0])
ratio1.SetLineWidth(2)
ratio1.SetLineColor(colors[1])


leg = ROOT.TLegend(0.3,0.65,0.89,0.85)
leg.SetFillColor(0)
leg.SetTextSize(0.03)
modelLeg = model
if model=="HTo2ZdTo2mu2x":
    modelLeg="h#rightarrow Z_{D}Z_{D}"
leg.SetHeader("%s (c#tau = %s mm)"%(modelLeg,ctau))
for l,limit in enumerate(sets):
    leg.AddEntry(limit,legnames[l],"L")
if compare:
    leg.AddEntry(cgexp,"Run 2 result (101 fb^{-1})","L")

if doRatio and compare:
    leg2 = ROOT.TLegend(0.13,0.65,0.5,0.89)
    leg2.SetFillColor(0)
    leg2.SetLineWidth(0)
    leg2.SetTextSize(0.08)
    leg2.SetMargin(0.1)
    leg2.AddEntry(ratio,"Run 3 data (run 2 grid) / Run 3 data (run 3 grid)","L")
    leg2.AddEntry(ratio0,"Run 3 data (run 3 grid) / Run 2 result","L")
    leg2.AddEntry(ratio1,"Run 3 data (run 2 grid) / Run 2 result","L")

## Axis histogram
miny=5e-6
maxy=1e-1
haxis = ROOT.TH1D("haxis","",10,sets[0].GetX()[0]-0.05*massv[0],sets[0].GetX()[len(sets[0].GetX())-1]+0.05*sets[0].GetX()[len(sets[0].GetX())-1])
haxis.GetXaxis().SetTitle("Mass [GeV]")
haxis.GetYaxis().SetTitle("95% CL upper limit on #sigma/#sigma_{theory}")
haxis.GetYaxis().SetTitle(ylabel)
haxis.GetXaxis().SetLabelSize(0.03)
haxis.GetYaxis().SetLabelSize(0.03)
haxis.SetMinimum(miny)
haxis.SetMaximum(maxy)
haxis.GetYaxis().SetRangeUser(miny,maxy)

if doRatio:
    haxis_ratio = ROOT.TH1D("haxis_ratio","",10,sets[0].GetX()[0]-0.05*massv[0],sets[0].GetX()[len(sets[0].GetX())-1]+0.05*sets[0].GetX()[len(sets[0].GetX())-1])
    haxis_ratio.GetXaxis().SetTitle("")
    haxis_ratio.GetYaxis().SetTitle("Ratio")
    haxis_ratio.GetXaxis().SetLabelSize(0.08)
    haxis_ratio.GetYaxis().SetLabelSize(0.08)
    haxis_ratio.GetYaxis().SetRangeUser(0.0,8.0)

line = ROOT.TLine(massv[0],1.0,massv[len(massv)-1],1.0)
line.SetLineColor(ROOT.kGray)
line.SetLineStyle(2)

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetFrameLineWidth(2)
pads = []
if doRatio:
    can = ROOT.TCanvas("can","",700,700)
    pads.append(ROOT.TPad("1","1",0.0,0.3,1.0,1.0))
    pads.append(ROOT.TPad("2","2",0.0,0.01,1.0,0.3))
else:
    can = ROOT.TCanvas("can","",700,600)
    pads.append(ROOT.TPad("1","1",0,0,1,1))

#can.cd()
pads[0].SetTicky()
pads[0].SetLogy()
pads[0].SetLogx()
pads[0].Draw()
if doRatio:
    pads[1].Draw()
    pads[1].SetLogx()
#ROOT.gPad.SetTicky()
#ROOT.gPad.SetLogy()
#ROOT.gPad.SetLogx()

if doRatio:
    pads[1].cd()
    haxis_ratio.Draw()
pads[0].cd()
haxis.Draw()


for l,limit in enumerate(sets):
    pads[0].cd()
    limit.SetLineColor(colors[l])
    limit.SetLineStyle(7)
    limit.Draw("L")
if compare:
    cgexp.Draw("L")


if maskSMResonances:

    brho = ROOT.TGraph(4)
    brho.SetPoint(0, 0.7, miny+0.02*miny)
    brho.SetPoint(1, 0.89, miny+0.02*miny)
    brho.SetPoint(2, 0.89, maxy-0.02*maxy)
    brho.SetPoint(3, 0.7, maxy-0.02*maxy)
    brho.SetFillColorAlpha(ROOT.kGray, 1.0)
    brho.Draw("f")

    bphi = ROOT.TGraph(4)
    bphi.SetPoint(0, 0.91, miny+0.02*miny)
    bphi.SetPoint(1, 1.2, miny+0.02*miny)
    bphi.SetPoint(2, 1.2, maxy-0.02*maxy)
    bphi.SetPoint(3, 0.91, maxy-0.02*maxy)
    bphi.SetFillColorAlpha(ROOT.kGray, 1.0)
    bphi.Draw("f")

    bjpsi = ROOT.TGraph(4)
    bjpsi.SetPoint(0, 2.8, miny+0.02*miny)
    bjpsi.SetPoint(1, 4.1, miny+0.02*miny)
    bjpsi.SetPoint(2, 4.1, maxy-0.02*maxy)
    bjpsi.SetPoint(3, 2.8, maxy-0.02*maxy)
    bjpsi.SetFillColorAlpha(ROOT.kGray, 1.0)
    bjpsi.Draw("f, same")

    bups = ROOT.TGraph(4)
    bups.SetPoint(0, 8.6, miny+0.02*miny)
    bups.SetPoint(1, 11.0, miny+0.02*miny)
    bups.SetPoint(2, 11.0, maxy-0.02*maxy)
    bups.SetPoint(3, 8.6, maxy-0.02*maxy)
    bups.SetFillColorAlpha(ROOT.kGray, 1.0)
    bups.Draw("f, same")


#line.Draw("same")

leg.Draw("same")
pads[1].cd()
leg2.Draw("same")

if doRatio:
    pads[1].cd()
    ratio.Draw("L")
    ratio0.Draw("L")
    ratio1.Draw("L")
           
ROOT.gPad.RedrawAxis()
can.Update()
if scaleToFullLumi:
    can.SaveAs("compare_%s_ctau%s_%s_scaled.png"%(model,ctau,typeOfLimit))
else:
    can.SaveAs("compare_%s_ctau%s_%s.png"%(model,ctau,typeOfLimit))
