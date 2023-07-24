import ROOT
import os,sys,json
from datetime import date    
import numpy as np
from DataFormats.FWLite import Events, Handle
sys.path.append('utils')
import plotUtils

# Functions
def applyMuonSelection(muVec):
    selected = True
    # Add your cut below
    #selected = (v.M() < 5.0)
    return selected

isData = True
removeDuplicates = False
if not isData:
    removeDuplicates = False
oncondor = False
MUON_MASS = 0.10566

if "Signal" in sys.argv[2]:
    isData = False

user = os.environ.get("USER")
today= date.today().strftime("%b-%d-%Y")
outdir = 'plots_%s'%today
if not os.path.exists(outdir):
    os.makedirs(outdir)
#os.system('cp '+os.environ.get("PWD")+'/utils/index.php '+outdir)

files = []
if not oncondor:
    for f in os.listdir("/ceph/cms%s"%sys.argv[1]):
        if len(sys.argv)>=3:
            thisfile="_%s.root"%sys.argv[2]
            if "output" in f and ".root" in f and thisfile in f and os.path.isfile("/ceph/cms%s/%s"%(sys.argv[1],f)):
                files.append("/ceph/cms%s/%s"%(sys.argv[1],f))
        else:
            if "output" in f and ".root" in f and os.path.isfile("/ceph/cms%s/%s"%(sys.argv[1],f)):
                files.append("/ceph/cms%s/%s"%(sys.argv[1],f))
else:
    infile = open('infiles.txt',"r")
    for l in infile.readlines():
        files.append("davs://redirector.t2.ucsd.edu:1095%s/%s"%(sys.argv[1],l.replace("\n","")))

index = -1
pace = 250000
if len(sys.argv)>=4:
    index = int(sys.argv[3])

t = ROOT.TChain("tout")
for f in files:
    t.Add(f)

# Histograms:
h1d = []
h2d = []

# Displaced vertices:
h_nsv = ROOT.TH1D("h_nsv","",10,0,10)
h_nsv.GetXaxis().SetTitle("Number of valid SVs")
h_nsv.GetYaxis().SetTitle("Events")
h_nsv.SetLineColor(1)
h1d.append(h_nsv)

h_nsvsel = ROOT.TH1D("h_nsvsel","",10,0,10)
h_nsvsel.GetXaxis().SetTitle("Number of selected valid SVs")
h_nsvsel.GetYaxis().SetTitle("Events")
h_nsvsel.SetLineColor(4)
h1d.append(h_nsvsel)

hsv_chi2ndof = ROOT.TH1D("hsv_chi2ndof","",100,0.0,10.0)
hsv_chi2ndof.GetXaxis().SetTitle("SV #chi^{2}/ndof")
hsv_chi2ndof.GetYaxis().SetTitle("Events / 0.1")
hsv_chi2ndof.SetLineColor(1)
h1d.append(hsv_chi2ndof)

hsv_chi2prob = ROOT.TH1D("hsv_chi2prob","",100,0.0,1.0)
hsv_chi2prob.GetXaxis().SetTitle("SV #chi^{2} probability")
hsv_chi2prob.GetYaxis().SetTitle("Events / 0.01")
hsv_chi2prob.SetLineColor(1)
h1d.append(hsv_chi2prob)

hsv_xerr = ROOT.TH1D("hsv_xerr","",50,0,0.5)
hsv_xerr.GetXaxis().SetTitle("SV x error [cm]")
hsv_xerr.GetYaxis().SetTitle("Events / 0.01 cm")
hsv_xerr.SetLineColor(1)
h1d.append(hsv_xerr)

hsv_yerr = ROOT.TH1D("hsv_yerr","",50,0,0.5)
hsv_yerr.GetXaxis().SetTitle("SV y error [cm]")
hsv_yerr.GetYaxis().SetTitle("Events / 0.01 cm")
hsv_yerr.SetLineColor(1)
h1d.append(hsv_yerr)

hsv_zerr = ROOT.TH1D("hsv_zerr","",50,0,0.5)
hsv_zerr.GetXaxis().SetTitle("SV z error [cm]")
hsv_zerr.GetYaxis().SetTitle("Events / 0.01 cm")
hsv_zerr.SetLineColor(1)
h1d.append(hsv_zerr)

hsv_lxy = ROOT.TH1D("hsv_lxy","",500,0,50)
hsv_lxy.GetXaxis().SetTitle("l_{xy} (from PV) [cm]")
hsv_lxy.GetYaxis().SetTitle("Events / 0.1 cm")
hsv_lxy.SetLineColor(1)
h1d.append(hsv_lxy)

hsvsel_lxy = ROOT.TH1D("hsvsel_lxy","",500,0,50)
hsvsel_lxy.GetXaxis().SetTitle("l_{xy} (from PV) [cm]")
hsvsel_lxy.GetYaxis().SetTitle("Events / 0.1 cm")
hsvsel_lxy.SetLineColor(4)
h1d.append(hsvsel_lxy)

hsv_l3d = ROOT.TH1D("hsv_l3d","",1000,0,100)
hsv_l3d.GetXaxis().SetTitle("l_{3D} (from PV) [cm]")
hsv_l3d.GetYaxis().SetTitle("Events / 0.1 cm")
hsv_l3d.SetLineColor(1)
h1d.append(hsv_l3d)

hsvsel_l3d = ROOT.TH1D("hsvsel_l3d","",1000,0,100)
hsvsel_l3d.GetXaxis().SetTitle("l_{3D} (from PV) [cm]")
hsvsel_l3d.GetYaxis().SetTitle("Events / 0.1 cm")
hsvsel_l3d.SetLineColor(4)
h1d.append(hsvsel_l3d)

hsvsel_mindx = ROOT.TH1D("hsvsel_mindx","",100,0,1)
hsvsel_mindx.GetXaxis().SetTitle("min D_{x} (SV_{i}, SV_{j}) [cm]")
hsvsel_mindx.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_mindx.SetLineColor(4)
h1d.append(hsvsel_mindx)

hsvsel_mindy = ROOT.TH1D("hsvsel_mindy","",100,0,1)
hsvsel_mindy.GetXaxis().SetTitle("min D_{y} (SV_{i}, SV_{j}) [cm]")
hsvsel_mindy.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_mindy.SetLineColor(4)
h1d.append(hsvsel_mindy)

hsvsel_mindz = ROOT.TH1D("hsvsel_mindz","",100,0,1)
hsvsel_mindz.GetXaxis().SetTitle("min d_{z} (SV_{i}, SV_{j}) [cm]")
hsvsel_mindz.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_mindz.SetLineColor(4)
h1d.append(hsvsel_mindz)

hsvsel_maxdx = ROOT.TH1D("hsvsel_maxdx","",100,0,1)
hsvsel_maxdx.GetXaxis().SetTitle("MAX D_{x} (SV_{i}, SV_{j}) [cm]")
hsvsel_maxdx.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_maxdx.SetLineColor(4)
h1d.append(hsvsel_maxdx)

hsvsel_maxdy = ROOT.TH1D("hsvsel_maxdy","",100,0,1)
hsvsel_maxdy.GetXaxis().SetTitle("MAX D_{y} (SV_{i}, SV_{j}) [cm]")
hsvsel_maxdy.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_maxdy.SetLineColor(4)
h1d.append(hsvsel_maxdy)

hsvsel_maxdz = ROOT.TH1D("hsvsel_maxdz","",100,0,1)
hsvsel_maxdz.GetXaxis().SetTitle("MAX D_{z} (SV_{i}, SV_{j}) [cm]")
hsvsel_maxdz.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_maxdz.SetLineColor(4)
h1d.append(hsvsel_maxdz)

# Muons
h_nmuonsass = ROOT.TH1D("h_nmuonsass","",10,0,10)
h_nmuonsass.GetXaxis().SetTitle("Number of muons (from any SV)")
h_nmuonsass.GetYaxis().SetTitle("Events")
h_nmuonsass.SetLineColor(4)
h1d.append(h_nmuonsass)

h_nmuonsassoverlap = ROOT.TH1D("h_nmuonsassoverlap","",10,0,10)
h_nmuonsassoverlap.GetXaxis().SetTitle("Number of muons (from overlapping SVs)")
h_nmuonsassoverlap.GetYaxis().SetTitle("Events")
h_nmuonsassoverlap.SetLineColor(4)
h1d.append(h_nmuonsassoverlap)

h_nmuonssel = ROOT.TH1D("h_nmuonssel","",10,0,10)
h_nmuonssel.GetXaxis().SetTitle("Number of muons (p_{T}>3 GeV, |#eta|<2.4, from any SV)")
h_nmuonssel.GetYaxis().SetTitle("Events")
h_nmuonssel.SetLineColor(4)
h1d.append(h_nmuonssel)

hmuon_pt = ROOT.TH1D("hmuon_pt","",50,0.0,150.0)
hmuon_pt.GetXaxis().SetTitle("Muon p_{T} [GeV]")
hmuon_pt.GetYaxis().SetTitle("Events / 3 GeV")
hmuon_pt.SetLineColor(2)
h1d.append(hmuon_pt)

hmuon_eta = ROOT.TH1D("hmuon_eta","",48,-2.4,2.4)
hmuon_eta.GetXaxis().SetTitle("Muon #eta")
hmuon_eta.GetYaxis().SetTitle("Events / 0.1")
hmuon_eta.SetLineColor(2)
h1d.append(hmuon_eta)

hmuon_phi = ROOT.TH1D("hmuon_phi","",64,-3.2,3.2)
hmuon_phi.GetXaxis().SetTitle("Muon #phi [rad]")
hmuon_phi.GetYaxis().SetTitle("Events / 0.1 rad")
hmuon_phi.SetLineColor(2)
h1d.append(hmuon_phi)

hmuon_ch = ROOT.TH1D("hmuon_ch","",3,-1.5,1.5)
hmuon_ch.GetXaxis().SetTitle("Muon charge")
hmuon_ch.GetYaxis().SetTitle("Events")
hmuon_ch.SetLineColor(2)
h1d.append(hmuon_ch)

hmuon_mindr = ROOT.TH1D("hmuon_mindr","",100,0.0,5.0)
hmuon_mindr.GetXaxis().SetTitle("min #DeltaR(#mu_{i}, #mu_{j})")
hmuon_mindr.GetYaxis().SetTitle("Events / 0.05")
hmuon_mindr.SetLineColor(2)
h1d.append(hmuon_mindr)

hmuon_maxdr = ROOT.TH1D("hmuon_maxdr","",100,0.0,5.0)
hmuon_maxdr.GetXaxis().SetTitle("MAX #DeltaR(#mu_{i}, #mu_{j})")
hmuon_maxdr.GetYaxis().SetTitle("Events / 0.05")
hmuon_maxdr.SetLineColor(2)
h1d.append(hmuon_maxdr)

hmuon_normchi2 = ROOT.TH1D("hmuon_normchi2","",50,0.0,5.0)
hmuon_normchi2.GetXaxis().SetTitle("Muon #chi^{2}/ndof")
hmuon_normchi2.GetYaxis().SetTitle("Events")
hmuon_normchi2.SetLineColor(2)
h1d.append(hmuon_normchi2)

mtypes = ["G & T","G & !T","!G & T","SA & !G & !T","!SA & !G & !T"]
hmuon_type = ROOT.TH1D("hmuon_type","",len(mtypes),0,len(mtypes))
hmuon_type.GetXaxis().SetTitle("Muon type")
hmuon_type.GetYaxis().SetTitle("Events")
hmuon_type.SetLineColor(2)
for b in range(1,hmuon_type.GetNbinsX()+1):
    hmuon_type.GetXaxis().SetBinLabel(b,mtypes[b-1])
h1d.append(hmuon_type)

hmuon_ecaliso = ROOT.TH1D("hmuon_ecaliso","",40,0,20)
hmuon_ecaliso.GetXaxis().SetTitle("Muon ECAL isolation [GeV]")
hmuon_ecaliso.GetYaxis().SetTitle("Events / 0.5 GeV")
hmuon_ecaliso.SetLineColor(2)
h1d.append(hmuon_ecaliso)

hmuon_ecalreliso = ROOT.TH1D("hmuon_ecalreliso","",40,0,2.0)
hmuon_ecalreliso.GetXaxis().SetTitle("Muon ECAL isolation / p_{T}")
hmuon_ecalreliso.GetYaxis().SetTitle("Events / 0.05")
hmuon_ecalreliso.SetLineColor(2)
h1d.append(hmuon_ecalreliso)

hmuon_hcaliso = ROOT.TH1D("hmuon_hcaliso","",40,0,20)
hmuon_hcaliso.GetXaxis().SetTitle("Muon HCAL isolation [GeV]")
hmuon_hcaliso.GetYaxis().SetTitle("Events / 0.5 GeV")
hmuon_hcaliso.SetLineColor(2)
h1d.append(hmuon_hcaliso)

hmuon_hcalreliso = ROOT.TH1D("hmuon_hcalreliso","",40,0,2.0)
hmuon_hcalreliso.GetXaxis().SetTitle("Muon HCAL isolation / p_{T}")
hmuon_hcalreliso.GetYaxis().SetTitle("Events / 0.05")
hmuon_hcalreliso.SetLineColor(2)
h1d.append(hmuon_hcalreliso)

hmuon_trackiso = ROOT.TH1D("hmuon_trackiso","",20,0,10)
hmuon_trackiso.GetXaxis().SetTitle("Muon track isolation [GeV]")
hmuon_trackiso.GetYaxis().SetTitle("Events / 0.5 GeV")
hmuon_trackiso.SetLineColor(2)
h1d.append(hmuon_trackiso)

hmuon_trackreliso = ROOT.TH1D("hmuon_trackreliso","",20,0,1.0)
hmuon_trackreliso.GetXaxis().SetTitle("Muon track isolation / p_{T}")
hmuon_trackreliso.GetYaxis().SetTitle("Events / 0.05")
hmuon_trackreliso.SetLineColor(2)
h1d.append(hmuon_trackreliso)

hmuon_mindrjet = ROOT.TH1D("hmuon_mindrjet","",100,0,5.0)
hmuon_mindrjet.GetXaxis().SetTitle("min #DeltaR(#mu, PF jet)")
hmuon_mindrjet.GetYaxis().SetTitle("Events / 0.05")
hmuon_mindrjet.SetLineColor(2)
h1d.append(hmuon_mindrjet)

hmuon_dxy = ROOT.TH1D("hmuon_dxy","",100,0,5.0)
hmuon_dxy.GetXaxis().SetTitle("Muon |d_{xy}| [cm]")
hmuon_dxy.GetYaxis().SetTitle("Events / 0.05 cm")
hmuon_dxy.SetLineColor(2)
h1d.append(hmuon_dxy)

hmuon_dxysig = ROOT.TH1D("hmuon_dxysig","",100,0,50)
hmuon_dxysig.GetXaxis().SetTitle("Muon |d_{xy}|/#sigma_{xy}")
hmuon_dxysig.GetYaxis().SetTitle("Events / 0.5")
hmuon_dxysig.SetLineColor(2)
h1d.append(hmuon_dxysig)

hmuon_dz = ROOT.TH1D("hmuon_dz","",300,0,30)
hmuon_dz.GetXaxis().SetTitle("Muon |d_{z}| [cm]")
hmuon_dz.GetYaxis().SetTitle("Events / 0.1 cm")
hmuon_dz.SetLineColor(2)
h1d.append(hmuon_dz)

hmuon_dzsig = ROOT.TH1D("hmuon_dzsig","",100,0,50)
hmuon_dzsig.GetXaxis().SetTitle("Muon |d_{z}|/#sigma_{z}")
hmuon_dzsig.GetYaxis().SetTitle("Events / 0.5")
hmuon_dzsig.SetLineColor(2)
h1d.append(hmuon_dzsig)

hmuon_nsahits = ROOT.TH1D("hmuon_nsahits","",50,0,50)
hmuon_nsahits.GetXaxis().SetTitle("Number of valid SA muon hits")
hmuon_nsahits.GetYaxis().SetTitle("Events")
hmuon_nsahits.SetLineColor(2)
h1d.append(hmuon_nsahits)

hmuon_nsamatchedstats = ROOT.TH1D("hmuon_nsamatchedstats","",10,0,10)
hmuon_nsamatchedstats.GetXaxis().SetTitle("Number of SA muon matched stations")
hmuon_nsamatchedstats.GetYaxis().SetTitle("Events")
hmuon_nsamatchedstats.SetLineColor(2)
h1d.append(hmuon_nsamatchedstats)

hmuon_nmuhits = ROOT.TH1D("hmuon_nmuhits","",50,0,50)
hmuon_nmuhits.GetXaxis().SetTitle("Number of valid muon hits")
hmuon_nmuhits.GetYaxis().SetTitle("Events")
hmuon_nmuhits.SetLineColor(2)
h1d.append(hmuon_nmuhits)

hmuon_nmuchambs = ROOT.TH1D("hmuon_nmuchambs","",25,0,25)
hmuon_nmuchambs.GetXaxis().SetTitle("Number of muon chambers")
hmuon_nmuchambs.GetYaxis().SetTitle("Events")
hmuon_nmuchambs.SetLineColor(2)
h1d.append(hmuon_nmuchambs)

hmuon_nmuchambsCSCorDT = ROOT.TH1D("hmuon_nmuchambsCSCorDT","",20,0,20)
hmuon_nmuchambsCSCorDT.GetXaxis().SetTitle("Number of muon chambers (CSC or DT)")
hmuon_nmuchambsCSCorDT.GetYaxis().SetTitle("Events")
hmuon_nmuchambsCSCorDT.SetLineColor(2)
h1d.append(hmuon_nmuchambsCSCorDT)

hmuon_nmumatches = ROOT.TH1D("hmuon_nmumatches","",10,0,10)
hmuon_nmumatches.GetXaxis().SetTitle("Number of muon matches")
hmuon_nmumatches.GetYaxis().SetTitle("Events")
hmuon_nmumatches.SetLineColor(2)
h1d.append(hmuon_nmumatches)

hmuon_nmumatchedstats = ROOT.TH1D("hmuon_nmumatchedstats","",10,0,10)
hmuon_nmumatchedstats.GetXaxis().SetTitle("Number of muon matched stations")
hmuon_nmumatchedstats.GetYaxis().SetTitle("Events")
hmuon_nmumatchedstats.SetLineColor(2)
h1d.append(hmuon_nmumatchedstats)

hmuon_nmuexpmatchedstats = ROOT.TH1D("hmuon_nmuexpmatchedstats","",10,0,10)
hmuon_nmuexpmatchedstats.GetXaxis().SetTitle("Number of muon expected matched stations")
hmuon_nmuexpmatchedstats.GetYaxis().SetTitle("Events")
hmuon_nmuexpmatchedstats.SetLineColor(2)
h1d.append(hmuon_nmuexpmatchedstats)

hmuon_nmumatchedstatsmexp = ROOT.TH1D("hmuon_nmumatchedstatsmexp","",20,-10,10)
hmuon_nmumatchedstatsmexp.GetXaxis().SetTitle("Number of muon matched stations - expected")
hmuon_nmumatchedstatsmexp.GetYaxis().SetTitle("Events")
hmuon_nmumatchedstatsmexp.SetLineColor(2)
h1d.append(hmuon_nmumatchedstatsmexp)

hmuon_nmumatchedRPClayers = ROOT.TH1D("hmuon_nmumatchedRPClayers","",10,0,10)
hmuon_nmumatchedRPClayers.GetXaxis().SetTitle("Number of muon matched RPC layers")
hmuon_nmumatchedRPClayers.GetYaxis().SetTitle("Events")
hmuon_nmumatchedRPClayers.SetLineColor(2)
h1d.append(hmuon_nmumatchedRPClayers)

hmuon_npixelhits = ROOT.TH1D("hmuon_npixelhits","",10,0,10)
hmuon_npixelhits.GetXaxis().SetTitle("Number of pixel hits")
hmuon_npixelhits.GetYaxis().SetTitle("Events")
hmuon_npixelhits.SetLineColor(2)
h1d.append(hmuon_npixelhits)

hmuon_npixellayers = ROOT.TH1D("hmuon_npixellayers","",10,0,10)
hmuon_npixellayers.GetXaxis().SetTitle("Number of pixel layers")
hmuon_npixellayers.GetYaxis().SetTitle("Events")
hmuon_npixellayers.SetLineColor(2)
h1d.append(hmuon_npixellayers)

hmuon_nstriphits = ROOT.TH1D("hmuon_nstriphits","",30,0,30)
hmuon_nstriphits.GetXaxis().SetTitle("Number of strip hits")
hmuon_nstriphits.GetYaxis().SetTitle("Events")
hmuon_nstriphits.SetLineColor(2)
h1d.append(hmuon_nstriphits)

hmuon_ntrackerlayers = ROOT.TH1D("hmuon_ntrackerlayers","",25,0,25)
hmuon_ntrackerlayers.GetXaxis().SetTitle("Number of tracker layers")
hmuon_ntrackerlayers.GetYaxis().SetTitle("Events")
hmuon_ntrackerlayers.SetLineColor(2)
h1d.append(hmuon_ntrackerlayers)

hmuon_mindrpfc = ROOT.TH1D("hmuon_mindrpfc","",100,0,0.1)
hmuon_mindrpfc.GetXaxis().SetTitle("min #DeltaR(#mu, PF candidate)")
hmuon_mindrpfc.GetYaxis().SetTitle("Events / 0.001")
hmuon_mindrpfc.SetLineColor(2)
h1d.append(hmuon_mindrpfc)

hmuon_pfiso0p3chg = ROOT.TH1D("hmuon_pfiso0p3chg","",40,0,20)
hmuon_pfiso0p3chg.GetXaxis().SetTitle("Muon PF charged isolation [GeV]")
hmuon_pfiso0p3chg.GetYaxis().SetTitle("Events / 0.5 GeV")
hmuon_pfiso0p3chg.SetLineColor(2)
h1d.append(hmuon_pfiso0p3chg)

hmuon_pfreliso0p3chg = ROOT.TH1D("hmuon_pfreliso0p3chg","",40,0,2.0)
hmuon_pfreliso0p3chg.GetXaxis().SetTitle("Muon PF charged isolation / p_{T}")
hmuon_pfreliso0p3chg.GetYaxis().SetTitle("Events / 0.05")
hmuon_pfreliso0p3chg.SetLineColor(2)
h1d.append(hmuon_pfreliso0p3chg)

hmuon_pfiso0p3all = ROOT.TH1D("hmuon_pfiso0p3all","",40,0,20)
hmuon_pfiso0p3all.GetXaxis().SetTitle("Muon PF all isolation (#delta#Beta) [GeV]")
hmuon_pfiso0p3all.GetYaxis().SetTitle("Events / 0.5 GeV")
hmuon_pfiso0p3all.SetLineColor(2)
h1d.append(hmuon_pfiso0p3all)

hmuon_pfreliso0p3all = ROOT.TH1D("hmuon_pfreliso0p3all","",40,0,2.0)
hmuon_pfreliso0p3all.GetXaxis().SetTitle("Muon PF all isolation (#delta#Beta) / p_{T}")
hmuon_pfreliso0p3all.GetYaxis().SetTitle("Events / 0.05")
hmuon_pfreliso0p3all.SetLineColor(2)
h1d.append(hmuon_pfreliso0p3all)

#drt = [0.05,0.1,0.2,0.3,0.4,0.5,0.8,1.1,1.5,1e10]
#drtLabels = []
#for drtn in range(len(drt)-1):
#    drtLabels.append("<%.2f"%drt[drtn])
#drtLabels.append(">%.2f"%drt[len(drt)-2])
#h_nmuons_vs_dr = ROOT.TH2D("h_nmuons_vs_dr","",len(drt),0,len(drt),10,0,10)
#h_nmuons_vs_dr.GetXaxis().SetTitle("#DeltaR(#mu_{i}, #mu_{j})")
#h_nmuons_vs_dr.GetYaxis().SetTitle("Number of other muons within (exclusive) cone")
#h_nmuons_vs_dr.GetZaxis().SetTitle("Events")
#for b in range(1,h_nmuons_vs_dr.GetNbinsX()+1):
#    h_nmuons_vs_dr.GetXaxis().SetBinLabel(b,drtLabels[b-1])
#h2d.append(h_nmuons_vs_dr)

# Di-muon

hdimuon_mass = ROOT.TH1D("hdimuon_mass","",200,0.0,100.0)
hdimuon_mass.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
hdimuon_mass.GetYaxis().SetTitle("Events / 0.5 GeV")
hdimuon_mass.SetLineColor(2)
h1d.append(hdimuon_mass)

hdimuon_lxy = ROOT.TH1D("hdimuon_lxy","",500,0.0,50.0)
hdimuon_lxy.GetXaxis().SetTitle("l_{xy} (from PV) [cm]")
hdimuon_lxy.GetYaxis().SetTitle("Events / 0.1 cm")
hdimuon_lxy.SetLineColor(2)
h1d.append(hdimuon_lxy)

hdimuon_pt = ROOT.TH1D("hdimuon_pt","",100,0.0,500.0)
hdimuon_pt.GetXaxis().SetTitle("p_{T}^{#mu#mu} [GeV]")
hdimuon_pt.GetYaxis().SetTitle("Events / 5 GeV")
hdimuon_pt.SetLineColor(2)
h1d.append(hdimuon_pt)

hdimuon_dr = ROOT.TH1D("hdimuon_dr","",100,0.0,5.0)
hdimuon_dr.GetXaxis().SetTitle("#DeltaR(#mu, #mu)")
hdimuon_dr.GetYaxis().SetTitle("Events / 0.05")
hdimuon_dr.SetLineColor(2)
h1d.append(hdimuon_dr)

hdimuon_dphi = ROOT.TH1D("hdimuon_dphi","",32,0.0,3.2)
hdimuon_dphi.GetXaxis().SetTitle("|#Delta#phi(#mu, #mu)| [rad]")
hdimuon_dphi.GetYaxis().SetTitle("Events / 0.1 rad")
hdimuon_dphi.SetLineColor(2)
h1d.append(hdimuon_dphi)

hdimuon_deta = ROOT.TH1D("hdimuon_deta","",50,0.0,5.0)
hdimuon_deta.GetXaxis().SetTitle("|#Delta#eta(#mu, #mu)|")
hdimuon_deta.GetYaxis().SetTitle("Events / 0.1")
hdimuon_deta.SetLineColor(2)
h1d.append(hdimuon_deta)

hdimuon_detaOverdphi = ROOT.TH1D("hdimuon_detaOverdphi","",20,-5.0,5.0)
hdimuon_detaOverdphi.GetXaxis().SetTitle("log_{10}|#Delta#eta(#mu, #mu) / #Delta#phi(#mu, #mu)|")
hdimuon_detaOverdphi.GetYaxis().SetTitle("Events / 0.1")
hdimuon_detaOverdphi.SetLineColor(2)
h1d.append(hdimuon_detaOverdphi)

hdimuon_3dangle = ROOT.TH1D("hdimuon_3dangle","",64,0.0,6.4)
hdimuon_3dangle.GetXaxis().SetTitle("|3D angle(#mu, #mu)|")
hdimuon_3dangle.GetYaxis().SetTitle("Events / 0.1")
hdimuon_3dangle.SetLineColor(2)
h1d.append(hdimuon_3dangle)

hdimuon_dphisv = ROOT.TH1D("hdimuon_dphisv","",32,0.0,3.2)
hdimuon_dphisv.GetXaxis().SetTitle("|#Delta#phi(#vec{(#mu, #mu)}, #vec{SV})| [rad]")
hdimuon_dphisv.GetYaxis().SetTitle("Events / 0.1 rad")
hdimuon_dphisv.SetLineColor(2)
h1d.append(hdimuon_dphisv)

hdimuon_detasv = ROOT.TH1D("hdimuon_detasv","",50,0.0,5.0)
hdimuon_detasv.GetXaxis().SetTitle("|#Delta#eta(#vec{(#mu, #mu)}, #vec{SV})|")
hdimuon_detasv.GetYaxis().SetTitle("Events / 0.1")
hdimuon_detasv.SetLineColor(2)
h1d.append(hdimuon_detasv)

hdimuon_detasvodphisv = ROOT.TH1D("hdimuon_detasvodphisv","",20,-5.0,5.0)
hdimuon_detasvodphisv.GetXaxis().SetTitle("log_{10}(|#Delta#eta(#vec{(#mu, #mu)}, #vec{SV})|/|#Delta#phi(#vec{(#mu, #mu)}, #vec{SV})|)")
hdimuon_detasvodphisv.GetYaxis().SetTitle("Events / 0.5")
hdimuon_detasvodphisv.SetLineColor(2)
h1d.append(hdimuon_detasvodphisv)

hdimuon_3danglesv = ROOT.TH1D("hdimuon_3danglesv","",64,0.0,6.4)
hdimuon_3danglesv.GetXaxis().SetTitle("|3D angle(#vec{(#mu, #mu)}, #vec{SV})|")
hdimuon_3danglesv.GetYaxis().SetTitle("Events / 0.1")
hdimuon_3danglesv.SetLineColor(2)
h1d.append(hdimuon_3danglesv)

# Selected muons

hselmuon_pt = ROOT.TH1D("hselmuon_pt","",50,0.0,150.0)
hselmuon_pt.GetXaxis().SetTitle("Muon p_{T} [GeV]")
hselmuon_pt.GetYaxis().SetTitle("Events / 3 GeV")
hselmuon_pt.SetLineColor(2)
h1d.append(hselmuon_pt)

hselmuon_eta = ROOT.TH1D("hselmuon_eta","",48,-2.4,2.4)
hselmuon_eta.GetXaxis().SetTitle("Muon #eta")
hselmuon_eta.GetYaxis().SetTitle("Events / 0.1")
hselmuon_eta.SetLineColor(2)
h1d.append(hselmuon_eta)

hselmuon_phi = ROOT.TH1D("hselmuon_phi","",64,-3.2,3.2)
hselmuon_phi.GetXaxis().SetTitle("Muon #phi [rad]")
hselmuon_phi.GetYaxis().SetTitle("Events / 0.1 rad")
hselmuon_phi.SetLineColor(2)
h1d.append(hselmuon_phi)

hselmuon_normchi2 = ROOT.TH1D("hselmuon_normchi2","",50,0.0,5.0)
hselmuon_normchi2.GetXaxis().SetTitle("Muon #chi^{2}/ndof")
hselmuon_normchi2.GetYaxis().SetTitle("Events")
hselmuon_normchi2.SetLineColor(2)
h1d.append(hselmuon_normchi2)

mtypes = ["G & T","G & !T","!G & T","SA & !G & !T","!SA & !G & !T"]
hselmuon_type = ROOT.TH1D("hselmuon_type","",len(mtypes),0,len(mtypes))
hselmuon_type.GetXaxis().SetTitle("Muon type")
hselmuon_type.GetYaxis().SetTitle("Events")
hselmuon_type.SetLineColor(2)
for b in range(1,hselmuon_type.GetNbinsX()+1):
    hselmuon_type.GetXaxis().SetBinLabel(b,mtypes[b-1])
h1d.append(hselmuon_type)

hselmuon_ecaliso = ROOT.TH1D("hselmuon_ecaliso","",40,0,20)
hselmuon_ecaliso.GetXaxis().SetTitle("Muon ECAL isolation [GeV]")
hselmuon_ecaliso.GetYaxis().SetTitle("Events / 0.5 GeV")
hselmuon_ecaliso.SetLineColor(2)
h1d.append(hselmuon_ecaliso)

hselmuon_ecalreliso = ROOT.TH1D("hselmuon_ecalreliso","",40,0,2.0)
hselmuon_ecalreliso.GetXaxis().SetTitle("Muon ECAL isolation / p_{T}")
hselmuon_ecalreliso.GetYaxis().SetTitle("Events / 0.05")
hselmuon_ecalreliso.SetLineColor(2)
h1d.append(hselmuon_ecalreliso)

hselmuon_hcaliso = ROOT.TH1D("hselmuon_hcaliso","",40,0,20)
hselmuon_hcaliso.GetXaxis().SetTitle("Muon HCAL isolation [GeV]")
hselmuon_hcaliso.GetYaxis().SetTitle("Events / 0.5 GeV")
hselmuon_hcaliso.SetLineColor(2)
h1d.append(hselmuon_hcaliso)

hselmuon_hcalreliso = ROOT.TH1D("hselmuon_hcalreliso","",40,0,2.0)
hselmuon_hcalreliso.GetXaxis().SetTitle("Muon HCAL isolation / p_{T}")
hselmuon_hcalreliso.GetYaxis().SetTitle("Events / 0.05")
hselmuon_hcalreliso.SetLineColor(2)
h1d.append(hselmuon_hcalreliso)

hselmuon_trackiso = ROOT.TH1D("hselmuon_trackiso","",20,0,10)
hselmuon_trackiso.GetXaxis().SetTitle("Muon track isolation [GeV]")
hselmuon_trackiso.GetYaxis().SetTitle("Events / 0.5 GeV")
hselmuon_trackiso.SetLineColor(2)
h1d.append(hselmuon_trackiso)

hselmuon_trackreliso = ROOT.TH1D("hselmuon_trackreliso","",20,0,1.0)
hselmuon_trackreliso.GetXaxis().SetTitle("Muon track isolation / p_{T}")
hselmuon_trackreliso.GetYaxis().SetTitle("Events / 0.05")
hselmuon_trackreliso.SetLineColor(2)
h1d.append(hselmuon_trackreliso)

hselmuon_mindrjet = ROOT.TH1D("hselmuon_mindrjet","",100,0,5.0)
hselmuon_mindrjet.GetXaxis().SetTitle("min #DeltaR(#mu, PF jet)")
hselmuon_mindrjet.GetYaxis().SetTitle("Events / 0.05")
hselmuon_mindrjet.SetLineColor(2)
h1d.append(hselmuon_mindrjet)

hselmuon_dxy = ROOT.TH1D("hselmuon_dxy","",100,0,5.0)
hselmuon_dxy.GetXaxis().SetTitle("Muon |d_{xy}| [cm]")
hselmuon_dxy.GetYaxis().SetTitle("Events / 0.05 cm")
hselmuon_dxy.SetLineColor(2)
h1d.append(hselmuon_dxy)

hselmuon_dxysig = ROOT.TH1D("hselmuon_dxysig","",100,0,50)
hselmuon_dxysig.GetXaxis().SetTitle("Muon |d_{xy}|/#sigma_{xy}")
hselmuon_dxysig.GetYaxis().SetTitle("Events / 0.5")
hselmuon_dxysig.SetLineColor(2)
h1d.append(hselmuon_dxysig)

hselmuon_dz = ROOT.TH1D("hselmuon_dz","",300,0,30)
hselmuon_dz.GetXaxis().SetTitle("Muon |d_{z}| [cm]")
hselmuon_dz.GetYaxis().SetTitle("Events / 0.1 cm")
hselmuon_dz.SetLineColor(2)
h1d.append(hselmuon_dz)

hselmuon_dzsig = ROOT.TH1D("hselmuon_dzsig","",100,0,50)
hselmuon_dzsig.GetXaxis().SetTitle("Muon |d_{z}|/#sigma_{z}")
hselmuon_dzsig.GetYaxis().SetTitle("Events / 0.5")
hselmuon_dzsig.SetLineColor(2)
h1d.append(hselmuon_dzsig)

hselmuon_nsahits = ROOT.TH1D("hselmuon_nsahits","",50,0,50)
hselmuon_nsahits.GetXaxis().SetTitle("Number of valid SA muon hits")
hselmuon_nsahits.GetYaxis().SetTitle("Events")
hselmuon_nsahits.SetLineColor(2)
h1d.append(hselmuon_nsahits)

hselmuon_nsamatchedstats = ROOT.TH1D("hselmuon_nsamatchedstats","",10,0,10)
hselmuon_nsamatchedstats.GetXaxis().SetTitle("Number of SA muon matched stations")
hselmuon_nsamatchedstats.GetYaxis().SetTitle("Events")
hselmuon_nsamatchedstats.SetLineColor(2)
h1d.append(hselmuon_nsamatchedstats)

hselmuon_nmuhits = ROOT.TH1D("hselmuon_nmuhits","",50,0,50)
hselmuon_nmuhits.GetXaxis().SetTitle("Number of valid muon hits")
hselmuon_nmuhits.GetYaxis().SetTitle("Events")
hselmuon_nmuhits.SetLineColor(2)
h1d.append(hselmuon_nmuhits)

hselmuon_nmuchambs = ROOT.TH1D("hselmuon_nmuchambs","",25,0,25)
hselmuon_nmuchambs.GetXaxis().SetTitle("Number of muon chambers")
hselmuon_nmuchambs.GetYaxis().SetTitle("Events")
hselmuon_nmuchambs.SetLineColor(2)
h1d.append(hselmuon_nmuchambs)

hselmuon_nmuchambsCSCorDT = ROOT.TH1D("hselmuon_nmuchambsCSCorDT","",20,0,20)
hselmuon_nmuchambsCSCorDT.GetXaxis().SetTitle("Number of muon chambers (CSC or DT)")
hselmuon_nmuchambsCSCorDT.GetYaxis().SetTitle("Events")
hselmuon_nmuchambsCSCorDT.SetLineColor(2)
h1d.append(hselmuon_nmuchambsCSCorDT)

hselmuon_nmumatches = ROOT.TH1D("hselmuon_nmumatches","",10,0,10)
hselmuon_nmumatches.GetXaxis().SetTitle("Number of muon matches")
hselmuon_nmumatches.GetYaxis().SetTitle("Events")
hselmuon_nmumatches.SetLineColor(2)
h1d.append(hselmuon_nmumatches)

hselmuon_nmumatchedstats = ROOT.TH1D("hselmuon_nmumatchedstats","",10,0,10)
hselmuon_nmumatchedstats.GetXaxis().SetTitle("Number of muon matched stations")
hselmuon_nmumatchedstats.GetYaxis().SetTitle("Events")
hselmuon_nmumatchedstats.SetLineColor(2)
h1d.append(hselmuon_nmumatchedstats)

hselmuon_nmuexpmatchedstats = ROOT.TH1D("hselmuon_nmuexpmatchedstats","",10,0,10)
hselmuon_nmuexpmatchedstats.GetXaxis().SetTitle("Number of muon expected matched stations")
hselmuon_nmuexpmatchedstats.GetYaxis().SetTitle("Events")
hselmuon_nmuexpmatchedstats.SetLineColor(2)
h1d.append(hselmuon_nmuexpmatchedstats)

hselmuon_nmumatchedstatsmexp = ROOT.TH1D("hselmuon_nmumatchedstatsmexp","",20,-10,10)
hselmuon_nmumatchedstatsmexp.GetXaxis().SetTitle("Number of muon matched stations - expected")
hselmuon_nmumatchedstatsmexp.GetYaxis().SetTitle("Events")
hselmuon_nmumatchedstatsmexp.SetLineColor(2)
h1d.append(hselmuon_nmumatchedstatsmexp)

hselmuon_nmumatchedRPClayers = ROOT.TH1D("hselmuon_nmumatchedRPClayers","",10,0,10)
hselmuon_nmumatchedRPClayers.GetXaxis().SetTitle("Number of muon matched RPC layers")
hselmuon_nmumatchedRPClayers.GetYaxis().SetTitle("Events")
hselmuon_nmumatchedRPClayers.SetLineColor(2)
h1d.append(hselmuon_nmumatchedRPClayers)

hselmuon_npixelhits = ROOT.TH1D("hselmuon_npixelhits","",10,0,10)
hselmuon_npixelhits.GetXaxis().SetTitle("Number of pixel hits")
hselmuon_npixelhits.GetYaxis().SetTitle("Events")
hselmuon_npixelhits.SetLineColor(2)
h1d.append(hselmuon_npixelhits)

hselmuon_npixellayers = ROOT.TH1D("hselmuon_npixellayers","",10,0,10)
hselmuon_npixellayers.GetXaxis().SetTitle("Number of pixel layers")
hselmuon_npixellayers.GetYaxis().SetTitle("Events")
hselmuon_npixellayers.SetLineColor(2)
h1d.append(hselmuon_npixellayers)

hselmuon_nstriphits = ROOT.TH1D("hselmuon_nstriphits","",30,0,30)
hselmuon_nstriphits.GetXaxis().SetTitle("Number of strip hits")
hselmuon_nstriphits.GetYaxis().SetTitle("Events")
hselmuon_nstriphits.SetLineColor(2)
h1d.append(hselmuon_nstriphits)

hselmuon_ntrackerlayers = ROOT.TH1D("hselmuon_ntrackerlayers","",25,0,25)
hselmuon_ntrackerlayers.GetXaxis().SetTitle("Number of tracker layers")
hselmuon_ntrackerlayers.GetYaxis().SetTitle("Events")
hselmuon_ntrackerlayers.SetLineColor(2)
h1d.append(hselmuon_ntrackerlayers)

hselmuon_mindrpfc = ROOT.TH1D("hselmuon_mindrpfc","",100,0,0.1)
hselmuon_mindrpfc.GetXaxis().SetTitle("min #DeltaR(#mu, PF candidate)")
hselmuon_mindrpfc.GetYaxis().SetTitle("Events / 0.001")
hselmuon_mindrpfc.SetLineColor(2)
h1d.append(hselmuon_mindrpfc)

hselmuon_pfiso0p3chg = ROOT.TH1D("hselmuon_pfiso0p3chg","",40,0,20)
hselmuon_pfiso0p3chg.GetXaxis().SetTitle("Muon PF charged isolation [GeV]")
hselmuon_pfiso0p3chg.GetYaxis().SetTitle("Events / 0.5 GeV")
hselmuon_pfiso0p3chg.SetLineColor(2)
h1d.append(hselmuon_pfiso0p3chg)

hselmuon_pfreliso0p3chg = ROOT.TH1D("hselmuon_pfreliso0p3chg","",40,0,2.0)
hselmuon_pfreliso0p3chg.GetXaxis().SetTitle("Muon PF charged isolation / p_{T}")
hselmuon_pfreliso0p3chg.GetYaxis().SetTitle("Events / 0.05")
hselmuon_pfreliso0p3chg.SetLineColor(2)
h1d.append(hselmuon_pfreliso0p3chg)

hselmuon_pfiso0p3all = ROOT.TH1D("hselmuon_pfiso0p3all","",40,0,20)
hselmuon_pfiso0p3all.GetXaxis().SetTitle("Muon PF all isolation (#delta#Beta) [GeV]")
hselmuon_pfiso0p3all.GetYaxis().SetTitle("Events / 0.5 GeV")
hselmuon_pfiso0p3all.SetLineColor(2)
h1d.append(hselmuon_pfiso0p3all)

hselmuon_pfreliso0p3all = ROOT.TH1D("hselmuon_pfreliso0p3all","",40,0,2.0)
hselmuon_pfreliso0p3all.GetXaxis().SetTitle("Muon PF all isolation (#delta#Beta) / p_{T}")
hselmuon_pfreliso0p3all.GetYaxis().SetTitle("Events / 0.05")
hselmuon_pfreliso0p3all.SetLineColor(2)
h1d.append(hselmuon_pfreliso0p3all)

# Four-muon

hfourmuon_mass = ROOT.TH1D("hfourmuon_mass","",200,0.0,100.0)
hfourmuon_mass.GetXaxis().SetTitle("m_{4#mu} [GeV]")
hfourmuon_mass.GetYaxis().SetTitle("Events / 0.5 GeV")
hfourmuon_mass.SetLineColor(2)
h1d.append(hfourmuon_mass)

hfourmuon_lxy = ROOT.TH1D("hfourmuon_lxy","",500,0.0,50.0)
hfourmuon_lxy.GetXaxis().SetTitle("l_{xy} (from PV) [cm]")
hfourmuon_lxy.GetYaxis().SetTitle("Events / 0.1 cm")
hfourmuon_lxy.SetLineColor(2)
h1d.append(hfourmuon_lxy)

hfourmuon_pt = ROOT.TH1D("hfourmuon_pt","",100,0.0,500.0)
hfourmuon_pt.GetXaxis().SetTitle("p_{T}^{4#mu} [GeV]")
hfourmuon_pt.GetYaxis().SetTitle("Events / 5 GeV")
hfourmuon_pt.SetLineColor(2)
h1d.append(hfourmuon_pt)

hfourmuon_mindr = ROOT.TH1D("hfourmuon_mindr","",100,0.0,5.0)
hfourmuon_mindr.GetXaxis().SetTitle("min #DeltaR(#mu, #mu)")
hfourmuon_mindr.GetYaxis().SetTitle("Events / 0.05")
hfourmuon_mindr.SetLineColor(2)
h1d.append(hfourmuon_mindr)

hfourmuon_mindphi = ROOT.TH1D("hfourmuon_mindphi","",32,0.0,3.2)
hfourmuon_mindphi.GetXaxis().SetTitle("min |#Delta#phi(#mu, #mu)| [rad]")
hfourmuon_mindphi.GetYaxis().SetTitle("Events / 0.1 rad")
hfourmuon_mindphi.SetLineColor(2)
h1d.append(hfourmuon_mindphi)

hfourmuon_mindeta = ROOT.TH1D("hfourmuon_mindeta","",50,0.0,5.0)
hfourmuon_mindeta.GetXaxis().SetTitle("min |#Delta#eta(#mu, #mu)|")
hfourmuon_mindeta.GetYaxis().SetTitle("Events / 0.1")
hfourmuon_mindeta.SetLineColor(2)
h1d.append(hfourmuon_mindeta)

hfourmuon_maxdr = ROOT.TH1D("hfourmuon_maxdr","",100,0.0,5.0)
hfourmuon_maxdr.GetXaxis().SetTitle("max #DeltaR(#mu, #mu)")
hfourmuon_maxdr.GetYaxis().SetTitle("Events / 0.05")
hfourmuon_maxdr.SetLineColor(2)
h1d.append(hfourmuon_maxdr)

hfourmuon_maxdphi = ROOT.TH1D("hfourmuon_maxdphi","",32,0.0,3.2)
hfourmuon_maxdphi.GetXaxis().SetTitle("max |#Delta#phi(#mu, #mu)| [rad]")
hfourmuon_maxdphi.GetYaxis().SetTitle("Events / 0.1 rad")
hfourmuon_maxdphi.SetLineColor(2)
h1d.append(hfourmuon_maxdphi)

hfourmuon_maxdeta = ROOT.TH1D("hfourmuon_maxdeta","",50,0.0,5.0)
hfourmuon_maxdeta.GetXaxis().SetTitle("max |#Delta#eta(#mu, #mu)|")
hfourmuon_maxdeta.GetYaxis().SetTitle("Events / 0.1")
hfourmuon_maxdeta.SetLineColor(2)
h1d.append(hfourmuon_maxdeta)

hfourmuon_mindetaOvermindphi = ROOT.TH1D("hfourmuon_mindetaOvermindphi","",20,-5.0,5.0)
hfourmuon_mindetaOvermindphi.GetXaxis().SetTitle("log_{10}( min|#Delta#eta(#mu, #mu)| / min|#Delta#phi(#mu, #mu)| )")
hfourmuon_mindetaOvermindphi.GetYaxis().SetTitle("Events / 0.1")
hfourmuon_mindetaOvermindphi.SetLineColor(2)
h1d.append(hfourmuon_mindetaOvermindphi)

hfourmuon_maxdetaOvermindphi = ROOT.TH1D("hfourmuon_maxdetaOvermindphi","",20,-5.0,5.0)
hfourmuon_maxdetaOvermindphi.GetXaxis().SetTitle("log_{10}( max|#Delta#eta(#mu, #mu)| / min|#Delta#phi(#mu, #mu)| )")
hfourmuon_maxdetaOvermindphi.GetYaxis().SetTitle("Events / 0.1")
hfourmuon_maxdetaOvermindphi.SetLineColor(2)
h1d.append(hfourmuon_maxdetaOvermindphi)

hfourmuon_mindetaOvermaxdphi = ROOT.TH1D("hfourmuon_mindetaOvermaxdphi","",20,-5.0,5.0)
hfourmuon_mindetaOvermaxdphi.GetXaxis().SetTitle("log_{10}( min|#Delta#eta(#mu, #mu)| / max|#Delta#phi(#mu, #mu)| )")
hfourmuon_mindetaOvermaxdphi.GetYaxis().SetTitle("Events / 0.1")
hfourmuon_mindetaOvermaxdphi.SetLineColor(2)
h1d.append(hfourmuon_mindetaOvermaxdphi)

hfourmuon_maxdetaOvermaxdphi = ROOT.TH1D("hfourmuon_maxdetaOvermaxdphi","",20,-5.0,5.0)
hfourmuon_maxdetaOvermaxdphi.GetXaxis().SetTitle("log_{10}( max|#Delta#eta(#mu, #mu)| / max|#Delta#phi(#mu, #mu)| )")
hfourmuon_maxdetaOvermaxdphi.GetYaxis().SetTitle("Events / 0.1")
hfourmuon_maxdetaOvermaxdphi.SetLineColor(2)
h1d.append(hfourmuon_maxdetaOvermaxdphi)

hfourmuon_min3dangle = ROOT.TH1D("hfourmuon_min3dangle","",64,0.0,6.4)
hfourmuon_min3dangle.GetXaxis().SetTitle("min |3D angle(#mu, #mu)|")
hfourmuon_min3dangle.GetYaxis().SetTitle("Events / 0.1")
hfourmuon_min3dangle.SetLineColor(2)
h1d.append(hfourmuon_min3dangle)

hfourmuon_max3dangle = ROOT.TH1D("hfourmuon_max3dangle","",64,0.0,6.4)
hfourmuon_max3dangle.GetXaxis().SetTitle("max |3D angle(#mu, #mu)|")
hfourmuon_max3dangle.GetYaxis().SetTitle("Events / 0.1")
hfourmuon_max3dangle.SetLineColor(2)
h1d.append(hfourmuon_max3dangle)

hfourmuon_dphisv = ROOT.TH1D("hfourmuon_dphisv","",32,0.0,3.2)
hfourmuon_dphisv.GetXaxis().SetTitle("|#Delta#phi(#vec{(4#mu)}, #vec{SV})| [rad]")
hfourmuon_dphisv.GetYaxis().SetTitle("Events / 0.1 rad")
hfourmuon_dphisv.SetLineColor(2)
h1d.append(hfourmuon_dphisv)

hfourmuon_detasv = ROOT.TH1D("hfourmuon_detasv","",50,0.0,5.0)
hfourmuon_detasv.GetXaxis().SetTitle("|#Delta#eta(#vec{(4#mu)}, #vec{SV})|")
hfourmuon_detasv.GetYaxis().SetTitle("Events / 0.1")
hfourmuon_detasv.SetLineColor(2)
h1d.append(hfourmuon_detasv)

hfourmuon_detasvodphisv = ROOT.TH1D("hfourmuon_detasvodphisv","",20,-5.0,5.0)
hfourmuon_detasvodphisv.GetXaxis().SetTitle("log_{10}(|#Delta#eta(#vec{(4#mu)}, #vec{SV})|/|#Delta#phi(#vec{(4#mu)}, #vec{SV})|)")
hfourmuon_detasvodphisv.GetYaxis().SetTitle("Events / 0.5")
hfourmuon_detasvodphisv.SetLineColor(2)
h1d.append(hfourmuon_detasvodphisv)

hfourmuon_3danglesv = ROOT.TH1D("hfourmuon_3danglesv","",64,0.0,6.4)
hfourmuon_3danglesv.GetXaxis().SetTitle("|3D angle(#vec{(4#mu)}, #vec{SV})|")
hfourmuon_3danglesv.GetYaxis().SetTitle("Events / 0.1")
hfourmuon_3danglesv.SetLineColor(2)
h1d.append(hfourmuon_3danglesv)

###

for h in h1d:
    h.Sumw2(ROOT.kFALSE)
for h in h2d:
    h.Sumw2(ROOT.kFALSE)

###

elist = []
print("Starting loop over %d events"%t.GetEntries())
firste = 0
laste  = t.GetEntries()
if index>=0:
    firste = index*pace
    laste  = min((index+1)*pace,t.GetEntries())
if firste >= t.GetEntries():
    exit()
print("From event %d to event %d"%(firste,laste))
for e in range(firste,laste):
    t.GetEntry(e)
    #if e>=1000:
    #    break
    #if e%1000==0:
    #    print("At entry %d"%e)
    if removeDuplicates:
        # Run, event number & lumisection
        eid  = t.evtn
        run  = t.run
        lumi = t.lumi
        if (run,lumi,eid) in elist:
            continue
        else:
            elist.append((run,lumi,eid))

    # Muons
    nmuonsass=0
    nmuonssel=0

    nSV = len(t.SV_index)
    if nSV<1:
        continue
    h_nsv.Fill(nSV)

    nSVsel = 0
    for v in range(nSV):
        hsv_chi2ndof.Fill(t.SV_chi2Ndof[v])
        hsv_chi2prob.Fill(t.SV_prob[v])
        hsv_xerr.Fill(t.SV_xe[v])
        hsv_yerr.Fill(t.SV_ye[v])
        hsv_zerr.Fill(t.SV_ze[v])
        hsv_lxy.Fill(t.SV_lxy[v])
        hsv_l3d.Fill(t.SV_l3d[v])
        if not t.SV_selected[v]:
            continue
        nSVsel = nSVsel+1
        hsvsel_lxy.Fill(t.SV_lxy[v])
        hsvsel_l3d.Fill(t.SV_l3d[v])
        hsvsel_mindx.Fill(t.SV_mindx[v])
        hsvsel_maxdx.Fill(t.SV_maxdx[v])
        hsvsel_mindy.Fill(t.SV_mindy[v])
        hsvsel_maxdy.Fill(t.SV_maxdy[v])
        #hsvsel_mindxy.Fill(t.SV_mindxy[v])
        #hsvsel_maxdxy.Fill(t.SV_maxdxy[v])
        hsvsel_mindz.Fill(t.SV_mindz[v])
        hsvsel_maxdz.Fill(t.SV_maxdz[v])
        #hsvsel_mind3d.Fill(t.SV_mind3d[v])
        #hsvsel_maxd3d.Fill(t.SV_maxd3d[v])
    h_nsvsel.Fill(nSVsel)

    nMu = len(t.Muon_selected)
    nMuAss = 0
    nMuAssOverlap = 0
    nMuAssSel = 0
    muselidxs = []
    for m in range(nMu):
        if t.Muon_bestAssocSVOverlapIdx[m]>-1:
            nMuAss = nMuAss+1
            nMuAssOverlap = nMuAssOverlap+1
        elif t.Muon_bestAssocSVIdx[m]>-1:
            nMuAss = nMuAss+1
        else:
            continue
        hmuon_pt.Fill(t.Muon_pt[m])
        hmuon_eta.Fill(t.Muon_eta[m])
        hmuon_phi.Fill(t.Muon_phi[m])
        if not t.Muon_selected[m]:
            continue
        nMuAssSel = nMuAssSel+1
        muselidxs.append(m)
        hmuon_ch.Fill(t.Muon_ch[m])
        hmuon_normchi2.Fill(t.Muon_chi2Ndof[m])
        if t.Muon_isGlobal[m] and t.Muon_isTracker[m]:
            hmuon_type.Fill(0.5)
        elif t.Muon_isGlobal[m] and not t.Muon_isTracker[m]:
            hmuon_type.Fill(1.5)
        elif not t.Muon_isGlobal[m] and t.Muon_isTracker[m]:
            hmuon_type.Fill(2.5)
        elif not t.Muon_isGlobal[m] and not t.Muon_isTracker[m] and t.Muon_isStandAlone[m]:
            hmuon_type.Fill(3.5)
        elif not t.Muon_isGlobal[m] and not t.Muon_isTracker[m] and not t.Muon_isStandAlone[m]:
            hmuon_type.Fill(4.5)
        hmuon_ecaliso.Fill(t.Muon_ecalIso[m])
        hmuon_ecalreliso.Fill(t.Muon_ecalRelIso[m])
        hmuon_hcaliso.Fill(t.Muon_hcalIso[m])
        hmuon_hcalreliso.Fill(t.Muon_hcalRelIso[m])
        hmuon_trackiso.Fill(t.Muon_trackIso[m])
        hmuon_trackreliso.Fill(t.Muon_trackRelIso[m])
        hmuon_mindrjet.Fill(t.Muon_mindrJet[m])
        hmuon_mindrpfc.Fill(t.Muon_mindrPF[m])
        hmuon_pfiso0p3chg.Fill(t.Muon_PFIsoChg[m])
        hmuon_pfreliso0p3chg.Fill(t.Muon_PFRelIsoChg[m])
        hmuon_pfiso0p3all.Fill(t.Muon_PFIsoAll[m])
        hmuon_pfreliso0p3all.Fill(t.Muon_PFRelIsoAll[m])
        hmuon_dxy.Fill(abs(t.Muon_dxy[m]))
        hmuon_dxysig.Fill(abs(t.Muon_dxysig[m]))
        hmuon_dz.Fill(abs(t.Muon_dz[m]))
        hmuon_dzsig.Fill(abs(t.Muon_dzsig[m]))
        hmuon_nsahits.Fill(t.Muon_saHits[m])
        hmuon_nsamatchedstats.Fill(t.Muon_saMatchedStats[m])
        hmuon_nmuhits.Fill(t.Muon_muHits[m])
        hmuon_nmuchambs.Fill(t.Muon_muChambs[m])
        hmuon_nmuchambsCSCorDT.Fill(t.Muon_muCSCDT[m])
        hmuon_nmumatches.Fill(t.Muon_muMatch[m])
        hmuon_nmumatchedstats.Fill(t.Muon_muMatchedStats[m])
        hmuon_nmuexpmatchedstats.Fill(t.Muon_muExpMatchedStats[m])
        hmuon_nmumatchedstatsmexp.Fill(t.Muon_muMatchedStats[m]-t.Muon_muExpMatchedStats[m])
        hmuon_nmumatchedRPClayers.Fill(t.Muon_muMatchedRPC[m])
        hmuon_npixelhits.Fill(t.Muon_pixHits[m])
        hmuon_npixellayers.Fill(t.Muon_pixLayers[m])
        hmuon_nstriphits.Fill(t.Muon_stripHits[m])
        hmuon_ntrackerlayers.Fill(t.Muon_trkLayers[m])
        hmuon_mindr.Fill(t.Muon_mindr[m])
        hmuon_maxdr.Fill(t.Muon_maxdr[m])
    h_nmuonsass.Fill(nMuAss)
    h_nmuonsassoverlap.Fill(nMuAssOverlap)
    h_nmuonssel.Fill(nMuAssSel)
    if nMuAssSel<2:
        continue

    dmuvec = []
    svvec = []
    svidx = []
    dmuidxs = []
    dmuvec_osv = []
    dmuidxs_osv = []
    osvvec = []
    osvidx = []
    qmuvec_osv = []
    qmuidxs_osv = []
    osvvec_qmu = []
    osvidx_qmu = []
    for m in muselidxs:
        chg   = t.Muon_ch[m]
        ovidx = -1
        vidx  = -1
        ovpos = -1
        vpos = -1
        if t.Muon_bestAssocSVOverlapIdx[m]>-1:
            ovidx = t.Muon_bestAssocSVOverlapIdx[m]
            ovpos = ovidx
        elif t.Muon_bestAssocSVIdx[m]>-1:
            vidx = t.Muon_bestAssocSVIdx[m]
            for v in range(len(t.SV_index)):
                if t.SV_index[v]==vidx:
                    vpos = v
                    break
        for mm in muselidxs:
            if mm==m:
                continue
            if abs(chg+t.Muon_ch[mm])>0:
                continue
            if ovidx>-1 and t.Muon_bestAssocSVOverlapIdx[mm]==ovidx:
                if not (m in dmuidxs_osv or mm in dmuidxs_osv):
                    dmuvec_osv.append(t.Muon_vec[m])
                    dmuvec_osv[len(dmuvec_osv)-1] = dmuvec_osv[len(dmuvec_osv)-1] + t.Muon_vec[mm]
                    dmuidxs_osv.append(m)
                    dmuidxs_osv.append(mm)
                    osvvec.append(ROOT.TVector3())
                    osvvec[len(osvvec)-1].SetXYZ(t.SVOverlap_x[ovpos]-t.PV_x, t.SVOverlap_y[ovpos]-t.PV_y, t.SVOverlap_z[ovpos]-t.PV_z)
                    osvidx.append(ovpos)
            elif vidx>-1 and t.Muon_bestAssocSVIdx[mm]==vidx:
                if not (m in dmuidxs or mm in dmuidxs):
                    dmuvec.append(t.Muon_vec[m])
                    dmuvec[len(dmuvec)-1] = dmuvec[len(dmuvec)-1] + t.Muon_vec[mm]
                    dmuidxs.append(m)
                    dmuidxs.append(mm)
                    svvec.append(ROOT.TVector3())
                    svvec[len(svvec)-1].SetXYZ(t.SV_x[vpos]-t.PV_x, t.SV_y[vpos]-t.PV_y, t.SV_z[vpos]-t.PV_z)
                    svidx.append(vpos)
    if len(dmuidxs_osv)>2 and float(len(dmuidxs_osv))/float(len(set(osvidx)))>2:
        for m in range(len(dmuidxs_osv)):
            if m%2>0:
                continue
            if len(qmuidxs_osv)>3:
                break
            for mm in range(m+2,len(dmuidxs_osv)):
                if mm%2>0:
                    continue
                if len(qmuidxs_osv)>3:
                    break
                if osvidx[int(m/2)]==osvidx[int(mm/2)]:
                    if (dmuidxs_osv[m] in qmuidxs_osv) or (dmuidxs_osv[mm] in qmuidxs_osv):
                        continue
                    else:
                        qmuidxs_osv.append(dmuidxs_osv[m])
                        qmuidxs_osv.append(dmuidxs_osv[m+1])
                        qmuidxs_osv.append(dmuidxs_osv[mm])
                        qmuidxs_osv.append(dmuidxs_osv[mm+1])
                        qmuvec_osv.append(dmuvec_osv[int(m/2)])
                        qmuvec_osv[len(qmuvec_osv)-1] = qmuvec_osv[len(qmuvec_osv)-1]+dmuvec_osv[int(mm/2)]
                        osvvec_qmu.append(ROOT.TVector3())
                        osvvec_qmu[len(osvvec_qmu)-1].SetXYZ(t.SVOverlap_x[osvidx[int(m/2)]]-t.PV_x, t.SVOverlap_y[osvidx[int(m/2)]]-t.PV_y, t.SVOverlap_z[osvidx[int(m/2)]]-t.PV_z)
                        osvidx_qmu.append(osvidx[int(m/2)])
    seldmuidxs = []
    for vn,v in enumerate(dmuvec):
        if not applyMuonSelection(v):
            continue
        mass = v.M()
        seldmuidxs.append(dmuidxs[int(vn*2)])
        seldmuidxs.append(dmuidxs[int(vn*2)+1])
        pt   = v.Pt()
        drmm = t.Muon_vec[dmuidxs[int(vn*2)]].DeltaR(t.Muon_vec[dmuidxs[int(vn*2)+1]])
        dpmm = abs(t.Muon_vec[dmuidxs[int(vn*2)]].DeltaPhi(t.Muon_vec[dmuidxs[int(vn*2)+1]]))
        demm = abs(t.Muon_vec[dmuidxs[int(vn*2)]].Eta()-t.Muon_vec[dmuidxs[int(vn*2)+1]].Eta())
        a3dmm = abs(t.Muon_vec[dmuidxs[int(vn*2)]].Angle(t.Muon_vec[dmuidxs[int(vn*2)+1]].Vect()))
        lxy  = t.SV_lxy[svidx[vn]]
        dphisv = abs(v.Vect().DeltaPhi(svvec[vn]))
        detasv = abs(v.Vect().Eta()-svvec[vn].Eta())
        detasvodphisv = ROOT.TMath.Log10(detasv/dphisv)
        a3dsv  = abs(v.Vect().Angle(svvec[vn]))
        hdimuon_mass.Fill(mass)
        hdimuon_pt.Fill(pt)
        hdimuon_dr.Fill(drmm)
        hdimuon_dphi.Fill(dpmm)
        hdimuon_deta.Fill(demm)
        if demm == 0.0:
            hdimuon_detaOverdphi.Fill(-999.0)
        elif dpmm == 0.0:
            hdimuon_detaOverdphi.Fill(999.0)
        else:
            hdimuon_detaOverdphi.Fill(ROOT.TMath.Log10(demm/dpmm))
        hdimuon_3dangle.Fill(a3dmm)
        hdimuon_lxy.Fill(lxy)
        hdimuon_dphisv.Fill(dphisv)
        hdimuon_detasv.Fill(detasv)
        hdimuon_detasvodphisv.Fill(detasvodphisv)
        hdimuon_3danglesv.Fill(a3dsv)

    seldmuidxs_osv = []
    for vn,v in enumerate(dmuvec_osv):
        if not applyMuonSelection(v):
            continue
        mass = v.M()
        seldmuidxs_osv.append(dmuidxs_osv[int(vn*2)])
        seldmuidxs_osv.append(dmuidxs_osv[int(vn*2)+1])
        pt   = v.Pt()
        drmm = t.Muon_vec[dmuidxs_osv[int(vn*2)]].DeltaR(t.Muon_vec[dmuidxs_osv[int(vn*2)+1]])
        dpmm = abs(t.Muon_vec[dmuidxs_osv[int(vn*2)]].DeltaPhi(t.Muon_vec[dmuidxs_osv[int(vn*2)+1]]))
        demm = abs(t.Muon_vec[dmuidxs_osv[int(vn*2)]].Eta()-t.Muon_vec[dmuidxs_osv[int(vn*2)+1]].Eta())
        a3dmm = abs(t.Muon_vec[dmuidxs_osv[int(vn*2)]].Angle(t.Muon_vec[dmuidxs_osv[int(vn*2)+1]].Vect()))
        lxy  = t.SVOverlap_lxy[osvidx[vn]]
        dphisv = abs(v.Vect().DeltaPhi(osvvec[vn]))
        detasv = abs(v.Vect().Eta()-osvvec[vn].Eta())
        detasvodphisv = ROOT.TMath.Log10(detasv/dphisv)
        a3dsv  = abs(v.Vect().Angle(osvvec[vn]))
        hdimuon_mass.Fill(mass)
        hdimuon_pt.Fill(pt)
        hdimuon_dr.Fill(drmm)
        hdimuon_dphi.Fill(dpmm)
        hdimuon_deta.Fill(demm)
        if demm == 0.0:
            hdimuon_detaOverdphi.Fill(-999.0)
        elif dpmm == 0.0:
            hdimuon_detaOverdphi.Fill(999.0)
        else:
            hdimuon_detaOverdphi.Fill(ROOT.TMath.Log10(demm/dpmm))
        hdimuon_3dangle.Fill(a3dmm)
        hdimuon_lxy.Fill(lxy)
        hdimuon_dphisv.Fill(dphisv)
        hdimuon_detasv.Fill(detasv)
        hdimuon_detasvodphisv.Fill(detasvodphisv)
        hdimuon_3danglesv.Fill(a3dsv)

    mindrmm, mindpmm, mindemm, mina3dmm = 1e6, 1e6, 1e6, 1e6
    maxdrmm, maxdpmm, maxdemm, maxa3dmm = -1., -1., -1., -1
    for vn,v in enumerate(qmuvec_osv):
        if not applyMuonSelection(v):
            continue
        mass = v.M()
        pt   = v.Pt()
        for m in range(vn*4,vn*4+4):
            for mm in range(m+1,vn*4+4):
                drmm = t.Muon_vec[m].DeltaR(t.Muon_vec[mm])
                dpmm = abs(t.Muon_vec[m].DeltaPhi(t.Muon_vec[mm]))
                demm = abs(t.Muon_vec[m].Eta()-t.Muon_vec[mm].Eta())
                a3dmm = abs(t.Muon_vec[m].Angle(t.Muon_vec[mm].Vect()))
                if drmm<mindrmm:
                    mindrmm = drmm
                if drmm>maxdrmm:
                    maxdrmm = drmm
                if dpmm<mindpmm:
                    mindpmm = dpmm
                if dpmm>maxdpmm:
                    maxdpmm = dpmm
                if demm<mindemm:
                    mindemm = demm
                if demm>maxdemm:
                    maxdemm = demm
                if a3dmm<mina3dmm:
                    mina3dmm = a3dmm
                if a3dmm>maxa3dmm:
                    maxa3dmm = a3dmm                    
        lxy  = t.SVOverlap_lxy[osvidx_qmu[vn]]
        dphisv = abs(v.Vect().DeltaPhi(osvvec_qmu[vn]))
        detasv = abs(v.Vect().Eta()-osvvec_qmu[vn].Eta())
        detasvodphisv = ROOT.TMath.Log10(detasv/dphisv)
        a3dsv  = abs(v.Vect().Angle(osvvec_qmu[vn]))
        hfourmuon_mass.Fill(mass)
        hfourmuon_pt.Fill(pt)
        hfourmuon_mindr.Fill(mindrmm)
        hfourmuon_mindphi.Fill(mindpmm)
        hfourmuon_mindeta.Fill(mindemm)
        hfourmuon_min3dangle.Fill(mina3dmm)
        hfourmuon_maxdr.Fill(maxdrmm)
        hfourmuon_maxdphi.Fill(maxdpmm)
        hfourmuon_maxdeta.Fill(maxdemm)
        if mindemm == 0.0:
            hfourmuon_mindetaOvermindphi.Fill(-999.0)
            hfourmuon_mindetaOvermaxdphi.Fill(-999.0)
        else:
            if mindpmm == 0.0:
                hfourmuon_mindetaOvermindphi.Fill(999.0)
            else:
                hfourmuon_mindetaOvermindphi.Fill(ROOT.TMath.Log10(mindemm/mindpmm))
            if maxdpmm == 0.0:
                hfourmuon_mindetaOvermaxdphi.Fill(999.0)
            else:
                hfourmuon_mindetaOvermaxdphi.Fill(ROOT.TMath.Log10(mindemm/maxdpmm))
        if maxdemm == 0.0:
            hfourmuon_maxdetaOvermindphi.Fill(-999.0)
            hfourmuon_maxdetaOvermaxdphi.Fill(-999.0)
        else:
            if mindpmm == 0.0:
                hfourmuon_maxdetaOvermindphi.Fill(999.0)
            else:
                hfourmuon_maxdetaOvermindphi.Fill(ROOT.TMath.Log10(mindemm/mindpmm))
            if maxdpmm == 0.0:
                hfourmuon_maxdetaOvermaxdphi.Fill(999.0)
            else:
                hfourmuon_maxdetaOvermaxdphi.Fill(ROOT.TMath.Log10(mindemm/maxdpmm))
        hfourmuon_max3dangle.Fill(maxa3dmm)
        hfourmuon_lxy.Fill(lxy)
        hfourmuon_dphisv.Fill(dphisv)
        hfourmuon_detasv.Fill(detasv)
        hfourmuon_detasvodphisv.Fill(detasvodphisv)
        hfourmuon_3danglesv.Fill(a3dsv)

    selmuidxs = seldmuidxs+seldmuidxs_osv
    for m in selmuidxs:
        hselmuon_pt.Fill(t.Muon_pt[m])
        hselmuon_eta.Fill(t.Muon_eta[m])
        hselmuon_phi.Fill(t.Muon_phi[m])
        hselmuon_normchi2.Fill(t.Muon_chi2Ndof[m])
        if t.Muon_isGlobal[m] and t.Muon_isTracker[m]:
            hselmuon_type.Fill(0.5)
        elif t.Muon_isGlobal[m] and not t.Muon_isTracker[m]:
            hselmuon_type.Fill(1.5)
        elif not t.Muon_isGlobal[m] and t.Muon_isTracker[m]:
            hselmuon_type.Fill(2.5)
        elif not t.Muon_isGlobal[m] and not t.Muon_isTracker[m] and t.Muon_isStandAlone[m]:
            hselmuon_type.Fill(3.5)
        elif not t.Muon_isGlobal[m] and not t.Muon_isTracker[m] and not t.Muon_isStandAlone[m]:
            hselmuon_type.Fill(4.5)
        hselmuon_ecaliso.Fill(t.Muon_ecalIso[m])
        hselmuon_ecalreliso.Fill(t.Muon_ecalRelIso[m])
        hselmuon_hcaliso.Fill(t.Muon_hcalIso[m])
        hselmuon_hcalreliso.Fill(t.Muon_hcalRelIso[m])
        hselmuon_trackiso.Fill(t.Muon_trackIso[m])
        hselmuon_trackreliso.Fill(t.Muon_trackRelIso[m])
        hselmuon_mindrjet.Fill(t.Muon_mindrJet[m])
        hselmuon_mindrpfc.Fill(t.Muon_mindrPF[m])
        hselmuon_pfiso0p3chg.Fill(t.Muon_PFIsoChg[m])
        hselmuon_pfreliso0p3chg.Fill(t.Muon_PFRelIsoChg[m])
        hselmuon_pfiso0p3all.Fill(t.Muon_PFIsoAll[m])
        hselmuon_pfreliso0p3all.Fill(t.Muon_PFRelIsoAll[m])
        hselmuon_dxy.Fill(abs(t.Muon_dxy[m]))
        hselmuon_dxysig.Fill(abs(t.Muon_dxysig[m]))
        hselmuon_dz.Fill(abs(t.Muon_dz[m]))
        hselmuon_dzsig.Fill(abs(t.Muon_dzsig[m]))
        hselmuon_nsahits.Fill(t.Muon_saHits[m])
        hselmuon_nsamatchedstats.Fill(t.Muon_saMatchedStats[m])
        hselmuon_nmuhits.Fill(t.Muon_muHits[m])
        hselmuon_nmuchambs.Fill(t.Muon_muChambs[m])
        hselmuon_nmuchambsCSCorDT.Fill(t.Muon_muCSCDT[m])
        hselmuon_nmumatches.Fill(t.Muon_muMatch[m])
        hselmuon_nmumatchedstats.Fill(t.Muon_muMatchedStats[m])
        hselmuon_nmuexpmatchedstats.Fill(t.Muon_muExpMatchedStats[m])
        hselmuon_nmumatchedstatsmexp.Fill(t.Muon_muMatchedStats[m]-t.Muon_muExpMatchedStats[m])
        hselmuon_nmumatchedRPClayers.Fill(t.Muon_muMatchedRPC[m])
        hselmuon_npixelhits.Fill(t.Muon_pixHits[m])
        hselmuon_npixellayers.Fill(t.Muon_pixLayers[m])
        hselmuon_nstriphits.Fill(t.Muon_stripHits[m])
        hselmuon_ntrackerlayers.Fill(t.Muon_trkLayers[m])

doPlots = True
if len(sys.argv)>2:
    doPlots = False
# Draw histograms
unityArea = True
doLogy = True

# Labels
yearenergy = "1.46 fb^{-1} (Run2022D, 13.6 TeV)"
cmsExtra = "Preliminary"
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
can = ROOT.TCanvas("can","",600,600)
for h in h1d:
    if not doPlots:
        break
    if "_type" not in h.GetName():
        plotUtils.PutUnderflowInFirstBin(h)
        plotUtils.PutOverflowInLastBin(h)
    h.SetBinErrorOption(ROOT.TH1.kPoisson)
    h.GetYaxis().SetLabelSize(0.025)
    h.GetYaxis().SetMaxDigits(3)
    h.GetYaxis().SetRangeUser(0.0, 1.1*h.GetMaximum())
    h.SetLineWidth(2)
    ytitle = h.GetYaxis().GetTitle()
    if doLogy:
        h.GetYaxis().SetRangeUser(0.9, 2.0*h.GetMaximum())
    if unityArea:
        ytitle = ytitle.replace("Events", "Fraction of events")
        if h.Integral(0,-1)>0:
            h.Scale(1.0/h.Integral(0,-1))
        h.GetYaxis().SetRangeUser(0.0,1.1*h.GetMaximum())
        if doLogy:
            h.GetYaxis().SetRangeUser(0.9/h.GetEntries(), 2.0*h.GetMaximum())
    if "hsv" in h.GetName():
        ytitle = ytitle.replace("Events", "Number of SVs")
        ytitle = ytitle.replace("events", "SVs")
    if "hmuon" in h.GetName():
        ytitle = ytitle.replace("Events", "Number of muons")
        ytitle = ytitle.replace("events", "muons")
    h.GetYaxis().SetTitle(ytitle)
    can.cd()
    if doLogy:
        ROOT.gPad.SetLogy()
    h.Draw("E")
    can.Update()
    latex.DrawLatex(0.90, 0.91, yearenergy);
    if drawCMSOnTop:
        latexCMS.DrawLatex(0.1,0.91,"CMS");
        latexCMSExtra.DrawLatex(0.19,0.91, cmsExtra);
    else:
        latexCMS.DrawLatex(0.13,0.86,"CMS");
        latexCMSExtra.DrawLatex(0.13,0.815, cmsExtra);
    ROOT.gPad.RedrawAxis()
    outname = "%s"%h.GetName()
    if "h_" in outname:
        outname = outname.replace("h_","")
    elif "hsv" in outname or "hmuon" in outname:
        outname = outname.replace("hsv","sv")
        outname = outname.replace("hmuon","muon")
    can.SaveAs("%s/%s.png"%(outdir,outname))
    can.Clear()

for h in h2d:
    if not doPlots:
        break
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
    outname = "%s"%h.GetName()
    if "h_" in outname:
        outname = outname.replace("h_","")
    can.SaveAs("%s/%s.png"%(outdir,outname))
    can.Clear()

del can

foname = "%s/histograms_all.root"%outdir
if len(sys.argv)>2:
    if len(sys.argv)>3:
        foname = "%s/histograms_file%s_%s.root"%(outdir,sys.argv[2],sys.argv[3])
    else:
        foname = "%s/histograms_file%s.root"%(outdir,sys.argv[2])
fout = ROOT.TFile(foname,"RECREATE")
fout.cd()
for h in h1d:
    h.Write()
for h in h2d:
    h.Write()
fout.Close()
