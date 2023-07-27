import ROOT
import os,sys,json
from datetime import date    
import numpy as np
from DataFormats.FWLite import Events, Handle

MUON_MASS = 0.10566
user = os.environ.get("USER")
today= date.today().strftime("%b-%d-%Y")

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--inDir", default="/ceph/cms/store/user/"+user+"/Run3ScoutingOutput/looperOutput_"+today, help="Choose input directory. Default: '/ceph/cms/store/user/"+user+"/Run3ScoutingOutput/looperOutput_"+today+"'")
parser.add_argument("--inSample", default="*", help="Choose sample; for all samples in input directory, choose '*'")
parser.add_argument("--inFile", default="*", help="Choose input file by index (for debug); for all files in input directory, choose '*'")
parser.add_argument("--outDir", default=os.environ.get("PWD")+"/outputGENHistograms_"+today, help="Choose output directory. Default: '"+os.environ.get("PWD")+"/outputGENHistograms_"+today+"'")
parser.add_argument("--condor", default=False, action="store_true", help="Run on condor")
parser.add_argument("--data", default=False, action="store_true", help="Process data")
parser.add_argument("--signal", default=False, action="store_true", help="Process signal")
parser.add_argument("--year", default="2022", help="Year to be processes. Default: 2022")
parser.add_argument("--partialUnblinding", default=False, action="store_true", help="Process x% (default: x=50) of available data")
parser.add_argument("--partialUnblindingFraction", default="0.5", help="Fraction of available data to be processed")
parser.add_argument("--removeDuplicates", default=False, action="store_true", help="Check for and remove duplicates")
parser.add_argument("--splitIndex", default="-1", help="Split index")
parser.add_argument("--splitPace", default="250000", help="Split pace")
parser.add_argument("--dimuonMassSel", default=[], nargs="+", help="Selection on dimuon mass: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonPtSel", default=[], nargs="+", help="Selection on dimuon pT: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--fourmuonMassSel", default=[], nargs="+", help="Selection on four-muon mass: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--fourmuonPtSel", default=[], nargs="+", help="Selection on four-muon pT: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--lxySel", default=[], nargs="+", help="Selection on lxy: first (or only) value is lower cut, second (optional) value is upper cut")
args = parser.parse_args()

indir  = args.inDir.replace("/ceph/cms","")
outdir = args.outDir

if not os.path.exists(outdir):
    os.makedirs(outdir)

isData = args.data
if not args.signal:
    isData = True

if isData:
    print("Can NOT do GEN-level analysis in data. Exiting.")
    exit()

files = []
prependtodir = ""
if not args.condor:
    prependtodir = "/ceph/cms"
else:
    prependtodir = "davs://redirector.t2.ucsd.edu:1095"
if args.inFile!="*" and args.inSample!="*":
    thisfile="output_%s_%s_%s.root"%(args.inSample,args.year,args.inFile)
    if os.path.isfile("/ceph/cms%s/%s"%(indir,thisfile)):
        files.append("%s%s/%s"%(prependtodir,indir,thisfile))
elif args.inSample!="*":
    for f in os.listdir("/ceph/cms%s"%indir):
        if (args.inSample in f) and (args.year in f) and (".root" in f) and os.path.isfile("/ceph/cms%s/%s"%(indir,f)):
            files.append("%s%s/%s"%(prependtodir,indir,f))
else:
    for f in os.listdir("/ceph/cms%s"%indir):
        if (args.year in f) and (".root" in f) and os.path.isfile("/ceph/cms%s/%s"%(indir,f)):
            files.append("%s%s/%s"%(prependtodir,indir,f))

index = int(args.splitIndex)
pace  = int(args.splitPace)

t = ROOT.TChain("tout")
for f in files:
    t.Add(f)

# Histograms:
h1d = []
h2d = []

h_nGenMuon = ROOT.TH1D("h_nGenMuon","",20,0,20)
h_nGenMuon.GetXaxis().SetTitle("Number of GEN #mu")
h_nGenMuon.GetYaxis().SetTitle("Events")
h1d.append(h_nGenMuon)

h_nGenMuonFromDP = ROOT.TH1D("h_nGenMuonFromDP","",20,0,20)
h_nGenMuonFromDP.GetXaxis().SetTitle("Number of GEN #mu from A'")
h_nGenMuonFromDP.GetYaxis().SetTitle("Events")
h1d.append(h_nGenMuonFromDP)

h_nGenDP = ROOT.TH1D("h_nGenDP","",20,0,20)
h_nGenDP.GetXaxis().SetTitle("Number of GEN A'")
h_nGenDP.GetYaxis().SetTitle("Events")
h1d.append(h_nGenDP)

h_gendimuonmass = ROOT.TH1D("h_gendimuonmass","",50,0,5)
h_gendimuonmass.GetXaxis().SetTitle("GEN m_{#mu#mu} (from A') [GeV]")
h_gendimuonmass.GetYaxis().SetTitle("Events / 0.1 GeV")
h1d.append(h_gendimuonmass)

h_gendimuonpt = ROOT.TH1D("h_gendimuonpt","",20,0,100)
h_gendimuonpt.GetXaxis().SetTitle("GEN p_{T}^{#mu#mu} (from A') [GeV]")
h_gendimuonpt.GetYaxis().SetTitle("Events / 5 GeV")
h1d.append(h_gendimuonpt)

h_gendimuonlxy = ROOT.TH1D("h_gendimuolxy","",500,0,50)
h_gendimuonlxy.GetXaxis().SetTitle("GEN l_{xy} [cm]")
h_gendimuonlxy.GetYaxis().SetTitle("Events / 0.1 cm")
h1d.append(h_gendimuonlxy)

h_genmuonptleading = ROOT.TH1D("h_genmuonptleading","",20,0,100)
h_genmuonptleading.GetXaxis().SetTitle("GEN p_{T}^{leading #mu} (from A') [GeV]")
h_genmuonptleading.GetYaxis().SetTitle("Events / 5 GeV")
h1d.append(h_genmuonptleading)

h_genmuonpttrailing = ROOT.TH1D("h_genmuonpttrailing","",20,0,100)
h_genmuonpttrailing.GetXaxis().SetTitle("GEN p_{T}^{trailing #mu} (from A') [GeV]")
h_genmuonpttrailing.GetYaxis().SetTitle("Events / 5 GeV")
h1d.append(h_genmuonpttrailing)

h_genmuonetaleading = ROOT.TH1D("h_genmuonetaleading","",50,-2.5,2.5)
h_genmuonetaleading.GetXaxis().SetTitle("GEN #eta^{leading #mu} (from A')")
h_genmuonetaleading.GetYaxis().SetTitle("Events / 0.1")
h1d.append(h_genmuonetaleading)

h_genmuonetatrailing = ROOT.TH1D("h_genmuonetatrailing","",50,-2.5,2.5)
h_genmuonetatrailing.GetXaxis().SetTitle("GEN #eta^{trailing #mu} (from A')")
h_genmuonetatrailing.GetYaxis().SetTitle("Events / 0.1")
h1d.append(h_genmuonetatrailing)

###

for h in h1d:
    h.Sumw2(ROOT.kFALSE)
for h in h2d:
    h.Sumw2(ROOT.kFALSE)

###

mid  = int(13)
dpid = int(999999)

elist = []
print("Starting loop over %d events"%t.GetEntries())
for e in range(t.GetEntries()):
    t.GetEntry(e)

    nGP  = len(t.GenPart_index)
    nGMu = 0
    nGMuFromDP = 0
    nGDP = 0
    Glxy = []
    Gmass = []
    Gptmm = []
    Gmmid = []
    Gptml = []
    Gptmt = []
    Getaml = []
    Getamt = []
    for p in range(nGP):
        if abs(t.GenPart_pdgId[p])!=mid and abs(t.GenPart_pdgId[p])!=dpid:
            continue
        if abs(t.GenPart_pdgId[p])==mid:
            nGMu = nGMu+1
            if abs(t.GenPart_motherPdgId[p])==dpid:
                nGMuFromDP = nGMuFromDP+1
                for pp in range(nGP):
                    if pp==p:
                        continue
                    if abs(t.GenPart_pdgId[pp])!=mid:
                        continue
                    if (t.GenPart_index[p] in Gmmid) or (t.GenPart_index[pp] in Gmmid):
                        continue
                    if (t.GenPart_pdgId[pp]+t.GenPart_pdgId[p])==0 and abs(t.GenPart_motherPdgId[pp])==dpid and t.GenPart_motherIndex[pp]==t.GenPart_motherIndex[p]:
                        Gmmid.append(t.GenPart_index[p])
                        Gmmid.append(t.GenPart_index[pp])
                        Glxy.append(t.GenPart_lxy[p])
                        tvec0 = ROOT.TLorentzVector()
                        tvec0.SetPtEtaPhiM(t.GenPart_pt[p], t.GenPart_eta[p], t.GenPart_phi[p], t.GenPart_m[p])
                        tvec1 = ROOT.TLorentzVector()
                        tvec1.SetPtEtaPhiM(t.GenPart_pt[pp], t.GenPart_eta[pp], t.GenPart_phi[pp], t.GenPart_m[pp])
                        Gmass.append((tvec0+tvec1).M())
                        Gptmm.append((tvec0+tvec1).Pt())
                        if tvec0.Pt()>tvec1.Pt():
                            Gptml.append(tvec0.Pt())
                            Gptmt.append(tvec1.Pt())
                            Getaml.append(tvec0.Eta())
                            Getamt.append(tvec1.Eta())
                        else:
                            Gptml.append(tvec1.Pt())
                            Gptmt.append(tvec0.Pt())
                            Getaml.append(tvec1.Eta())
                            Getamt.append(tvec0.Eta())
        elif abs(t.GenPart_pdgId[p])==dpid:
            nGDP = nGDP+1

    h_nGenDP.Fill(nGDP)
    h_nGenMuon.Fill(nGMu)
    h_nGenMuonFromDP.Fill(nGMuFromDP)
    for i in range(len(Glxy)):
        h_gendimuonmass.Fill(Gmass[i])
        h_gendimuonpt.Fill(Gptmm[i])
        h_gendimuonlxy.Fill(Glxy[i])
        h_genmuonptleading.Fill(Gptml[i])
        h_genmuonpttrailing.Fill(Gptmt[i])
        h_genmuonetaleading.Fill(Getaml[i])
        h_genmuonetatrailing.Fill(Getamt[i])

### Write histograms
foname = "%s/histograms_GEN_%s_all.root"%(outdir,args.year)
if args.inSample!="*":
    if args.inFile!="*":
        foname = "%s/histograms_GEN_file%s_%s_%s"%(outdir,args.inFile,args.inSample,args.year)
    else:
        foname = "%s/histograms_GEN_%s_%s"%(outdir,args.inSample,args.year)
if index>=0:
    foname = foname+("_%d"%index)
fout = ROOT.TFile(foname+".root","RECREATE")
fout.cd()
for h in h1d:
    h.Write()
for h in h2d:
    h.Write()
fout.Close()
