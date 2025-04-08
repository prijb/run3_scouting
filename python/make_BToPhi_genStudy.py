import ROOT
import os,sys,json
from datetime import date    
import numpy as np
import argparse
from tqdm import tqdm
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
parser.add_argument("--lxySel", default=[], nargs="+", help="Selection on lxy: first (or only) value is lower cut, second (optional) value is upper cut")
args = parser.parse_args()

#indir  = args.inDir.replace("/ceph/cms","")
#outdir = args.outDir
#
#if not os.path.exists(outdir):
#    os.makedirs(outdir)

#files = []
#prependtodir = ""
#if not args.condor:
#    prependtodir = "/ceph/cms"
#else:
#    prependtodir = "davs://redirector.t2.ucsd.edu:1095"
#if args.inFile!="*" and args.inSample!="*":
#    thisfile="output_%s_%s_%s.root"%(args.inSample,args.year,args.inFile)
#    if os.path.isfile("/ceph/cms%s/%s"%(indir,thisfile)):
#        files.append("%s%s/%s"%(prependtodir,indir,thisfile))
#elif args.inSample!="*":
#    for f in os.listdir("/ceph/cms%s"%indir):
#        if (args.inSample in f) and (args.year in f) and (".root" in f) and os.path.isfile("/ceph/cms%s/%s"%(indir,f)):
#            files.append("%s%s/%s"%(prependtodir,indir,f))
#else:
#    for f in os.listdir("/ceph/cms%s"%indir):
#        if (args.year in f) and (".root" in f) and os.path.isfile("/ceph/cms%s/%s"%(indir,f)):
#            files.append("%s%s/%s"%(prependtodir,indir,f))

index = int(args.splitIndex)
pace  = int(args.splitPace)

# Load the tree
indir = '/ceph/cms/store/group/Run3Scouting/GENScouting_noFilter_v1/BToPhi_MPhi-2p0_ctau-1mm-pythia8/private-Run3Summer22_noFilter_/250207_052508/0000/'
files = os.listdir(indir)
print(f'Number of files: {len(files)}')
#['/ceph/cms/store/group/Run3Scouting/GENScouting_noFilter_v1/BToPhi_MPhi-2p0_ctau-1mm-pythia8/private-Run3Summer22_noFilter_/250207_052508/0000/output_gensim_773.root']
events = Events(['%s/%s'%(indir,f) for f in files])

# Handler
ghandle = Handle('std::vector<reco::GenParticle>')
glabel  = ("genParticles")


# Histograms:
h1d = []

h_nbhadron = ROOT.TH1D("h_nbhadron","",10,0,10)
h_nbhadron.GetXaxis().SetTitle("Number of b-hadrons")
h_nbhadron.GetYaxis().SetTitle("Events")
h1d.append(h_nbhadron)

h_bhadron_pt = ROOT.TH1D("h_bhadron_pt","",20,0,20)
h_bhadron_pt.GetXaxis().SetTitle("b-hadron pt")
h_bhadron_pt.GetYaxis().SetTitle("Events")
h1d.append(h_bhadron_pt)

###
#
for h in h1d:
    h.Sumw2(ROOT.kFALSE)

###
#
print("Starting loop over %d events"%events.size())
for en,e in tqdm(enumerate(events), total=events.size(), desc="Processing", unit="event"):
    #
    e.getByLabel(glabel,ghandle)
    genParticles = ghandle.product()
    #
    bhadrons = []
    for gn,g in enumerate(genParticles):
        # get b-hadron
        if not (g.isLastCopy()):
            continue
        if abs(g.pdgId()) not in [521, 511, 531, 541, 5122]:
            continue
        # get b-hadron decaying to phi
        hasPhi = False
        for dn in range(g.numberOfDaughters()):
            if (abs(g.daughter(dn).pdgId())==6000211):
                hasPhi = True
        if hasPhi:
            bhadrons.append(g) 
        #else:
        #    print('This b-hadron doesnt have a phi')
    nbhadron = len(bhadrons)
    #
    if nbhadron == 1:
        h_bhadron_pt.Fill(bhadrons[0].pt())

#### Write histograms
#foname = "%s/histograms_GEN_%s_all.root"%(outdir,args.year)
#if args.inSample!="*":
#    if args.inFile!="*":
#        foname = "%s/histograms_GEN_file%s_%s_%s"%(outdir,args.inFile,args.inSample,args.year)
#    else:
#        foname = "%s/histograms_GEN_%s_%s"%(outdir,args.inSample,args.year)
#if index>=0:
#    foname = foname+("_%d"%index)
#fout = ROOT.TFile(foname+".root","RECREATE")
#fout.cd()
#for h in h1d:
#    h.Write()
#for h in h2d:
#    h.Write()
#fout.Close()
