import ROOT
import os,sys,json
import argparse
from datetime import date    
import numpy as np
from DataFormats.FWLite import Events, Handle
sys.path.append('utils')
import histDefinition

ROOT.EnableImplicitMT(2)

user = os.environ.get("USER")
today= date.today().strftime("%b-%d-%Y")

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--inDir", default="/ceph/cms/store/user/"+user+"/Run3ScoutingOutput/looperOutput_"+today, help="Choose input directory. Default: '/ceph/cms/store/user/"+user+"/Run3ScoutingOutput/looperOutput_"+today+"'")
parser.add_argument("--inSample", default="*", help="Choose sample; for all samples in input directory, choose '*'")
parser.add_argument("--inFile", default="*", help="Choose input file by index (for debug); for all files in input directory, choose '*'")
parser.add_argument("--outDir", default=os.environ.get("PWD")+"/outputHistograms_"+today, help="Choose output directory. Default: '"+os.environ.get("PWD")+"/outputHistograms_"+today+"'")
parser.add_argument("--outSuffix", default="", help="Choose output directory. Default: ''")
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
parser.add_argument("--dimuonMassSidebandSel", default=[], nargs="+", help="Selection on dimuon mass sidebands: first pair of values is left sideband, second (optional) is right sideband")
parser.add_argument("--dimuonPtSel", default=[], nargs="+", help="Selection on dimuon pT: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--fourmuonMassSel", default=[], nargs="+", help="Selection on four-muon mass: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--fourmuonPtSel", default=[], nargs="+", help="Selection on four-muon pT: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassSelForFourMuon", default=[], nargs="+", help="Selection on dimuon mass in four-muon system: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassDiffSelForFourMuon", default=[], nargs="+", help="Selection on dimuon mass difference / mean in four-muon system: first (or only) value is *upper* cut, second (optional) value is *lower* cut")
parser.add_argument("--dimuonPtSelForFourMuon", default=[], nargs="+", help="Selection on dimuon pT in four-muon system: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassSelForFourMuonOSV", default=[], nargs="+", help="Selection on dimuon mass in four-muon system from overlapping SV: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassDiffSelForFourMuonOSV", default=[], nargs="+", help="Selection on dimuon mass difference / mean in four-muon system from overlapping SV: first (or only) value is *upper* cut, second (optional) value is *lower* cut")
parser.add_argument("--dimuonPtSelForFourMuonOSV", default=[], nargs="+", help="Selection on dimuon pT in four-muon system from overlapping SV: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--lxySel", default=[], nargs="+", help="Selection on lxy: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--lxySelForFourMuon", default=[], nargs="+", help="Selection on lxy: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--noPreSel", default=False, action="store_true", help="Do not fill pre-selection/association histograms")
parser.add_argument("--noDiMuon", default=False, action="store_true", help="Do not fill dimuon histograms")
parser.add_argument("--noFourMuon", default=False, action="store_true", help="Do not fill four-muon histograms for four-muon systems")
parser.add_argument("--noFourMuonOSV", default=False, action="store_true", help="Do not fill four-muon histograms for four-muon systems from overlapping SVs")
args = parser.parse_args()

# Functions for selection
def applyDiMuonSelection(vec):
    selected = True
    if len(args.dimuonMassSel)>0:
        selected = selected and (vec.M() > float(args.dimuonMassSel[0]))
    if len(args.dimuonMassSel)>1:
        selected = selected and (vec.M() < float(args.dimuonMassSel[1]))
    if len(args.dimuonPtSel)>0:
        selected = selected and (vec.Pt() > float(args.dimuonPtSel[0]))
    if len(args.dimuonPtSel)>1:
        selected = selected and (vec.Pt() < float(args.dimuonPtSel[1]))
    if len(args.dimuonMassSidebandSel)>3:
        selected = selected and ((vec.M() > float(args.dimuonMassSidebandSel[0]) and vec.M() < float(args.dimuonMassSidebandSel[1])) or
                                 (vec.M() > float(args.dimuonMassSidebandSel[2]) and vec.M() < float(args.dimuonMassSidebandSel[3])))
    return selected

def applyFourMuonSelection(vec):
    selected = True
    if len(args.fourmuonMassSel)>0:
        selected = selected and (vec.M() > float(args.fourmuonMassSel[0]))
    if len(args.fourmuonMassSel)>1:
        selected = selected and (vec.M() < float(args.fourmuonMassSel[1]))
    if len(args.fourmuonPtSel)>0:
        selected = selected and (vec.Pt() > float(args.fourmuonPtSel[0]))
    if len(args.fourmuonPtSel)>1:
        selected = selected and (vec.Pt() < float(args.fourmuonPtSel[1]))
    return selected

def applyDiMuonSelectionForFourMuon(vecf,vecs):
    selected = True
    if len(args.dimuonMassSelForFourMuon)>0:
        selected = selected and (vecf.M() > float(args.dimuonMassSelForFourMuon[0])) and (vecs.M() > float(args.dimuonMassSelForFourMuon[0]))
    if len(args.dimuonMassSelForFourMuon)>1:
        selected = selected and (vecf.M() < float(args.dimuonMassSelForFourMuon[1])) and (vecs.M() < float(args.dimuonMassSelForFourMuon[1]))
    if len(args.dimuonPtSelForFourMuon)>0:
        selected = selected and (vecf.Pt() > float(args.dimuonPtSelForFourMuon[0])) and (vecs.Pt() > float(args.dimuonPtSelForFourMuon[0]))
    if len(args.dimuonPtSelForFourMuon)>1:
        selected = selected and (vecf.Pt() < float(args.dimuonPtSelForFourMuon[1])) and (vecs.Pt() < float(args.dimuonPtSelForFourMuon[1]))
    if len(args.dimuonMassDiffSelForFourMuon)>0:
        selected = selected and (abs(vecf.M()-vecs.M())*2.0/(vecf.M()+vecs.M()) < float(args.dimuonMassDiffSelForFourMuon[0]))
    if len(args.dimuonMassDiffSelForFourMuon)>1:
        selected = selected and (abs(vecf.M()-vecs.M())*2.0/(vecf.M()+vecs.M()) > float(args.dimuonMassDiffSelForFourMuon[1]))
    return selected

def applyDiMuonSelectionForFourMuonOSV(vecf,vecs):
    selected = True
    if len(args.dimuonMassSelForFourMuonOSV)>0:
        selected = selected and (vecf.M() > float(args.dimuonMassSelForFourMuonOSV[0])) and (vecs.M() > float(args.dimuonMassSelForFourMuonOSV[0]))
    if len(args.dimuonMassSelForFourMuonOSV)>1:
        selected = selected and (vecf.M() < float(args.dimuonMassSelForFourMuonOSV[1])) and (vecs.M() < float(args.dimuonMassSelForFourMuonOSV[1]))
    if len(args.dimuonPtSelForFourMuonOSV)>0:
        selected = selected and (vecf.Pt() > float(args.dimuonPtSelForFourMuonOSV[0])) and (vecs.Pt() > float(args.dimuonPtSelForFourMuonOSV[0]))
    if len(args.dimuonPtSelForFourMuonOSV)>1:
        selected = selected and (vecf.Pt() < float(args.dimuonPtSelForFourMuonOSV[1])) and (vecs.Pt() < float(args.dimuonPtSelForFourMuonOSV[1]))
    if len(args.dimuonMassDiffSelForFourMuonOSV)>0:
        selected = selected and (abs(vecf.M()-vecs.M())*2.0/(vecf.M()+vecs.M()) < float(args.dimuonMassDiffSelForFourMuonOSV[0]))
    if len(args.dimuonMassDiffSelForFourMuonOSV)>1:
        selected = selected and (abs(vecf.M()-vecs.M())*2.0/(vecf.M()+vecs.M()) > float(args.dimuonMassDiffSelForFourMuonOSV[1]))
    return selected

def applyLxySelection(lxy):
    selected = True
    if len(args.lxySel)>0:
        selected = selected and (lxy > float(args.lxySel[0]))
    if len(args.lxySel)>1:
        selected = selected and (lxy < float(args.lxySel[1]))
    return selected

def applyFourMuonLxySelection(lxymin,lxymax):
    selected = True
    if len(args.lxySelForFourMuon)>0:
        selected = selected and (lxymin > float(args.lxySel[0]))
    if len(args.lxySelForFourMuon)>1:
        selected = selected and (lxymax < float(args.lxySel[1]))
    return selected

# Muon type
def muonType(isGlobal, isTracker, isStandAlone):
    if isGlobal and isTracker:
        return 0.5
    elif isGlobal and not isTracker:
        return 1.5
    elif not isGlobal and isTracker:
        return 2.5
    elif not isGlobal and not isTracker and isStandAlone:
        return 3.5
    elif not isGlobal and not isTracker and not isStandAlone:
        return 4.5

indir  = args.inDir.replace("/ceph/cms","")
outdir = args.outDir
if args.outSuffix!="":
    outdir = outdir+"_"+args.outSuffix
if not os.path.exists(outdir):
    os.makedirs(outdir)

isData = args.data
if "Data" in args.inSample:
    isData = True
removeDuplicates = args.removeDuplicates
if not isData:
    removeDuplicates = False
MUON_MASS = 0.10566

if args.signal:
    isData = False

skimFraction = float(args.partialUnblindingFraction)
skimEvents = args.partialUnblinding
if not isData:
    skimEvents = False
rndm_partialUnblinding = ROOT.TRandom3(42)

files = []
prependtodir = ""
if not args.condor:
    prependtodir = "/ceph/cms"
else:
    prependtodir = "davs://redirector.t2.ucsd.edu:1095"
if not args.condor:
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
else:
    os.system('xrdfs redirector.t2.ucsd.edu:1095 ls %s > filein.txt'%indir)
    fin = open("filein.txt","r")
    for f in fin.readlines():
        f = f.strip("\n")
        if args.inFile!="*" and args.inSample!="*":
            thisfile="output_%s_%s_%s.root"%(args.inSample,args.year,args.inFile)
            if thisfile in f:
                files.append("%s/%s"%(prependtodir,thisfile))
        elif args.inSample!="*":
            if (args.inSample in f) and (args.year in f) and (".root" in f):
                files.append("%s/%s"%(prependtodir,f))
        else:
            if (args.year in f) and (".root" in f):
                files.append("%s/%s"%(prependtodir,f))
    fin.close()
    os.system('rm -f filein.txt')

index = int(args.splitIndex)
pace  = int(args.splitPace)

t = ROOT.TChain("tout")
for f in files:
    t.Add(f)

# Histograms:
h1d = []
variable1d = dict()
#
h2d = []
variable2d = dict()
#
h1d,variable1d,h2d,variable2d = histDefinition.histInitialization(not(args.noPreSel),not(args.noDiMuon),not(args.noFourMuon),not(args.noFourMuonOSV))
###
hall = h1d+h2d
for h in hall:
    if isData:
        h.Sumw2(ROOT.kFALSE)
    else:
        h.Sumw2()
###

elist = [] # for duplicate removal
print("Starting loop over %d events"%t.GetEntries())
firste = 0
laste  = t.GetEntries()
print(args.inSample, laste)
if index>=0:
    firste = index*pace
    laste  = min((index+1)*pace,t.GetEntries())
if firste >= t.GetEntries():
    exit()
print("From event %d to event %d"%(firste,laste))
for e in range(firste,laste):
    t.GetEntry(e)
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
    if skimEvents and skimFraction>0.0:
        if rndm_partialUnblinding.Rndm() > skimFraction:
            continue
    ### As advised in LUM POG TWiki 
    ### (https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun3),
    ### exclude runs 359571 + 359661
    if isData and t.run==359571 or t.run==359661:
            continue

    # Loop over SVs
    nSV = len(t.SV_index)
    if nSV<1:
        continue
    nSVsel = 0
    for v in range(nSV):
        if args.noPreSel:
            break
        if not t.SV_selected[v]:
            continue
        nSVsel = nSVsel+1
        lxy = t.SV_lxy[v]
        for h in h1d:
            tn = h.GetName()
            if "hsvsel_" not in tn:
                continue
            else:
                h.Fill(eval(variable1d[h.GetName()]))
        for h in h2d:
            tn = h.GetName()
            if "hsvsel_" not in tn:
                continue
            else:
                h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]))
    nSVs = nSVsel
    for h in h1d:
        if args.noPreSel:
            break
        tn = h.GetName()
        if "h_nsvsel" not in tn or "ass" in tn:
            continue
        else:
            h.Fill(eval(variable1d[h.GetName()]))

    # Loop over muons
    nMu = len(t.Muon_selected)
    nMuSel = 0
    nMuAss = 0
    nMuAssOverlap = 0
    muselidxs = []
    for m in range(nMu):
        if not t.Muon_selected[m]:
            continue
        nMuSel = nMuSel+1
        if t.Muon_bestAssocSVOverlapIdx[m]>-1:
            nMuAss = nMuAss+1
            nMuAssOverlap = nMuAssOverlap+1
        elif t.Muon_bestAssocSVIdx[m]>-1:
            nMuAss = nMuAss+1
        else:
            continue
        muselidxs.append(m)
        if args.noPreSel:
            continue
        for h in h1d:
            if args.noPreSel:
                break
            tn = h.GetName()
            if "hmuon_" not in tn:
                continue
            else:
                h.Fill(eval(variable1d[h.GetName()]))
        for h in h2d:
            if args.noPreSel:
                break
            tn = h.GetName()
            if "hmuon_" not in tn:
                continue
            else:
                h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]))
    for h in h1d:
        if args.noPreSel:
            break
        tn = h.GetName()
        if "h_nmuons" not in tn:
            continue
        else:
            h.Fill(eval(variable1d[h.GetName()]))
    # Select events witb at least two muons associated to a SV
    if nMuAss<2:
        continue

    # Muon pairing
    dmuvec = []
    svvec = []
    svidx = []
    dmuidxs = []
    dmuvec_osv = []
    dmuidxs_osv = []
    osvvec = []
    osvidx = []
    qmuvec_osv = []
    qmu_dmuvec_osv = []
    qmuidxs_osv = []
    qmuidxs_osv_sel = []
    osvvec_qmu = []
    osvidx_qmu = []
    qmuvec = []
    qmuidxs = []
    qmuidxs_sel = []
    qmu_dmuvecminlxy = []
    qmu_dmuvecmaxlxy = []
    svvecminlxy_qmu = []
    svidxminlxy_qmu = []
    svvecmaxlxy_qmu = []
    svidxmaxlxy_qmu = []
    for m in muselidxs:
        chg   = t.Muon_ch[m]
        ovidx = -1
        vidx  = -1
        ovpos = -1
        vpos = -1
        # First, identify muons from overlapping SVs
        if t.Muon_bestAssocSVOverlapIdx[m]>-1:
            ovidx = t.Muon_bestAssocSVOverlapIdx[m]
            ovpos = ovidx
        # Then, identify muons from non-overlapping SVs
        elif t.Muon_bestAssocSVIdx[m]>-1:
            vidx = t.Muon_bestAssocSVIdx[m]
            for v in range(len(t.SV_index)):
                if t.SV_index[v]==vidx:
                    vpos = v
                    break
        # Loop over muons, and do pairing
        for mm in muselidxs:
            if mm==m:
                continue
            if abs(chg+t.Muon_ch[mm])>0:
                continue
            # First, identify muon pairs from overlapping SVs
            if ovidx>-1 and t.Muon_bestAssocSVOverlapIdx[mm]==ovidx:
                if not (m in dmuidxs_osv or mm in dmuidxs_osv):
                    dmuvec_osv.append(t.Muon_vec[m])
                    dmuvec_osv[len(dmuvec_osv)-1] = dmuvec_osv[len(dmuvec_osv)-1] + t.Muon_vec[mm]
                    dmuidxs_osv.append(m)
                    dmuidxs_osv.append(mm)
                    osvvec.append(ROOT.TVector3())
                    osvvec[len(osvvec)-1].SetXYZ(t.SVOverlap_x[ovpos]-t.PV_x, t.SVOverlap_y[ovpos]-t.PV_y, t.SVOverlap_z[ovpos]-t.PV_z)
                    osvidx.append(ovpos)
            # Then, identify muon pairs from non-overlapping SVs
            elif vidx>-1 and t.Muon_bestAssocSVIdx[mm]==vidx:
                if not (m in dmuidxs or mm in dmuidxs):
                    dmuvec.append(t.Muon_vec[m])
                    dmuvec[len(dmuvec)-1] = dmuvec[len(dmuvec)-1] + t.Muon_vec[mm]
                    dmuidxs.append(m)
                    dmuidxs.append(mm)
                    svvec.append(ROOT.TVector3())
                    svvec[len(svvec)-1].SetXYZ(t.SV_x[vpos]-t.PV_x, t.SV_y[vpos]-t.PV_y, t.SV_z[vpos]-t.PV_z)
                    svidx.append(vpos)

    # If multiple muon pairs from overlapping SVs are found, create four-muon system
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
                        qmu_dmuvec_osv.append(dmuvec_osv[int(m/2)])
                        qmu_dmuvec_osv.append(dmuvec_osv[int(mm/2)])
                        osvvec_qmu.append(ROOT.TVector3())
                        osvvec_qmu[len(osvvec_qmu)-1].SetXYZ(t.SVOverlap_x[osvidx[int(m/2)]]-t.PV_x, t.SVOverlap_y[osvidx[int(m/2)]]-t.PV_y, t.SVOverlap_z[osvidx[int(m/2)]]-t.PV_z)
                        osvidx_qmu.append(osvidx[int(m/2)])

    dmuidxs_all = dmuidxs_osv+dmuidxs
    dmuvec_all = dmuvec_osv+dmuvec
    svidx_all = osvidx+svidx
    svvec_all = osvvec+svvec
    # If multiple muon pairs are found not from overlapping SVs, create four-muon system from non-overlapping SVs
    if len(dmuidxs_all)>=4:
        for m in range(len(dmuidxs_all)):
            if m%2>0:
                continue
            if len(qmuidxs)>3:
                break
            if m in dmuidxs_osv and m in qmuidxs_osv and len(qmuidxs_osv)>3:
                continue
            for mm in range(m+2,len(dmuidxs_all)):
                if mm%2>0:
                    continue
                if len(qmuidxs)>3:
                    break
                if mm in dmuidxs_osv and mm in qmuidxs_osv and len(qmuidxs_osv)>3:
                    continue
                if svidx_all[int(m/2)] in osvidx_qmu or svidx_all[int(mm/2)] in osvidx_qmu:
                    continue
                else:
                    qmuidxs.append(dmuidxs_all[m])
                    qmuidxs.append(dmuidxs_all[m+1])
                    qmuidxs.append(dmuidxs_all[mm])
                    qmuidxs.append(dmuidxs_all[mm+1])
                    qmuvec.append(dmuvec_all[int(m/2)])
                    qmuvec[len(qmuvec)-1] = qmuvec[len(qmuvec)-1]+dmuvec_all[int(mm/2)]
                    if svvec_all[int(m/2)].Perp() < svvec_all[int(mm/2)].Perp():
                        svvecminlxy_qmu.append(svvec_all[int(m/2)])
                        svidxminlxy_qmu.append(svidx_all[int(m/2)])
                        qmu_dmuvecminlxy.append(dmuvec_all[int(m/2)])
                        svvecmaxlxy_qmu.append(svvec_all[int(mm/2)])
                        svidxmaxlxy_qmu.append(svidx_all[int(mm/2)])
                        qmu_dmuvecmaxlxy.append(dmuvec_all[int(mm/2)])
                    else:
                        svvecminlxy_qmu.append(svvec_all[int(mm/2)])
                        svidxminlxy_qmu.append(svidx_all[int(mm/2)])
                        qmu_dmuvecminlxy.append(dmuvec_all[int(mm/2)])
                        svvecmaxlxy_qmu.append(svvec_all[int(m/2)])
                        svidxmaxlxy_qmu.append(svidx_all[int(m/2)])
                        qmu_dmuvecmaxlxy.append(dmuvec_all[int(m/2)])

    # Apply selections and fill histograms for muon pairs from non-overlapping SVs
    seldmuidxs = []
    seldmusvidxs = []
    for vn,v in enumerate(dmuvec):
        if args.noDiMuon:
            break
        if not applyDiMuonSelection(v):
            continue
        lxy  = t.SV_lxy[svidx[vn]]
        if not applyLxySelection(lxy):
            continue
        seldmuidxs.append(dmuidxs[int(vn*2)])
        seldmuidxs.append(dmuidxs[int(vn*2)+1])
        seldmusvidxs.append(svidx[vn])
        mass = v.M()
        pt   = v.Pt()
        drmm = t.Muon_vec[dmuidxs[int(vn*2)]].DeltaR(t.Muon_vec[dmuidxs[int(vn*2)+1]])
        dpmm = abs(t.Muon_vec[dmuidxs[int(vn*2)]].DeltaPhi(t.Muon_vec[dmuidxs[int(vn*2)+1]]))
        demm = abs(t.Muon_vec[dmuidxs[int(vn*2)]].Eta()-t.Muon_vec[dmuidxs[int(vn*2)+1]].Eta())
        dedpmm = 1e6
        if dpmm>0.0:
            dedpmm = demm/dpmm
        a3dmm = abs(t.Muon_vec[dmuidxs[int(vn*2)]].Angle(t.Muon_vec[dmuidxs[int(vn*2)+1]].Vect()))
        dphisv = abs(v.Vect().DeltaPhi(svvec[vn]))
        detasv = abs(v.Vect().Eta()-svvec[vn].Eta())
        detadphisv = 1e6
        if dphisv>0.0:
            detadphisv = detasv/dphisv
        a3dsv  = abs(v.Vect().Angle(svvec[vn]))
        for h in h1d:
            tn = h.GetName()
            if "hdimuon_" not in tn or "osv" in tn:
                continue
            else:
                h.Fill(eval(variable1d[h.GetName()]))
        for h in h2d:
            tn = h.GetName()
            if "hdimuon_" not in tn or "osv" in tn:
                continue
            else:
                h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]))

    # Apply selections and fill histograms for muon pairs from overlapping SVs
    seldmuidxs_osv = []
    seldmusvidxs_osv = []
    for vn,v in enumerate(dmuvec_osv):
        if args.noDiMuon:
            break
        if not applyDiMuonSelection(v):
            continue
        lxy  = t.SVOverlap_lxy[osvidx[vn]]
        if not applyLxySelection(lxy):
            continue
        seldmuidxs_osv.append(dmuidxs_osv[int(vn*2)])
        seldmuidxs_osv.append(dmuidxs_osv[int(vn*2)+1])
        seldmusvidxs_osv.append(t.SVOverlap_vtxIdxs[osvidx[vn]][0])
        mass = v.M()
        pt   = v.Pt()
        drmm = t.Muon_vec[dmuidxs_osv[int(vn*2)]].DeltaR(t.Muon_vec[dmuidxs_osv[int(vn*2)+1]])
        dpmm = abs(t.Muon_vec[dmuidxs_osv[int(vn*2)]].DeltaPhi(t.Muon_vec[dmuidxs_osv[int(vn*2)+1]]))
        demm = abs(t.Muon_vec[dmuidxs_osv[int(vn*2)]].Eta()-t.Muon_vec[dmuidxs_osv[int(vn*2)+1]].Eta())
        dedpmm = 1e6
        if dpmm>0.0:
            dedpmm = demm/dpmm
        a3dmm = abs(t.Muon_vec[dmuidxs_osv[int(vn*2)]].Angle(t.Muon_vec[dmuidxs_osv[int(vn*2)+1]].Vect()))
        dphisv = abs(v.Vect().DeltaPhi(osvvec[vn]))
        detasv = abs(v.Vect().Eta()-osvvec[vn].Eta())
        detadphisv = 1e6
        if dphisv>0.0:
            detadphisv = detasv/dphisv
        a3dsv  = abs(v.Vect().Angle(osvvec[vn]))
        for h in h1d:
            tn = h.GetName()
            if "hdimuon_" not in tn:
                continue
            else:
                h.Fill(eval(variable1d[h.GetName()]))
        for h in h2d:
            tn = h.GetName()
            if "hdimuon_" not in tn:
                continue
            else:
                h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]))

    # Fill histograms for selected SVs (with a selected muon pair)
    seldmusvidxs_all = seldmusvidxs_osv+seldmusvidxs
    nSVselass = 0
    nSVselass_osv = 0
    for v in set(seldmusvidxs_all):
        nSVselass = nSVselass+1
        lxy = t.SV_lxy[v]
        for h in h1d:
            tn = h.GetName()
            if "hsvselass_" not in tn or "osv" in hn:
                continue
            else:
                h.Fill(eval(variable1d[h.GetName()]))
        for h in h2d:
            tn = h.GetName()
            if "hsvselass_" not in tn or "osv" in hn:
                continue
            else:
                h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]))
        if v in seldmusvidxs_osv:
            nSVselass_osv = nSVselass_osv+1
            for h in h1d:
                tn = h.GetName()
                if "hsvselass_osv" not in tn:
                    continue
                else:
                    h.Fill(eval(variable1d[h.GetName()]))
            for h in h2d:
                tn = h.GetName()
                if "hsvselass_osv" not in tn:
                    continue
                else:
                    h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]))
    nSVs = nSVselass
    for h in h1d:
        if args.noDiMuon:
            break
        tn = h.GetName()
        if "h_nsvselass" not in tn or "osv" in tn:
            continue
        else:
            h.Fill(eval(variable1d[h.GetName()]))
    nSVs = nSVselass_osv
    for h in h1d:
        if args.noDiMuon:
            break
        tn = h.GetName()
        if "h_nsvselass_osv" not in tn:
            continue
        else:
            h.Fill(eval(variable1d[h.GetName()]))

    # Apply selections and fill histograms for four-muon systems from overlapping SVs
    selqmusvidxs_osv = []
    mindrmm, mindpmm, mindemm, mindedpmm, mina3dmm = 1e6, 1e6, 1e6, 1e6, 1e6
    maxdrmm, maxdpmm, maxdemm, maxdedpmm, maxa3dmm = -1., -1., -1., -1., -1.
    for vn,v in enumerate(qmuvec_osv):
        if args.noFourMuonOSV:
            break
        if not applyDiMuonSelectionForFourMuonOSV(qmu_dmuvec_osv[vn*2], qmu_dmuvec_osv[vn*2+1]):
            continue
        if not applyFourMuonSelection(v):
            continue
        lxy  = t.SVOverlap_lxy[osvidx_qmu[vn]]
        if not applyLxySelection(lxy):
            continue
        qmuidxs_osv_sel.append(qmuidxs_osv[vn*4])
        qmuidxs_osv_sel.append(qmuidxs_osv[vn*4+1])
        qmuidxs_osv_sel.append(qmuidxs_osv[vn*4+2])
        qmuidxs_osv_sel.append(qmuidxs_osv[vn*4+3])
        selqmusvidxs_osv.append(t.SVOverlap_vtxIdxs[osvidx_qmu[vn]][0])
        mass = v.M()
        pt   = v.Pt()
        for m in range(vn*4,vn*4+4):
            for mm in range(m+1,vn*4+4):
                drmm = t.Muon_vec[qmuidxs_osv[m]].DeltaR(t.Muon_vec[qmuidxs_osv[mm]])
                dpmm = abs(t.Muon_vec[qmuidxs_osv[m]].DeltaPhi(t.Muon_vec[qmuidxs_osv[mm]]))
                demm = abs(t.Muon_vec[qmuidxs_osv[m]].Eta()-t.Muon_vec[qmuidxs_osv[mm]].Eta())
                dedpmm = 1e6
                if dpmm>0.0:
                    dedpmm = demm/dpmm
                a3dmm = abs(t.Muon_vec[qmuidxs_osv[m]].Angle(t.Muon_vec[qmuidxs_osv[mm]].Vect()))
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
                if dedpmm<mindedpmm:
                    mindedpmm = dedpmm
                if dedpmm>maxdedpmm:
                    maxdedpmm = dedpmm
                if a3dmm<mina3dmm:
                    mina3dmm = a3dmm
                if a3dmm>maxa3dmm:
                    maxa3dmm = a3dmm
        dphisv = abs(v.Vect().DeltaPhi(osvvec_qmu[vn]))
        detasv = abs(v.Vect().Eta()-osvvec_qmu[vn].Eta())
        detadphisv = 1e6
        if dphisv>0.0:
            detadphisv = detasv/dphisv
        a3dsv  = abs(v.Vect().Angle(osvvec_qmu[vn]))
        for h in h1d:
            tn = h.GetName()
            if "hfourmuon_osv_" not in tn:
                continue
            else:
                h.Fill(eval(variable1d[h.GetName()]))
        for h in h2d:
            tn = h.GetName()
            if "hfourmuon_osv_" not in tn:
                continue
            else:
                h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]))

    # Fill histograms for selected overlapping SVs (with a selected four muon system)
    nSVselass_qmu_osv = 0
    for v in set(selqmusvidxs_osv):
        nSVselass_qmu_osv = nSVselass_qmu_osv+1
        lxy = t.SV_lxy[v]
        for h in h1d:
            tn = h.GetName()
            if "hsvselass_fourmu_osv" not in tn:
                continue
            else:
                h.Fill(eval(variable1d[h.GetName()]))
        for h in h2d:
            tn = h.GetName()
            if "hsvselass_fourmu_osv" not in tn:
                continue
            else:
                h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]))
    nSVs = nSVselass_qmu_osv
    for h in h1d:
        if args.noFourMuonOSV:
            break
        tn = h.GetName()
        if "h_nsvselass_fourmu_osv" not in tn:
            continue
        else:
            h.Fill(eval(variable1d[h.GetName()]))

    # Apply selections and fill histograms for four-muon systems from non-overlapping SVs
    selqmusvidxs = []
    mindrmm, mindpmm, mindemm, mindedpmm, mina3dmm = 1e6, 1e6, 1e6, 1e6, 1e6
    maxdrmm, maxdpmm, maxdemm, maxdedpmm, maxa3dmm = -1., -1., -1., -1., -1.
    for vn,v in enumerate(qmuvec):
        if args.noFourMuon:
            break
        if not applyDiMuonSelectionForFourMuon(qmu_dmuvecminlxy[vn], qmu_dmuvecmaxlxy[vn]):
            continue
        if not applyFourMuonSelection(v):
            continue
        minlxy  = svvecminlxy_qmu[vn].Perp()
        maxlxy  = svvecmaxlxy_qmu[vn].Perp()
        if not applyFourMuonLxySelection(minlxy,maxlxy):
            continue
        qmuidxs_sel.append(qmuidxs[vn*4])
        qmuidxs_sel.append(qmuidxs[vn*4+1])
        qmuidxs_sel.append(qmuidxs[vn*4+2])
        qmuidxs_sel.append(qmuidxs[vn*4+3])
        selqmusvidxs.append(svidxminlxy_qmu[vn])
        selqmusvidxs.append(svidxmaxlxy_qmu[vn])
        mass    = v.M()
        minmass = min(qmu_dmuvecminlxy[vn].M(), qmu_dmuvecmaxlxy[vn].M())
        maxmass = max(qmu_dmuvecminlxy[vn].M(), qmu_dmuvecmaxlxy[vn].M())
        avgmass = 0.5*(minmass+maxmass)
        reldmass= (maxmass-minmass)/avgmass
        pt    = v.Pt()
        minpt = min(qmu_dmuvecminlxy[vn].Pt(), qmu_dmuvecmaxlxy[vn].Pt())
        maxpt = max(qmu_dmuvecminlxy[vn].Pt(), qmu_dmuvecmaxlxy[vn].Pt())
        for m in range(vn*4,vn*4+4):
            if not m%2==0:
                continue
            drmm = t.Muon_vec[qmuidxs[m]].DeltaR(t.Muon_vec[qmuidxs[m+1]])
            dpmm = abs(t.Muon_vec[qmuidxs[m]].DeltaPhi(t.Muon_vec[qmuidxs[m+1]]))
            demm = abs(t.Muon_vec[qmuidxs[m]].Eta()-t.Muon_vec[qmuidxs[m+1]].Eta())
            dedpmm = 1e6
            if dpmm>0.0:
                dedpmm = demm/dpmm
            a3dmm = abs(t.Muon_vec[qmuidxs[m]].Angle(t.Muon_vec[qmuidxs[m+1]].Vect()))
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
            if dedpmm<mindedpmm:
                mindedpmm = dedpmm
            if dedpmm>maxdedpmm:
                maxdedpmm = dedpmm
            if a3dmm<mina3dmm:
                mina3dmm = a3dmm
            if a3dmm>maxa3dmm:
                maxa3dmm = a3dmm
        mindphisv = min(abs(qmu_dmuvecminlxy[vn].Vect().DeltaPhi(svvecminlxy_qmu[vn])),abs(qmu_dmuvecmaxlxy[vn].Vect().DeltaPhi(svvecmaxlxy_qmu[vn])))
        maxdphisv = max(abs(qmu_dmuvecminlxy[vn].Vect().DeltaPhi(svvecminlxy_qmu[vn])),abs(qmu_dmuvecmaxlxy[vn].Vect().DeltaPhi(svvecmaxlxy_qmu[vn])))
        mindetasv = min(abs(qmu_dmuvecminlxy[vn].Vect().Eta()-svvecminlxy_qmu[vn].Eta()),abs(qmu_dmuvecmaxlxy[vn].Vect().Eta()-svvecmaxlxy_qmu[vn].Eta()))
        maxdetasv = max(abs(qmu_dmuvecminlxy[vn].Vect().Eta()-svvecminlxy_qmu[vn].Eta()),abs(qmu_dmuvecmaxlxy[vn].Vect().Eta()-svvecmaxlxy_qmu[vn].Eta()))
        mindetadphisv = min(abs(qmu_dmuvecminlxy[vn].Vect().Eta()-svvecminlxy_qmu[vn].Eta())/abs(qmu_dmuvecminlxy[vn].Vect().DeltaPhi(svvecminlxy_qmu[vn])),
                            abs(qmu_dmuvecmaxlxy[vn].Vect().Eta()-svvecmaxlxy_qmu[vn].Eta())/abs(qmu_dmuvecmaxlxy[vn].Vect().DeltaPhi(svvecmaxlxy_qmu[vn])))
        maxdetadphisv = max(abs(qmu_dmuvecminlxy[vn].Vect().Eta()-svvecminlxy_qmu[vn].Eta())/abs(qmu_dmuvecminlxy[vn].Vect().DeltaPhi(svvecminlxy_qmu[vn])),
                            abs(qmu_dmuvecmaxlxy[vn].Vect().Eta()-svvecmaxlxy_qmu[vn].Eta())/abs(qmu_dmuvecmaxlxy[vn].Vect().DeltaPhi(svvecmaxlxy_qmu[vn])))
        mina3dsv  = min(abs(qmu_dmuvecminlxy[vn].Vect().Angle(svvecminlxy_qmu[vn])),abs(qmu_dmuvecmaxlxy[vn].Vect().Angle(svvecmaxlxy_qmu[vn])))
        maxa3dsv  = max(abs(qmu_dmuvecminlxy[vn].Vect().Angle(svvecminlxy_qmu[vn])),abs(qmu_dmuvecmaxlxy[vn].Vect().Angle(svvecmaxlxy_qmu[vn])))
        for h in h1d:
            tn = h.GetName()
            if "hfourmuon_" not in tn or "_osv" in tn:
                continue
            else:
                h.Fill(eval(variable1d[h.GetName()]))
        for h in h2d:
            tn = h.GetName()
            if "hfourmuon_" not in tn or "_osv" in tn:
                continue
            else:
                h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]))

    # Fill histograms for selected non-overlapping SVs (with a selected four muon system)
    nSVselass_qmu = 0
    for v in set(selqmusvidxs):
        nSVselass_qmu = nSVselass_qmu+1
        lxy = t.SV_lxy[v]
        for h in h1d:
            tn = h.GetName()
            if "hsvselass_fourmu" not in tn or "osv" in tn:
                continue
            else:
                h.Fill(eval(variable1d[h.GetName()]))
        for h in h2d:
            tn = h.GetName()
            if "hsvselass_fourmu" not in tn or "osv" in tn:
                continue
            else:
                h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]))
    nSVs = nSVselass_qmu
    for h in h1d:
        if args.noFourMuon:
            break
        tn = h.GetName()
        if "h_nsvselass_fourmu" not in tn or "osv" in tn:
            continue
        else:
            h.Fill(eval(variable1d[h.GetName()]))

    # Fill histograms for muons from selected dimuon and four-muon systems
    selmuidxs_dmu = seldmuidxs+seldmuidxs_osv
    for m in selmuidxs_dmu:
        for h in h1d:
            tn = h.GetName()
            if "hselmuon_" not in tn or "fourmu" in tn or "osv" in tn:
                continue
            else:
                h.Fill(eval(variable1d[h.GetName()]))
        for h in h2d:
            tn = h.GetName()
            if "hselmuon_" not in tn or "fourmu" in tn or "osv" in tn:
                continue
            else:
                h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]))
        if m in seldmuidxs_osv:
            for h in h1d:
                tn = h.GetName()
                if "hselmuon_osv" not in tn or "fourmu" in tn:
                    continue
                else:
                    h.Fill(eval(variable1d[h.GetName()]))
            for h in h2d:
                tn = h.GetName()
                if "hselmuon_osv" not in tn or "fourmu" in tn:
                    continue
                else:
                    h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]))

    if len(qmuidxs_sel)>3:
        for m in qmuidxs_sel:
            for h in h1d:
                tn = h.GetName()
                if "hselmuon_fourmu" not in tn or "osv" in tn:
                    continue
                else:
                    h.Fill(eval(variable1d[h.GetName()]))
            for h in h2d:
                tn = h.GetName()
                if "hselmuon_fourmu" not in tn or "osv" in tn:
                    continue
                else:
                    h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]))

    if len(qmuidxs_osv_sel)>3:
        for m in qmuidxs_osv_sel:
            for h in h1d:
                tn = h.GetName()
                if "hselmuon_fourmu_osv" not in tn:
                    continue
                else:
                    h.Fill(eval(variable1d[h.GetName()]))
            for h in h2d:
                tn = h.GetName()
                if "hselmuon_fourmu_osv" not in tn:
                    continue
                else:
                    h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]))

### Write histograms
foname = "%s/histograms_%s_all.root"%(outdir,args.year)
if args.inSample!="*":
    if args.inFile!="*":
        foname = "%s/histograms_file%s_%s_%s"%(outdir,args.inFile,args.inSample,args.year)
    else:
        foname = "%s/histograms_%s_%s"%(outdir,args.inSample,args.year)
if index>=0:
    foname = foname+("_%d"%index)
fout = ROOT.TFile(foname+".root","RECREATE")
fout.cd()
for h in hall:
    h.Write()
fout.Close()
