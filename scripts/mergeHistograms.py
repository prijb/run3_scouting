import ROOT
import os,sys
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--inDirs", default=[], help="Choose input directories with histograms to be merged.")
parser.add_argument("--inSample", default="*", help="Choose sample; for all samples in input directory, choose '*'")
parser.add_argument("--outDir", default=os.environ.get("PWD")+"/outputHistograms_"+today, help="Choose output directory. Default: '"+os.environ.get("PWD")+"/outputHistograms_"+today+"'")
parser.add_argument("--allHistos", default=False, action="store_true", help="Merge histograms binned in dimuon lxy")
parser.add_argument("--dimuonLxy", default=False, action="store_true", help="Merge histograms binned in dimuon lxy")
parser.add_argument("--dimuonMass", default=False, action="store_true", help="Merge histograms binned in dimuon mass")
args = parser.parse_args()

if len(sys.inDirs)<2:
    print("Please, specify at least two input directories containing histograms to be merged.")
    exit()

outdir = args.outDir
if not os.path.exists(outdir):
    os.makedirs(outdir)

din = args.inDir
fin = []
for f in os.listdir(din[0]):
    if args.inSample!="*" and args.inSample not in f:
        continue
    for d in din:
        if os.path.isfile("%s/%s"%(d,f)):
            fin.append(ROOT.TFile("%s/%s"%(d,f)))
    fout = ROOT.TFile("%s/%s"%(outdir,f), "RECREATE")
    nhistos = 0
    if len(fin)>0:
        listkeys = fin[0].GetListOfKeys()
        nhistos  = listkeys.GetSize()
    histos   = []
    for h in range(0,nhistos):
        if "TH" not in listkeys.At(h).GetClassName():
            continue
        histos.append(listkeys.At(h).GetName())
    for hn in histos:
        domerge = False
        if args.allHistos:
            domerge = True
        elif args.dimuonLxy:
            domerge = (("dimuon"  in hn) or 
                       ("fourmu"  in hn and "osv" in hn) or
                       ("selass"  in hn and ("fourmu" not in hn or "osv" in hn)) or
                       ("selmuon" in hn and ("fourmu" not in hn or "osv" in hn)))
        elif args.dimuonMass:
            domerge = (("dimuon" in hn or "selass"  in hn or "selmuon" in hn) and 
                       ("fourmu" not in hn))
        fout.cd()
        th = fin[0].Get(hn).Clone()
        for fn in range(1, len(fin)):
            oh = fin[fn].Get(hn).Clone("%s_%d"%(hn,fn))
            if domerge:
                th.Add(oh)
        th.Write()
        fout.Close()
