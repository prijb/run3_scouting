import ROOT
import os,sys

if len(sys.argv)<4:
    print("Please, specify two input directories containing histograms from each sideband, and name for output directory")
    exit()

outdir = sys.argv[3]
if not os.path.exists(outdir):
    os.makedirs(outdir)

din = [sys.argv[1], sys.argv[2]]
fin = []
for f in os.listdir(din[0]):
    if os.path.isfile("%s/%s"%(din[0],f)) and os.path.isfile("%s/%s"%(din[1],f)):
        fin.append(ROOT.TFile("%s/%s"%(din[0],f)))
        fin.append(ROOT.TFile("%s/%s"%(din[1],f)))
        fout = ROOT.TFile("%s/%s"%(outdir,f), "RECREATE")
        listkeys = fin[0].GetListOfKeys()
        nhistos  = listkeys.GetSize()
        histos   = []
        for h in range(0,nhistos):
            if "TH" not in listkeys.At(h).GetClassName():
                continue
            histos.append(listkeys.At(h).GetName())
        for hn in histos:
            fout.cd()
            th = fin[0].Get(hn).Clone()
            oh = fin[1].Get(hn).Clone("%s_other"%hn)
            if ( ("dimuon"  in hn or "selass"  in hn or "selmuon" in hn) and "fourmu" not in hn ):
                th.Add(oh)
                th.Write()
            else:
                th.Write()
        fout.Close()
