import os,sys
import ROOT
import numpy as np
from datetime import date
import mplhep as hep
import matplotlib.pyplot as plt

user = os.environ.get("USER")
today= date.today().strftime("%b-%d-%Y")

### if date ("MMM-DD-YYYY") is not specified, use today's date
inDate = today
outDate = today
thisDir = os.environ.get("PWD")
if len(sys.argv)>2:
    inDir=sys.argv[1]
    tag=sys.argv[2]
else:
    inDir  = "%s/datacards_all_Aug-12-2024_2022_Iso1HighPt"%(thisDir)
    tag = "Iso1HighPt"

fitDir = "%s/fitResults_2022"%(thisDir)

#limDir = "%s/limits_asymptotic_%s/"%(thisDir,inDate)
limDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/limits_Aug-01-2024_2022"

biasCombination = True
biasPerChannel = True
biasPerPDF = False
nToysPerChannel = 300
nToysPerCombination = 300
biasSummary = False
checkIndividualPDFs = False
plotEnvelope = False

useCategorizedSignal = True
useCategorizedBackground = True

outDir = "%s/impacts_%s/"%(thisDir,outDate)
if not os.path.exists(outDir):
    os.makedirs(outDir)

dNames = []
dNames.append("d_FourMu_sep")
dNames.append("d_FourMu_osv")
dNames.append("d_Dimuon_lxy0p0to0p2_iso0_ptlow")
dNames.append("d_Dimuon_lxy0p0to0p2_iso0_pthigh")
dNames.append("d_Dimuon_lxy0p0to0p2_iso1_ptlow")
dNames.append("d_Dimuon_lxy0p0to0p2_iso1_pthigh")
dNames.append("d_Dimuon_lxy0p2to1p0_iso0_ptlow")
dNames.append("d_Dimuon_lxy0p2to1p0_iso0_pthigh")
dNames.append("d_Dimuon_lxy0p2to1p0_iso1_ptlow")
dNames.append("d_Dimuon_lxy0p2to1p0_iso1_pthigh")
dNames.append("d_Dimuon_lxy1p0to2p4_iso0_ptlow")
dNames.append("d_Dimuon_lxy1p0to2p4_iso0_pthigh")
dNames.append("d_Dimuon_lxy1p0to2p4_iso1_ptlow")
dNames.append("d_Dimuon_lxy1p0to2p4_iso1_pthigh")
dNames.append("d_Dimuon_lxy2p4to3p1_iso0_ptlow")
dNames.append("d_Dimuon_lxy2p4to3p1_iso0_pthigh")
dNames.append("d_Dimuon_lxy2p4to3p1_iso1_ptlow")
dNames.append("d_Dimuon_lxy2p4to3p1_iso1_pthigh")
dNames.append("d_Dimuon_lxy3p1to7p0_iso0_ptlow")
dNames.append("d_Dimuon_lxy3p1to7p0_iso0_pthigh")
dNames.append("d_Dimuon_lxy3p1to7p0_iso1_ptlow")
dNames.append("d_Dimuon_lxy3p1to7p0_iso1_pthigh")
dNames.append("d_Dimuon_lxy7p0to11p0_iso0_ptlow")
dNames.append("d_Dimuon_lxy7p0to11p0_iso0_pthigh")
dNames.append("d_Dimuon_lxy7p0to11p0_iso1_ptlow")
dNames.append("d_Dimuon_lxy7p0to11p0_iso1_pthigh")
dNames.append("d_Dimuon_lxy11p0to16p0_iso0_ptlow")
dNames.append("d_Dimuon_lxy11p0to16p0_iso0_pthigh")
dNames.append("d_Dimuon_lxy11p0to16p0_iso1_ptlow")
dNames.append("d_Dimuon_lxy11p0to16p0_iso1_pthigh")
dNames.append("d_Dimuon_lxy16p0to70p0_iso0_ptlow")
dNames.append("d_Dimuon_lxy16p0to70p0_iso0_pthigh")
dNames.append("d_Dimuon_lxy16p0to70p0_iso1_ptlow")
dNames.append("d_Dimuon_lxy16p0to70p0_iso1_pthigh")
dNames.append("d_Dimuon_lxy0p0to0p2_non-pointing")
dNames.append("d_Dimuon_lxy0p2to1p0_non-pointing")
dNames.append("d_Dimuon_lxy1p0to2p4_non-pointing")
dNames.append("d_Dimuon_lxy2p4to3p1_non-pointing")
dNames.append("d_Dimuon_lxy3p1to7p0_non-pointing")
dNames.append("d_Dimuon_lxy7p0to11p0_non-pointing")
dNames.append("d_Dimuon_lxy11p0to16p0_non-pointing")
dNames.append("d_Dimuon_lxy16p0to70p0_non-pointing")


years = []
###
years.append("2022")

# Signals
sigModels = []
sigModels.append("HTo2ZdTo2mu2x")

sigMasses = []
#sigMasses = [0.5, 0.7, 1.5, 2.0, 2.5, 5.0, 6.0, 7.0, 8.0, 14.0, 16.0, 20.0, 22.0, 24.0, 30.0, 34.0, 40.0, 44.0, 50.0]
sigMasses = [5.0, 6.0, 7.0, 8.0, 14.0, 16.0, 20.0, 22.0, 24.0, 30.0, 34.0, 40.0, 44.0, 50.0]
#sigMasses = [20.0, 22.0, 24.0, 30.0, 34.0, 40.0, 50.0]
sigMasses = [2.5]
sigCTaus = [1, 10, 100]

# Injected for impacts
inj = [0, 1, 2]
rmin = -20 # Always in between -3 and -15

mean = 0.0
sigma = 0.0
for y in years:
    for s in sigModels:
        for m in sigMasses:
            for t in sigCTaus:
                for r in inj:
                    card = "%s/card_combined_%s_M%.1f_ctau%i_%s.root"%(inDir,s,m,t,y)
                    masked = []
                    for ch in range(1, 8+1): # Hardcoded (merged lxy's)
                        # Open datacard
                        f = ROOT.TFile(card)
                        print("> Opened file: %s"%(card))
                        # Retrieve workspace from file
                        w = f.Get('w')
                        print("> Workspace contents:")
                        #w.Print()
                        # Retrieve signal normalization
                        nSig = w.obj("n_exp_binch%i_proc_signal"%ch).getVal()
                        # Retrieve background normalization
                        nBG = w.obj("n_exp_final_binch%i_proc_background"%ch).getVal()
                        # Close input file with workspace
                        f.Close()
                        print('Found numbers for channel ch%i:'%ch)
                        print('> nSig = %f'%nSig)
                        print('> nBG = %i'%nBG)
                        if ((3*nSig > (nBG+r*nSig)) or nBG < 1):
                            masked.append(ch)
                            continue
                        rmin = min([-3, max([rmin, -(nBG+r*nSig)/(nSig)])])
                    print('Mask:')
                    print(masked)
                    # Define the mask:
                    mask = '--setParameters '
                    for ch in masked:
                        mask = mask + 'mask_ch%i=1,'%ch
                    mask = mask[:-1] # remove last coma
                    # Define workdir
                    subDir = '%s/%s_M%.1f_ctau%i_%s'%(outDir,s,m,t,y)
                    if not os.path.exists(subDir):
                        os.mkdir(subDir)
                    # Options
                    options = "--cminDefaultMinimizerStrategy 0 -t -1"
                    name = "_%s_M%.1f_ctau%i_%s_t%i"%(s,m,t,tag,r)
                    # Initial fit:
                    print("combine -M FitDiagnostics %s %s --expectSignal %i -n %s -m %.0f --forceRecreateNLL --rMin %.2f %s"%(card,options,r,name,m,rmin,mask))
                    os.system("combine -M FitDiagnostics %s %s --expectSignal %i -n %s -m %.0f --forceRecreateNLL --rMin %.2f %s"%(card,options,r,name,m,rmin,mask))
                    print("python3 $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a fitDiagnostics%s.root -g plots%s.root >& fitResults%s.txt"%(name,name,name))
                    os.system("python3 $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a fitDiagnostics%s.root -g plots%s.root >& fitResults%s.txt"%(name,name,name))
                    #
                    print("combineTool.py -M Impacts -d %s %s --expectSignal %i -n %s -m %.0f --doInitialFit --rMin %.2f %s"%(card,options,r,name,m,rmin,mask))
                    os.system("combineTool.py -M Impacts -d %s %s --expectSignal %i -n %s -m %.0f --doInitialFit --rMin %.2f %s"%(card,options,r,name,m,rmin,mask))
                    print("combineTool.py -M Impacts -d %s -o impacts%s.json %s --expectSignal %i -n %s -m %.0f --doFits --parallel 20 --task-name %s %s  >& /dev/null"%(card,name,options,r,name,m,name,mask))
                    os.system("combineTool.py -M Impacts -d %s -o impacts%s.json %s --expectSignal %i -n %s -m %.0f --rMin %.2f --doFits --parallel 20 --task-name %s %s >& /dev/null"%(card,name,options,r,name,m,rmin,name,mask))
                    print("combineTool.py -M Impacts -d %s -m %.0f -n %s -o impacts%s.json"%(card,m,name,name))
                    os.system("combineTool.py -M Impacts -d %s -m %.0f -n %s -o impacts%s.json"%(card,m,name,name))
                    print("plotImpacts.py -i impacts%s.json -o impacts%s"%(name,name))
                    os.system("plotImpacts.py -i impacts%s.json -o impacts%s"%(name,name))
                    # Save roots:
                    os.system('mv fitDiagnostics%s.root %s'%(name,subDir))
                    os.system('mv higgsCombine%s.FitDiagnostics.mH%.0f.root %s'%(name,m,subDir))
                    os.system('mv fitResults%s.txt %s'%(name,subDir))
                    os.system('mv higgsCombine_initialFit_%s.MultiDimFit.mH%.0f.root %s'%(name,m,subDir))
                    os.system('mv higgsCombine_paramFit_%s_*.MultiDimFit.mH%.0f.root %s'%(name,m,subDir))
                    os.system('mv plots%s.root %s'%(name,subDir))
                    os.system('mv impacts%s.json %s'%(name,subDir))
                    os.system('mv impacts%s.pdf %s'%(name,subDir))
                    os.system('cp utils/index.php %s'%(subDir))




