import ROOT

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


for d_,d in enumerate(dNames):
   if d=="d_FourMu_sep":
       binidx=1
   elif d=="d_FourMu_osv":
       binidx=2
   elif d=="d_Dimuon_lxy0p0to0p2_iso0_ptlow":
       binidx=3
   elif d=="d_Dimuon_lxy0p0to0p2_iso0_pthigh":
       binidx=4
   elif d=="d_Dimuon_lxy0p0to0p2_iso1_ptlow":
       binidx=5
   elif d=="d_Dimuon_lxy0p0to0p2_iso1_pthigh":
       binidx=6
   elif d=="d_Dimuon_lxy0p2to1p0_iso0_ptlow":
       binidx=7
   elif d=="d_Dimuon_lxy0p2to1p0_iso0_pthigh":
       binidx=8
   elif d=="d_Dimuon_lxy0p2to1p0_iso1_ptlow":
       binidx=9
   elif d=="d_Dimuon_lxy0p2to1p0_iso1_pthigh":
       binidx=10
   elif d=="d_Dimuon_lxy1p0to2p4_iso0_ptlow":
       binidx=11
   elif d=="d_Dimuon_lxy1p0to2p4_iso0_pthigh":
       binidx=12
   elif d=="d_Dimuon_lxy1p0to2p4_iso1_ptlow":
       binidx=13
   elif d=="d_Dimuon_lxy1p0to2p4_iso1_pthigh":
       binidx=14
   elif d=="d_Dimuon_lxy2p4to3p1_iso0_ptlow":
       binidx=15
   elif d=="d_Dimuon_lxy2p4to3p1_iso0_pthigh":
       binidx=16
   elif d=="d_Dimuon_lxy2p4to3p1_iso1_ptlow":
       binidx=17
   elif d=="d_Dimuon_lxy2p4to3p1_iso1_pthigh":
       binidx=18
   elif d=="d_Dimuon_lxy3p1to7p0_iso0_ptlow":
       binidx=19
   elif d=="d_Dimuon_lxy3p1to7p0_iso0_pthigh":
       binidx=20
   elif d=="d_Dimuon_lxy3p1to7p0_iso1_ptlow":
       binidx=21
   elif d=="d_Dimuon_lxy3p1to7p0_iso1_pthigh":
       binidx=22
   elif d=="d_Dimuon_lxy7p0to11p0_iso0_ptlow":
       binidx=23
   elif d=="d_Dimuon_lxy7p0to11p0_iso0_pthigh":
       binidx=24
   elif d=="d_Dimuon_lxy7p0to11p0_iso1_ptlow":
       binidx=25
   elif d=="d_Dimuon_lxy7p0to11p0_iso1_pthigh":
       binidx=26
   elif d=="d_Dimuon_lxy11p0to16p0_iso0_ptlow":
       binidx=27
   elif d=="d_Dimuon_lxy11p0to16p0_iso0_pthigh":
       binidx=28
   elif d=="d_Dimuon_lxy11p0to16p0_iso1_ptlow":
       binidx=29
   elif d=="d_Dimuon_lxy11p0to16p0_iso1_pthigh":
       binidx=30
   elif d=="d_Dimuon_lxy16p0to70p0_iso0_ptlow":
       binidx=31
   elif d=="d_Dimuon_lxy16p0to70p0_iso0_pthigh":
       binidx=32
   elif d=="d_Dimuon_lxy16p0to70p0_iso1_ptlow":
       binidx=33
   elif d=="d_Dimuon_lxy16p0to70p0_iso1_pthigh":
       binidx=34
   elif d=="d_Dimuon_lxy0p0to0p2_non-pointing":
       binidx=35
   elif d=="d_Dimuon_lxy0p2to1p0_non-pointing":
       binidx=36
   elif d=="d_Dimuon_lxy1p0to2p4_non-pointing":
       binidx=37
   elif d=="d_Dimuon_lxy2p4to3p1_non-pointing":
       binidx=38
   elif d=="d_Dimuon_lxy3p1to7p0_non-pointing":
       binidx=39
   elif d=="d_Dimuon_lxy7p0to11p0_non-pointing":
       binidx=40
   elif d=="d_Dimuon_lxy11p0to16p0_non-pointing":
       binidx=41
   elif d=="d_Dimuon_lxy16p0to70p0_non-pointing":
       binidx=42
   # Check workspace:
   year = 2022
   catExtS = ""
   catExtB = ""
   catExtS = "_ch%d_%s"%(binidx, year)
   catExtB = "_ch%d_%s"%(binidx, year)
   # Open input file with workspace
   finame = finame = "datacards_all_Oct-20-2024_2022/fitResults_2022/%s_Signal_HTo2ZdTo2mu2x_MZd-7p0_ctau-100mm_2022_workspace.root"%(d)
   f = ROOT.TFile(finame)
   # Retrieve workspace from file
   wfit = f.Get('wfit')
   # Retrieve signal and bkg normalization
   nSig = wfit.var("signalNorm%s"%catExtS).getValV()
   nBG  = wfit.data("data_obs%s"%catExtB).sumEntries()
   # Then check the output merged datacard
   combcard = 'datacards_all_Oct-20-2024_2022/card_combined_HTo2ZdTo2mu2x_M7.000_ctau100_2022.root'
   card = ROOT.TFile(combcard)
   w = card.Get('w')
   nSigCard = w.function("n_exp_binch%i_proc_signal"%(binidx)).getVal()
   nBGCard = w.function("n_exp_final_binch%i_proc_background"%(binidx)).getVal()
   nDataCard = w.data("data_obs").reduce("CMS_channel==%i"%(binidx-1)).sumEntries()
   # Print
   print(binidx, nSig, nSigCard)
   print(binidx, nBG, nBGCard)
   print(binidx, nBG, nDataCard)



