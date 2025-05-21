import os,sys,math
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import ROOT

def getLowerMask(lowerBound):
    rel_sigma = 0.018
    return lowerBound / (1 + 5*rel_sigma)

def getUpperMask(upperBound):
    rel_sigma = 0.018
    return upperBound / (1 - 5*rel_sigma)


ROOT.gROOT.SetBatch(1)
drawObserved = True
drawPoints = True
maskSMResonances = True
NORMCONST = 0.01 # Needs to be consistent with what it was put on make_datacards.py e.g. if signal was scaled by NORMCONST xsec should be divided by it
typeOfLimit = "BRH" # "r", "xsec", "xsecBR" "BRH" 
xsec_base = 1.0 # in pb, used to normalize the MC
xsec_h = 59.8 # higgs cross section in pb at 13.6 GeV, used to normalize the MC
scaleToFullLumi = False
compare = True
luminosity = 35
doSmartScaling = True

model = sys.argv[1]
limdir = sys.argv[2]
mass = sys.argv[3]
year = sys.argv[4]
if len(sys.argv) > 5:
    limtype = sys.argv[5]
else:
    limtype = 'asymptotic'

print("> Limits for model %s"%model)
print("> Selected mass %.3f"%float(mass))
print("> Using limits computed with %s"%limtype)

if 'Scenario' in model:
    xsec_h = 52.23
else:
    xsec_h = 59.8

if year=='2022':
    luminosity = 35
elif year=='2023':
    luminosity = 27
else:
    luminosity = 27 + 35

ctaul = []
obsl  = []
expl  = []
m2sl  = []
m1sl  = []
p1sl  = []
p2sl  = []
mext  = []
pext  = []

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
    if "Scenario" in model:
        ylabel = "95% CL upper limit on Br(h#rightarrow#Psi#Psi)"
elif typeOfLimit=="xsecBR":
    ylabel = "95% CL upper limit on #sigmaxBR"
    if model=="HTo2ZdTo2mu2x":
        ylabel = "95% CL upper limit on #sigma(h#rightarrowZ_{D}Z_{D})xB(Z_{D}#rightarrow#mu#mu) [pb]"

fin = open("%s/limits_%s_%s_%s.txt"%(limdir,model,limtype,year),"r")
for l in fin.readlines():
    if l.startswith("#"):
        continue
    ls = l.split(",")
    if ls[0]!=model:
        continue
    if "Scenario" not in model:
        if float(ls[1])!=float(mass):
            continue
    else:
        if float(ls[1])!=float(mass.split(',')[0]):
            continue
        if float(ls[2])!=float(mass.split(',')[1]):
            continue
    scale = 1.0
    # Smart scaling:
    if doSmartScaling:
        if model=="HTo2ZdTo2mu2x":
            scalingFile = '/ceph/cms/store/user/fernance/Run3ScoutingOutput/limits_HTo2ZdTo2mu2x_NormSmart_Apr-28-2025_vsCTau_asymptotic_allEras/limits_HTo2ZdTo2mu2x_allEras.txt'
            with open(scalingFile) as fin:
                for ll in fin.readlines():
                    if ll.startswith("#"):
                        continue
                    lls = ll.split(",")
                    if (float(lls[1])== float(ls[1])) and (float(lls[2])== float(ls[2])):
                        e2m = float(lls[5])
                        NORMCONST = 0.01 * ( e2m / 0.6 )
                        print(" -> Using smart scaling of %.2f for em2 limit of %.3f  to be 0.6" % (NORMCONST, e2m))
                        xsec = xsec_base*NORMCONST
                        break
    else:
        xsec = xsec_base*NORMCONST
    # Branching ratio:
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
    if typeOfLimit=="xsec":
        scale = xsec
    if typeOfLimit=="BRH":
        scale = xsec/xsec_h
    elif typeOfLimit=="xsecBR":
        if model=="HTo2ZdTo2mu2x":
            scale = xsec*BR
        print(float(ls[1]), float(ls[4]), BR)
    #if scaleToFullLumi:
    #    scale = scale / math.sqrt(10.0)
    #if model=="HTo2ZdTo2mu2x" and float(ctau) > 10 and float(ls[1]) < 1.5:
    #    continue
    if "Scenario" not in model:
        ctaul.append(float(ls[2]))
        obsl .append(scale*float(ls[3]))
        expl .append(scale*float(ls[4]))
        m2sl .append(scale*float(ls[5]))
        m1sl .append(scale*float(ls[6]))
        p1sl .append(scale*float(ls[7]))
        p2sl .append(scale*float(ls[8]))
    else:
        ctaul.append(float(ls[3]))
        obsl .append(scale*float(ls[4]))
        expl .append(scale*float(ls[5]))
        m2sl .append(scale*float(ls[6]))
        m1sl .append(scale*float(ls[7]))
        p1sl .append(scale*float(ls[8]))
        p2sl .append(scale*float(ls[9]))

    if scaleToFullLumi:
        mext .append(scale*float(ls[4])/10.0)
        pext .append(scale*float(ls[4])/math.sqrt(10.0))

fin.close()

print(ctaul)
ctauv = np.array(ctaul,"d")
obsv  = np.array(obsl ,"d")
expv  = np.array(expl ,"d")
m1sv  = np.array(m1sl ,"d")
p1sv  = np.array(p1sl ,"d")
m2sv  = np.array(m2sl ,"d")
p2sv  = np.array(p2sl ,"d")

# To cm
ctauv = ctauv * 0.1

# Comparisons
cfile = None
if compare:
    print("> Request to open files to compare")
    try:
        if model=="HTo2ZdTo2mu2x":
            if typeOfLimit=="BRH":
                if float(mass)==20.0:
                    cfile = ROOT.TFile("data/run-3/HEPData-ins2760892-v2-Figure_15_top_right.root")
                    cgobs = cfile.Get("Figure 15 top right/Graph1D_y1")
                    print(cgobs)
                if float(mass)==30.0:
                    cfile = ROOT.TFile("data/run-3/HEPData-ins2760892-v2-Figure_15_center_left.root")
                    cgobs = cfile.Get("Figure 15 center left/Graph1D_y1")
                    print(cgobs)
                if float(mass)==40.0:
                    cfile = ROOT.TFile("data/run-3/HEPData-ins2760892-v2-Figure_15_center_right.root")
                    cgobs = cfile.Get("Figure 15 center right/Graph1D_y1")
                    print(cgobs)
                if float(mass)==50.0:
                    cfile = ROOT.TFile("data/run-3/HEPData-ins2760892-v2-Figure_15_bottom_left.root")
                    cgobs = cfile.Get("Figure 15 bottom left/Graph1D_y1")
                    print(cgobs)
    except NameError:
        print("> Files not available")
        pass

# Plot the limit
plt.style.use(hep.style.CMS)
#
fig, ax = plt.subplots(figsize=(11, 8))
#
if drawObserved:
    ax.plot(ctauv, obsv, 'k-o', label='Observed', linewidth=2)
#
ax.plot(ctauv, expv, 'r--', label='Expected', linewidth=2)
#
ax.fill_between(ctauv, m2sv, p2sv, color='orange', label=r'$\pm 2\sigma$ expected')
ax.fill_between(ctauv, m1sv, p1sv, color='green', label=r'$\pm 1\sigma$ expected')
#

if compare:
    if cfile is not None:
        m_cgobs = np.array(cgobs.GetX())
        l_cgobs = np.array(cgobs.GetY())
        ax.plot(m_cgobs, l_cgobs, linestyle='solid', color='deepskyblue', label=r'Displaced dimuon with standard streams [JHEP 05 (2024) 047]', linewidth=2)

# Other details
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'Lifetime $c\tau$ [cm]')
if model=="HTo2ZdTo2mu2x":
    ax.set_ylabel(r'95% CL upper limit on Br($h \rightarrow Z_d Z_d$)')
elif "Scenario" in model:
    ax.set_ylabel(r'95% CL upper limit on Br($h \rightarrow \Psi\Psi$)')
ax.set_xlim(ctauv[0], ctauv[-1])
ax.set_ylim(0.5*min(m2sv), 50.0*max(p2sv))
#ax.set_ylim(1e-6, 1e-2)
#
if year!='allEras':
    hep.cms.label(loc=0, data=True, llabel="Preliminary", lumi=luminosity, year=year, com=13.6)
else:
    hep.cms.label(loc=0, data=True, llabel="Preliminary", lumi=luminosity, com=13.6)
#
if model=="HTo2ZdTo2mu2x":
    legend = ax.legend(loc='upper right', title=r"$H\rightarrow Z_DZ_D$ ($m_{Z_D} =$ %.1f GeV, Br($Z_D \rightarrow \mu\mu$) = %.3f)"%(float(mass), float(BR)), fontsize=15, title_fontsize=15, frameon = True)
elif model=="ScenarioA":
    legend = ax.legend(loc='upper right', title=r"Scenario A ($m_{\pi} =$ %.2f GeV, $m_{A} =$ %.2f GeV)"%(float(mass.split(',')[0]), float(mass.split(',')[1])), fontsize=15, title_fontsize=15, frameon = True)
elif model=="ScenarioB1":
    legend = ax.legend(loc='upper right', title=r"Scenario B ($m_{\pi} =$ %.2f GeV, $m_{A} =$ %.2f GeV)"%(float(mass.split(',')[0]), float(mass.split(',')[1])), fontsize=15, title_fontsize=15, frameon = True)
#legend.get_title().set_ha('left') 
legend.get_title().set_weight('bold')
legend._legend_box.align = "left"
#
if "Scenario" not in model:
    if drawObserved:
        fig.savefig("%s/limits_%s_mass%s_%s_%s_obs"%(limdir,model,("%.1f"%float(mass)).replace('.','p'),typeOfLimit, limtype), dpi=140)
    else:
        fig.savefig("%s/limits_%s_mass%s_%s_%s"%(limdir,model,("%.1f"%float(mass)).replace('.','p'),typeOfLimit, limtype), dpi=140)
else:
    if drawObserved:
        fig.savefig("%s/limits_%s_mass%s_%s_%s_%s_obs"%(limdir,model,mass.split(',')[0].replace('.','p'),mass.split(',')[1].replace('.','p'),typeOfLimit, limtype), dpi=140)
    else:
        fig.savefig("%s/limits_%s_mass%s_%s_%s_%s"%(limdir,model,mass.split(',')[0].replace('.','p'),mass.split(',')[1].replace('.','p'),typeOfLimit, limtype), dpi=140)


