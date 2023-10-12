import os
import sys
import math
import numpy as np
import pickle
#import matplotlib.pyplot as plt
from utils import parse_table

#
# Partial width for a fermion [Eq.(2.12)]
#
def computePartGamma(mZd, mf, gL, gR):
    Nc = 3.0 # Assuming 3 colors
    term1 = Nc / (24.0*math.pi*mZd)
    term2 = math.sqrt(1.0 - (4.0*mf**2)/(mZd**2))
    term3 = mZd**2 * (gL**2 + gR**2)
    term4 = mf**2 * (-6.0*gL*gR + gL**2 + gR**2)
    value = term1 * term2 * (term3 - term4)
    return value


# Constants
CHSLASH = 0.19733e-15 # GeV*m
ECOP = math.sqrt(4.0 * math.pi * 1./137.) # Taken alpha = 1/128 at Q^2 = mW^2

# fermion masses GeV
me = 0.511* 1e-3
mmu = 106 * 1e-3
mtau = 1.7768
mu = 2.2e-3
md = 4.7e-3
mc = 1.28
ms = 96e-3
mb = 4.18

# Get br table data for the dark photon decay
info = parse_table('model-tables/HiggsedDarkPhoton_BrTableData.txt')

########## File with mass / epsilon / ctau
#
# Init the file
model_file = open('mass_epsilon_gamma_ctau.txt', 'w')
model_file.write('mZd (GeV)     \tEpsilon    \tGammaZd (GeV)     \tctau (mm)\n')

## todo: Fix the table, looks not aligned due to differences of heading and values length

#brs_file.write('mZd (GeV)\t\tEpsilon\t\tBR_MUMU\t\tBR_EE\t\tBR_TAUTAU\t\tBR_UUBAR\t\tBR_DDBAR\t\t BR_SSBAR\t\tBR_CCBAR\t\tBR_BBBAR\t\tBR_MUNUMUNUBAR\t\tBR_ENUENUBAR\t\tBR_TAUNUTAUNUBAR\n')

#
# Init the grid
model_grid = [] # mass : [epsilon values]
model_grid.append([5, [1e-06, 5e-07, 1e-07, 3e-08]])
model_grid.append([10, [1e-06, 5e-07, 1e-07, 3e-08]])

#
# Init the loop
for model in model_grid:
    mZd = model[0]
    smZd = str(mZd)
    epsilon_list = model[1]
    true_ref = []
    for r,ref in enumerate(info):
        mZd_ref = ref[0]
        mZd_sref = mZd_ref[:]
        if '.' == mZd_sref[-1]:
            mZd_sref = mZd_sref[:-1]
        if smZd==mZd_ref or smZd==mZd_sref:
            true_ref = ref
            break
        if r==len(info)-1 and len(true_ref)<1:
            print('Mass %s not available, skipping...'% smZd)
    if len(true_ref) > 0:
        sgammaZd_over_eps = true_ref[1]
        gammaZd_over_eps = float(sgammaZd_over_eps)
        for epsilon in epsilon_list:
            # For model file
            gammaZd = gammaZd_over_eps * epsilon * epsilon
            sgammaZd = str(gammaZd)
            lifetime = CHSLASH * 1000. / gammaZd
            sepsilon = str(epsilon)
            slifetime = str(lifetime)
            outline = smZd + '\t' + sepsilon + '\t' + sgammaZd + '\t' +slifetime + '\n'
            model_file.write(outline)

model_file.close()

#
# Material for plot (matplotlib is not available in uaf?)
masses = ['0.7', '1.', '2.', '3.', '5.', '7.', '10.', '20']
ctau = np.logspace(-1, 4, 51)
epsilon = []
for smZd in masses:
    mZd = float(smZd)
    true_ref = []
    evalues = np.zeros(51)
    for r,ref in enumerate(info):
        mZd_ref = ref[0]
        if smZd==mZd_ref:
            true_ref = ref
            break
        if r==len(info)-1 and len(true_ref)<1:
            print('Mass %s not available, skipping...'% smZd)
    if len(true_ref) > 0:
        sgammaZd_over_eps = true_ref[1]
        gammaZd_over_eps = float(sgammaZd_over_eps)
        for c in range(0, len(ctau)):
            gammaZd = CHSLASH * 1000. / ctau[c]
            evalues[c] = (gammaZd / gammaZd_over_eps)**0.5
        epsilon.append(evalues)

data_to_save = [masses, ctau, epsilon]
# Dump the data to a file using pickle
with open('data.pkl', 'wb') as file:
    pickle.dump(data_to_save, file)
