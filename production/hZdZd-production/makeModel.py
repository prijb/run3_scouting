import os
import sys
from utils import parse_table

# Constants
CHSLASH = 0.19733e-15 # GeV*m


# Get br table data for the dark photon decay
info = parse_table('model-tables/HiggsedDarkPhoton_BrTableData.txt')

########## File with mass / epsilon / ctau
#
# Init the file
model_file = open('mass_epsilon_gamma_ctau.txt', 'w')
model_file.write('mZd (GeV)     \tEpsilon    \tGammaZd (GeV)     \tctau (mm)\n')

## todo: Fix the table, looks not aligned due to differences of heading and values length

#
# Init the grid
model_grid = [] # mass : [epsilon values]
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
            gammaZd = gammaZd_over_eps * epsilon * epsilon
            sgammaZd = str(gammaZd)
            lifetime = CHSLASH * 1000. / gammaZd
            sepsilon = str(epsilon)
            slifetime = str(lifetime)
            outline = smZd + '\t' + sepsilon + '\t' + sgammaZd + '\t' +slifetime + '\n'
            model_file.write(outline)

model_file.close()
