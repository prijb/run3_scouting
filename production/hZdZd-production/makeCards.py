import os
import sys

inputdir = "card-templates/"
tmp_customizecards = "LL_HAHM_MS_400_kappa_0p01_MZd_ZDMASS_eps_EPSILON_customizecards.dat"
tmp_extramodels    = "LL_HAHM_MS_400_kappa_0p01_MZd_ZDMASS_eps_EPSILON_extramodels.dat"
tmp_proc_card      = "LL_HAHM_MS_400_kappa_0p01_MZd_ZDMASS_eps_EPSILON_proc_card.dat"
tmp_run_card       = "LL_HAHM_MS_400_kappa_0p01_MZd_ZDMASS_eps_EPSILON_run_card.dat"

with open(inputdir + tmp_customizecards, 'r') as file_:
    txt_customizecards = file_.read()

with open(inputdir + tmp_extramodels, 'r') as file_:
    txt_extramodels = file_.read()

with open(inputdir + tmp_proc_card, 'r') as file_:
    txt_proc_card = file_.read()

with open(inputdir + tmp_run_card, 'r') as file_:
    txt_run_card = file_.read()

#
## Define the model grid: mass ; [epsilon]
model_grid = [] # mass : [epsilon values]
model_grid.append([5, [1e-06, 5e-07, 1e-07, 3e-08]])
model_grid.append([10, [1e-06, 5e-07, 1e-07, 3e-08]])

for point in model_grid:
    mass = str(point[0])
    for epsilon in point[1]:
        outdir = 'hZdZd/hZdZd_mZd_{ZDMASS}_eps_{EPSILON}/'.format(ZDMASS = mass.replace('.', 'p'), EPSILON = epsilon)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        customizecards = tmp_customizecards.replace('ZDMASS', '{ZDMASS}').replace('EPSILON', '{EPSILON}').format(ZDMASS = mass.replace('.', 'p'), EPSILON = epsilon)
        with open(outdir + customizecards, 'w') as file_:
            out_customizecards = txt_customizecards.format(ZDMASS = mass, EPSILON = epsilon) # not replace the '.'
            file_.write(out_customizecards)
            file_.close()
        extramodels = tmp_extramodels.replace('ZDMASS', '{ZDMASS}').replace('EPSILON', '{EPSILON}').format(ZDMASS = mass.replace('.', 'p'), EPSILON = epsilon)
        with open(outdir + extramodels, 'w') as file_:
            out_extramodels = txt_extramodels # No need to format (?)
            file_.write(out_extramodels)
            file_.close()
        proc_card = tmp_proc_card.replace('ZDMASS', '{ZDMASS}').replace('EPSILON', '{EPSILON}').format(ZDMASS = mass.replace('.', 'p'), EPSILON = epsilon)
        with open(outdir + proc_card, 'w') as file_:
            out_proc_card = txt_proc_card.format(ZDMASS = mass.replace('.', 'p'), EPSILON = epsilon)
            file_.write(out_proc_card)
            file_.close()
        run_card = tmp_run_card.replace('ZDMASS', '{ZDMASS}').replace('EPSILON', '{EPSILON}').format(ZDMASS = mass.replace('.', 'p'), EPSILON = epsilon)
        with open(outdir + run_card, 'w') as file_:
            out_run_card = txt_run_card # No need to format (?)
            file_.write(out_run_card)
            file_.close()



