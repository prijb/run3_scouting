import sys
import os

## Init the masses / epsilon (by default the same grid as in makeModel.py)
model_grid = [] # mass : [epsilon values]
model_grid.append([0.7, [5e-06, 2e-06, 5e-07, 2e-07]])
model_grid.append([1,   [5e-06, 2e-06, 5e-07, 2e-07]])
model_grid.append([2,   [3e-06, 1e-06, 3e-07, 1e-07]])
model_grid.append([3,   [3e-06, 1e-06, 3e-07, 1e-07]])
model_grid.append([5,   [1e-06, 6e-07, 2e-07, 7e-08]])
model_grid.append([7,   [1e-06, 6e-07, 2e-07, 7e-08]])
model_grid.append([10,  [1e-06, 5e-07, 1e-07, 3e-08]])


## TO BE RUN WITHIN MADGRAPH ENV
pwd = os.getcwd()
if not os.path.isdir('results_xsec/'):
    os.makedirs('results_xsec/')
outdir = pwd + '/results_xsec/'
os.system('wget https://raw.githubusercontent.com/cmstas/run3_scouting/hZdZdprod/production/hZdZd-production/card-templates/process_card.dat')
with open('process_card.dat', 'r') as f:
    process = f.read()
os.system('wget https://raw.githubusercontent.com/cmstas/run3_scouting/hZdZdprod/production/hZdZd-production/card-templates/param_card.dat')
with open('param_card.dat', 'r') as f:
    param = f.read()
for model in model_grid:
    mass = model[0]
    smass = str(model[0])
    epsilon_list = model[1]
    for epsilon in epsilon_list:
        ## Run process
        process_name = 'process_mZd_{ZDMASS}_eps_{EPSILON}_card.dat'.format(ZDMASS = smass.replace('.', 'p'), EPSILON = epsilon)
        process_txt = process.format(ZDMASS = smass.replace('.', 'p'), EPSILON = epsilon)
        process_dir = 'process_mZd_{ZDMASS}_eps_{EPSILON}'.format(ZDMASS = smass.replace('.', 'p'), EPSILON = epsilon)
        with open(process_name, 'w') as f:
            f.write(process_txt)
        os.system('./bin/mg5_aMC ' + process_name)
        ## Update run_card.dat
        get_run = 'wget -O {DIR}/Cards/run_card.dat https://raw.githubusercontent.com/cmstas/run3_scouting/hZdZdprod/production/hZdZd-production/card-templates/run_card.dat'.format(DIR = process_dir)
        os.system(get_run) 
        ## Update param_card.dat
        param_txt = param.format(mass, epsilon)
        with open(process_dir + '/Cards/param_card.dat', 'w') as f:
            f.write(param_txt)
        os.chdir(process_dir)
        os.system('./bin/generate_events -f')
        os.chdir(pwd)
        os.system(process_dir + '/HTML/run_01/results.html ' + outdir + process_dir + '.html')

