import os,sys
import json
import matplotlib.pyplot as plt
import numpy as np

useSignalMC = False

if len(sys.argv)<3:
    print("Please, specify model and limit directory.")
    exit(1)

model = sys.argv[1]
indir = sys.argv[2]
year = sys.argv[3]
outdir = 'output_prueba' 
if not os.path.exists(outdir):
    os.makedirs(outdir)

# We remove the 2
masses =  [0.5, 0.7, 1.5, 2.0, 2.5, 5.0, 6.0, 7.0, 8.0, 14.0, 16.0, 20.0, 22.0, 24.0, 30.0, 34.0, 40.0, 44.0, 50.0]
ctaus = [1, 10, 100]

### Goodness of fit (Background-only)
channel_dict = {}
channel_dict["3"] = [0,1]
channel_dict["4"] = [0,2]
channel_dict["5"] = [0,3]
channel_dict["6"] = [0,4]

channel_dict["7"] = [1,1]
channel_dict["8"] = [1,2]
channel_dict["9"] = [1,3]
channel_dict["10"] = [1,4]

channel_dict["11"] = [2,1]
channel_dict["12"] = [2,2]
channel_dict["13"] = [2,3]
channel_dict["14"] = [2,4]

channel_dict["15"] = [3,1]
channel_dict["16"] = [3,2]
channel_dict["17"] = [3,3]
channel_dict["18"] = [3,4]

channel_dict["19"] = [4,1]
channel_dict["20"] = [4,2]
channel_dict["21"] = [4,3]
channel_dict["22"] = [4,4]

channel_dict["23"] = [5,1]
channel_dict["24"] = [5,2]
channel_dict["25"] = [5,3]
channel_dict["26"] = [5,4]

channel_dict["27"] = [6,1]
channel_dict["28"] = [6,2]
channel_dict["29"] = [6,3]
channel_dict["30"] = [6,4]

channel_dict["31"] = [7,1]
channel_dict["32"] = [7,2]
channel_dict["33"] = [7,3]
channel_dict["34"] = [7,4]

cat = [r'non-isolated, $p_{T}^{\mu\mu} < 25$ GeV', r'non-isolated, $p_{T}^{\mu\mu} > 25$ GeV', r'isolated, $p_{T}^{\mu\mu} < 25$ GeV', r'isolated, $p_{T}^{\mu\mu} > 25$ GeV']
lxycat = [r'$l_{xy}\in[0.0, 0.2]$ cm, ', r'$l_{xy}\in[0.2, 1.0]$ cm, ', r'$l_{xy}\in[1.0, 2.4]$ cm, ', r'$l_{xy}\in[2.4, 3.1]$ cm, ', r'$l_{xy}\in[3.1, 7]$ cm, ', r'$l_{xy}\in[7, 11]$ cm, ', r'$l_{xy}\in[11, 16]$ cm, ', r'$l_{xy}\in[16, 7.0]$ cm, ']

fig, axes = plt.subplots(8, 4, figsize=(11, 14), sharex=True, sharey=True)

for ch in range(3, 35):
    channel = str(ch)
    p_values = []
    m_values = []
    for m in masses:
        m_values.append(m)
        mass = str(m)
        mass = mass.replace('.0', '')
        output_json = "%s/gof_%s_ch%s_M%s.json"%(outdir, model, channel,  mass)
        if not os.path.exists(output_json):
            seed = str(int(m))
            ### Generate summary
            command = "combineTool.py -M CollectGoodnessOfFit --input %s/higgsCombine_ch%s.GoodnessOfFit.mH%s.root %s/higgsCombine_ch%s.GoodnessOfFit.mH%s.%s.root -o %s -m %s"%(indir, channel, mass, indir, channel, mass, seed, output_json, mass)
            os.system(command)
        with open(output_json, 'r') as file:
            data = json.load(file)
            p_values.append(data[next(iter(data))]["p"])

    i = channel_dict[str(channel)][0]
    j = channel_dict[str(channel)][1]-1
    ax = axes[i, j]
    ax.scatter(np.array(m_values), np.array(p_values), color='red', s=4)
    ax.set_xscale('log')
    ax.set_xlim(0.4, 60.)
    ax.set_ylim(0, 1.25)
    y_ticks = [0, 0.25, 0.5, 0.75, 1.0]
    y_labels = ['0', '0.25', '0.5', '0.75', '1']
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    ax.text(0.5, 1.1, lxycat[i]+cat[j], fontsize=6)
    ax.grid(True)

fig.text(0.5, 0.04, r'Mass $m_{\mu\mu}$ (GeV)', ha='center', fontsize=14)
fig.text(0.04, 0.5, 'p-value', va='center', rotation='vertical', fontsize=14)

plt.subplots_adjust(hspace=0.0, wspace=0.0)
fig.savefig('pruebita.png', dpi=160)








