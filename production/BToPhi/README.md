# Instructions to generate B-hadron model

Note: 2023 production not implemented as of today (10-Mar-2024)

A set of predefined fragments are included in the ```fragment-examples/``` dir, and some gen cfg files are created out of them and saved in ```signal/``` with 2022 and 2022postEE conditions. There are two ways of running the samples: condor and crab. While condor generates the AOD samples in one job, crab needs to run GEN, RAWSIM and AOD in separate ones. However, crab is more stable and allows to publish the dataset, so it is recommended.

## Before running

Initialize a given release to produce the samples e.g. for 2022/2022postEE:
```
cmsrel CMSSW_12_4_16
cd CMSSW_12_4_16/src
cmsenv
cd ../../
```
and then eveything should work as normal.

## How to generate new fragments

Modify the parameters inside ```makeFragments.py```and then run:
```
python3 makeFragments.py
```
New fragments will be saved in ```outputFragments_Month-Day-Year```

## How to run with condor

Note: Not stable for bigger grids, returning a ```"Too many open files"``` exception in condor, also when using ProjectMetis. Not fixed.

Run production in condor by executing e.g.
```
sh condor/make_tarball.sh
sh condor/produce.sh BToPhi_MPhi-2p0_ctau-10mm_v1p0 2022 BToPhi_MPhi-2p0_ctau-10mm
```
with the predefined splitting, set in ```condor/produce.sub```. Make sure that your cfg files are in ```signal/``` dir.

## How to run with crab

1. Make sure that corresponding gen cfg files are saved in ```signal/```. Then initial GENSIM dataset can be launched like:
```
python crab/crab_gen_cfg.py 2022 BToPhi_MPhi-2p0_ctau-1mm
python crab/crab_gen_cfg.py 2022 BToPhi_MPhi-2p0_ctau-10mm
python crab/crab_gen_cfg.py 2022 BToPhi_MPhi-2p0_ctau-100mm
python crab/crab_gen_cfg.py 2022postEE BToPhi_MPhi-2p0_ctau-1mm
python crab/crab_gen_cfg.py 2022postEE BToPhi_MPhi-2p0_ctau-10mm
python crab/crab_gen_cfg.py 2022postEE BToPhi_MPhi-2p0_ctau-100mm
```
Note: If we go to higher private production, consider to automatize this.

2. Modify the ```crab/crab_rawsim_cfg.py``` file with the output GENSIM datasets and run:
```
python crab/crab_rawsim_cfg.py 2022
python crab/crab_rawsim_cfg.py 2022postEE
```

3. Modify the ```crab/crab_aodsim_cfg.py``` file with the output GENSIM datasets and run:
```
python crab/crab_aodsim_cfg.py 2022
python crab/crab_aodsim_cfg.py 2022postEE
```
