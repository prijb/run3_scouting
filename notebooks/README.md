# Notebooks

To plot the main figures for the paper and supplementary material.

To use in uaf-2 (or other remote machine) first tunnel through ssh in your local:
```
ssh -N -L 8893:localhost:8893 username@uaf-2.t2.ucsd.edu
```

and within this directory in your remote:
```
SCRAMARCH=el8_amd64_gcc10
CMSSWVERSION=CMSSW_12_6_0
source env/bin/activate
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/$SCRAMARCH/cms/cmssw/$CMSSWVERSION/src ; eval `scramv1 runtime -sh` ; cd -
jupyter notebook --no-browser --port=8893
```
source init.sh
```

## Main paper plots


## Fit plots


## Limit plot
