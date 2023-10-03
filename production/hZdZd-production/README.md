# HAHM production

The production takes as an example the existing HAHM hToZdZd simulation for Run 3 with mZd = 20 GeV and epsilon = 1e-8. Central fragment for this sample can be accessed through:

```
https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_fragment/EXO-Run3Summer22EEwmLHEGS-01017
```

Additional info:
McM record: https://cms-pdmv.cern.ch/mcm/requests?page=0&prepid=EXO-Run3Summer22EEwmLHEGS-01017


## Gridpack production

Template cards to make the gridpacks are found in ```card-templates/``` and are made from ones extracted from the central gridpack:
```
/cvmfs/cms.cern.ch/phys_generator/gridpacks/RunIII/13p6TeV/slc7_amd64_gcc10/madgraph/LL_HAHM_MS_400/LL_HAHM_MS_400_kappa_0p01_MZd_20_eps_1e-08_slc7_amd64_gcc10_CMSSW_12_4_8_tarball.tar.xz
```

To create the new cards for a predefined set of (mZd, epsilon) values fist add the chosen points in ```makeCards.py``` with the following structure:
```
model_grid = {}
model_grid["5.0"] = ['1e-8', '2e-7']
```
and then run:
```
python3 makeCards.py
``` 
The output cards will be available in the ```hZdZd/``` dir.

To create the gridpacks the central generation repository should be cloned:
```
git clone https://github.com/cms-sw/genproductions.git
```
The new cards should be copied within the ```genproductions/bin/MadGraph5_aMCatNLO/cards``` folder. Then, in ```genproductions/bin/MadGraph5_aMCatNLO/``` the gridpacks are created by running:
```
./gridpack_generation.sh LL_HAHM_MS_400_kappa_0p01_MZd_5p0_eps_1e-8 cards/hZdZd/hZdZd_mZd_5p0_eps_1e-8/
```

## Fragment production

(in progress)
