# Muon scouting

Clone repository from GitHub with:
``` shell
git clone --recursive https://github.com/cmstas/run3_scouting.git
```

## Skimming

The skimming code resides in `batch/`.

The skimming code runs on RAW data sets, and produces a skimmed RAW data set as output
with additional information, mainly:
- HLT and L1T flags;
- tracker information (number of expected muon tracker layers, position of displaced vertices relative to tracker modules).

## Looper (C++)

The C++ looper resides in `cpp/`.

The looper runs on skimmed RAW data sets, and produces a flat tree as output.
Please, refer to `README` in `cpp/` for further instructions.
Condor submission is set up in `condor/`.

## Histograms:

### Histogram filling:

Histograms are filled by `fillHistosScouting.py`, and written in a ROOT output file.
This PyROOT looper optionally applies selections on multi-muon system kinematics and displacement.
Condor submission is set up in `condor/`.

Histograms are defined in `utils/histDefinition.py`:
please, add your histograms there, following the existing structure.

### Histogram plotting:

For plotting output histograms: `plotHistosScouting.py`.
E.g.:
- in order to plot/compare different eras from directory with same selection:
``` shell
python3 plotHistosScouting.py --inSamples Data --inDir /ceph/cms/store/user/mmasciov/Run3ScoutingOutput/outputHistograms_Aug-02-2023_all/ --doRatio --shape --doRatio --logY
```
- in order to plot/compare data with different selection (i.e., from different directories), e.g., with different lxy selections on the dimuon system:
``` shell
python3 plotHistosScouting.py --inSamples Data --inMultiDir /ceph/cms/store/user/mmasciov/Run3ScoutingOutput/outputHistograms_Aug-02-2023_dimuonMass2p95to3p25_relaxedSVselection_onlyDiMuon/ /ceph/cms/store/user/mmasciov/Run3ScoutingOutput/outputHistograms_Aug-02-2023_dimuonMass2p5to2p95-3p25to3p4_relaxedSVselection_onlyDiMuon/ --inMultiLeg "J/#psi" "J/#psi sidebands" --doRatio --relaxedSVSel --dimuonMassSel 2.5 3.4 --shape --doRatio --logY  --noPreSel --noFourMuon --noFourMuonOSV --outSuffix JPsi
```

In order to get help with all (optional and required) arguments, just execute the script with no argument:
``` shell
python3 plotHistosScouting.py
```

### Other histogram manipulation scripts:
Other scripts to `hadd` histogram files from a directory, or to merge histograms with different selections (i.e., in different directories) are available in `scripts/`:
- `addHistosScouting.py`: to `hadd` histogram files from a directory;
- `mergeHistograms.py`: to merge histograms with different selections.

`bash` scripts to locally submit jobs are avaiable in the same directory.
E.g., to `hadd` histograms within a directory for all data sub-samples:
``` shell
source scripts/submitLocalAddHistosScouting.sh
```

## Draft analysis code with uproot and coffea

This is a draft of some potential analysis code, based on uproot and coffea.
To install on the uaf, run `source bootstrap.sh` (only required once).
Then, `./shell` will start the singularity container.

Inside `scouting/`, run `python minimal.py`.

