# To run the looper

Under the `cpp` folder, run:

1. `source setup.sh`: This sets up CMSSW (`12_6_0` for simplicity is used but the [differences wrt. to `12_6_0_patch1`](https://github.com/cms-sw/cmssw/compare/CMSSW_12_6_0...CMSSW_12_6_0_patch1) are not affecting the libraries used in the looper Makefile), and the proxy (hence the password prompt).
2. `make rootDict`: This sets up the proper ROOT dictionaries for output classes. Usually, it only needs to be run once. It must be remade only if a variable of a non-standard type is added to the output.
3. `make`: Compiles all of the relevant source code. It needs to be rerun each time a `.cc` file is modified.
4. `./main.exe /path/to/output/dir YEAR SAMPLE NUMBER_OF_EVENT_TO_START_RUNNING_FROM NUMBER_OF_EVENTS_TO_RUN (INCONDOR) (FROMCRAB)`: Runs the looper. Examples of arguments can be founf in the `condor/runScoutingLooper_onCondor.sub` file and details can be understood from the `main.cc` file.

## Examples

To read a sample locally e.g. HTo2ZdTo2mu2x:
```
./main.exe test 2022 Signal_HTo2ZdTo2mu2x_MZd-24p0_ctau-1mm_2022 0 100 0 1
```
just make sure to match the sample name and give the correct sample location in ```cpp/input/centralDatasets.txt```. 

To read a dataset from crab e.g. the BToPhi sample:
```
./main.exe test 2022 Signal_BToPhi_MPhi-2p0_ctau-100mm_2022postEE 0 100 0 1
```
The last 1 indicates that the input should be accessed by crab using Xrootd Service (AAA) for Remote Data Access (documented [here](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookXrootdService#ReDirector)). Just make sure to match the sample name and give the dataset name in ```cpp/input/centralDatasets.txt```.

**Note: By default it will look for privately produced datasets in 'phys03' dbs, which makes sense because the output of the skimming is private. No central samples should be given to the looper.**
