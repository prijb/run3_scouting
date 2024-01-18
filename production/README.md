# Sample production

## Signal sample production

### Installing CMSSW with custom pythia version using ALMA-8 (recommended)
``` shell
cmssw-el8
cmsrel CMSSW_12_4_12
cd CMSSW_12_4_12/src/
cmsenv
cd ../..

mkdir pythia8
cd pythia8/
export MYPYTHIA=`pwd`
wget https://pythia.org/download/pythia83/pythia8309.tgz
tar -xf pythia8309.tgz
cd pythia8309

./configure --prefix=$MYPYTHIA --enable-shared --with-hepmc2=/cvmfs/cms.cern.ch/el8_amd64_gcc10/external/hepmc/2.06.10-2c9cd2f87f463ffd02a99f594368a84f --with-lhapdf6=/cvmfs/cms.cern.ch/el8_amd64_gcc10/external/lhapdf/6.4.0-ae790c99b90d02ddcd723a0f776517df

cd include/Pythia8Plugins/
rm JetMatching.h
wget http://amcatnlo.web.cern.ch/amcatnlo/JetMatching.h
cd ../..

make -j64
make install

sed -i "/<environment name=\"PYTHIA8_BASE\" default=/c\\ \ \ \ <environment name=\"PYTHIA8_BASE\" default=\"$MYPYTHIA\"/>" $CMSSW_BASE/config/toolbox/el8_amd64_gcc10/tools/selected/pythia8.xml

cd $CMSSW_BASE/src
scram setup pythia8
cmsenv

git cms-addpkg GeneratorInterface/Pythia8Interface
scram tool remove evtgen
cp ../../signal/remove_evtgen.patch .
git apply remove_evtgen.patch
scram b -j 12
scram tool info pythia8 # check whether it is correctly linked to the standalone pythia8
```


### Testing that the setup works, creating GENSIM

``` shell
mkdir -p Configuration/GenProduction/python/
cp ../../signal/*.py Configuration/GenProduction/python/
scram b

cmsDriver.py Configuration/GenProduction/python/darkshower_fragment.py --python_filename darkshower_cfg.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM --fileout file:output_gensim.root --conditions 124X_mcRun3_2022_realistic_v12 --beamspot Realistic25ns13p6TeVEarly2022Collision --customise_commands process.source.numberEventsInLuminosityBlock="cms.untracked.uint32(100)" --step GEN,SIM --geometry DB:Extended --era Run3 --no_exec --mc -n 5
```
Make sure you remove the following lines from the cmsrun cfg:

``` python
            'ParticleDecays:limitTau0 = on',
            'ParticleDecays:tau0Max = 10',
```

``` shell
cmsRun darkshower_cfg.py
```

### Running RAWSIM and AOD

psets for RAWSIM and AODSIM can be found in `psets/` and run in the same CMSSW release.
Please run CMSSW in the ALMA-8 container, `cmssw-el8`.


## Batch submission (WIP)

All of this is outside any container or CMS environment.
Make sure the `ProjectMetis` submodule is checked out.
Run `setup.sh` to add metis to your path.
Create the tarball to be shipped to the condor workers: `source make_tarball.sh`.
Run `python3 submitter.py`.


## Baby making
You should already have a CMSSW_12_4_12 release setup, using the cmssw-el8 container.
If not, please follow the first lines of instructions for the pythia installation above.

Very basic baby maker based on [Nick's version](https://github.com/aminnj/scouting/blob/master/batch/babymaker.py) ported to python3 and Run3.
This is still WIP, usage at own risk.
Example:
``` shell
cmssw-el8 --bind /ceph
cd CMSSW_12_4_12/src/; cmsenv; cd -
ipython -i babymaker.py /ceph/cms/store/user/isuarez/ProjectMetis/DarkShower_ScenarioA_default_Run3Summer22GS_v0p29_AODSIM_v0p29/output_82.root -- -o test.root -n -1 -y 2018
```


## Legacy instructions

### Installing Pythia8.309 locally outside of a container (not recommended)

``` shell
cmsrel CMSSW_12_4_12
cd CMSSW_12_4_12/src/
cmsenv

cd ../..

mkdir pythia8
cd pythia8/
export MYPYTHIA=`pwd`
wget https://pythia.org/download/pythia83/pythia8309.tgz
tar -xf pythia8309.tgz
cd pythia8309

./configure --prefix=$MYPYTHIA --enable-shared --with-hepmc2=/cvmfs/cms.cern.ch/slc7_amd64_gcc10/external/hepmc/2.06.10-llifpc --with-lhapdf6=/cvmfs/cms.cern.ch/slc7_amd64_gcc10/external/lhapdf/6.4.0-ae790c99b90d02ddcd723a0f776517df

cd include/Pythia8Plugins/
rm JetMatching.h
wget http://amcatnlo.web.cern.ch/amcatnlo/JetMatching.h
cd ../..

make -j64
make install

sed -i "/<environment name=\"PYTHIA8_BASE\" default=/c\\ \ \ \ <environment name=\"PYTHIA8_BASE\" default=\"$MYPYTHIA\"/>" $CMSSW_BASE/config/toolbox/slc7_amd64_gcc10/tools/selected/pythia8.xml

cd $CMSSW_BASE/src
scram setup pythia8
cmsenv

git cms-addpkg GeneratorInterface/Pythia8Interface
scram tool remove evtgen
cp ../../signal/remove_evtgen.patch .
git apply remove_evtgen.patch
scram b -j 12
scram tool info pythia8 # check whether it is correctly linked to the standalone pythia8
```

