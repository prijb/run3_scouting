# Sample production

## Signal sample production

### Installing Pythia8.309

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
patch remove_evtgen.patch
scram b -j 12
scram tool info pythia8 # check whether it is correctly linked to the standalone pythia8
```

### Testing that the setup works

``` shell
mkdir -p Configuration/GenProduction/python/
cp ../../signal/darkshower_fragment.py Configuration/GenProduction/python/
scram b

cmsDriver.py Configuration/GenProduction/python/darkshower_fragment.py --python_filename darkshower_cfg.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM --fileout file:darkshower.root --conditions 124X_mcRun3_2022_realistic_v12 --beamspot Realistic25ns13p6TeVEarly2022Collision --customise_commands process.source.numberEventsInLuminosityBlock="cms.untracked.uint32(100)" --step GEN,SIM --geometry DB:Extended --era Run3 --no_exec --mc -n 5

cmsRun darkshower_cfg.py
```

