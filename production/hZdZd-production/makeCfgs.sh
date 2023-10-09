#!/bin/bash

export X509_USER_PROXY=$(voms-proxy-info -path)
export SCRAM_ARCH=el8_amd64_gcc10
export FRAGMENTS=$1
WORK=$PWD

source /cvmfs/cms.cern.ch/cmsset_default.sh
if [ -r CMSSW_12_4_11_patch3/src ] ; then
  echo release CMSSW_12_4_11_patch3 already exists
else
  scram p CMSSW CMSSW_12_4_11_patch3
fi
cd CMSSW_12_4_11_patch3/src
eval `scram runtime -sh`
mkdir -p Configuration/GenProduction/python
cd ../../
cp $FRAGMENTS/*.py CMSSW_12_4_11_patch3/src/Configuration/GenProduction/python/
cd CMSSW_12_4_11_patch3/src
scram b
cmsenv

EVENTS=10
SEED=$(($(date +%s) % 100 + 1));

for FILE in Configuration/GenProduction/python/*
do
    echo $FILE;
    if [ $FILE == "Configuration/GenProduction/python/__init__.py" ]; then
	echo "Continueeee"
        continue;
    fi
    prefix=${FILE%.*};
    suffix="_cfg.py"
    CFG="$prefix$suffix"; 
    ROOT="file:$prefix.root"
    cmsDriver.py $FILE --python_filename $CFG --eventcontent RAWSIM,LHE --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM,LHE --fileout ROOT --conditions 124X_mcRun3_2022_realistic_postEE_v1 --beamspot Realistic25ns13p6TeVEarly2022Collision --customise_commands process.RandomNumberGeneratorService.externalLHEProducer.initialSeed="int(${SEED})" --step LHE,GEN,SIM --geometry DB:Extended --era Run3 --no_exec --mc -n $EVENTS;
    mv $FILE $WORK
done
