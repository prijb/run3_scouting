#!/bin/bash

# This needs to run in the appropriate singularity container!

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc10

startTime=`date +%s`

OUTPUTDIR=$1
OUTPUTNAME=$2
INPUTFILENAMES=$3
IFILE=$4
CMSSWVERSION=$5
SCRAMARCH=$6
FRAGMENT=$7

# Make sure OUTPUTNAME doesn't have .root since we add it manually
OUTPUTNAME=$(echo $OUTPUTNAME | sed 's/\.root//')

export SCRAM_ARCH=${SCRAMARCH}

function getjobad {
    grep -i "^$1" "$_CONDOR_JOB_AD" | cut -d= -f2- | xargs echo
}

function setup_chirp {
    if [ -e ./condor_chirp ]; then
    # Note, in the home directory
        mkdir chirpdir
        mv condor_chirp chirpdir/
        export PATH="$PATH:$(pwd)/chirpdir"
        echo "[chirp] Found and put condor_chirp into $(pwd)/chirpdir"
    elif [ -e /usr/libexec/condor/condor_chirp ]; then
        export PATH="$PATH:/usr/libexec/condor"
        echo "[chirp] Found condor_chirp in /usr/libexec/condor"
    else
        echo "[chirp] No condor_chirp :("
    fi
}

function chirp {
    # Note, $1 (the classad name) must start with Chirp
    condor_chirp set_job_attr_delayed $1 $2
    ret=$?
    echo "[chirp] Chirped $1 => $2 with exit code $ret"
}

function stageout {
    COPY_SRC=$1
    COPY_DEST=$2
    retries=0
    COPY_STATUS=1
    until [ $retries -ge 3 ]
    do
        echo "Stageout attempt $((retries+1)): env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-copy -p -f -t 7200 --verbose --checksum ADLER32 ${COPY_SRC} ${COPY_DEST}"
        env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-copy -p -f -t 7200 --verbose --checksum ADLER32 ${COPY_SRC} ${COPY_DEST}
        COPY_STATUS=$?
        if [ $COPY_STATUS -ne 0 ]; then
            echo "Failed stageout attempt $((retries+1))"
        else
            echo "Successful stageout with $retries retries"
            break
        fi
        retries=$[$retries+1]
        echo "Sleeping for 30m"
        sleep 30m
    done
    if [ $COPY_STATUS -ne 0 ]; then
        echo "Removing output file because gfal-copy crashed with code $COPY_STATUS"
        env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-rm --verbose ${COPY_DEST}
        REMOVE_STATUS=$?
        if [ $REMOVE_STATUS -ne 0 ]; then
            echo "Uhh, gfal-copy crashed and then the gfal-rm also crashed with code $REMOVE_STATUS"
            echo "You probably have a corrupt file sitting on ceph now."
            exit 1
        fi
    fi
}

function setup_environment {
    if [ -r "$OSGVO_CMSSW_Path"/cmsset_default.sh ]; then
        echo "sourcing environment: source $OSGVO_CMSSW_Path/cmsset_default.sh"
        source "$OSGVO_CMSSW_Path"/cmsset_default.sh
    elif [ -r "$OSG_APP"/cmssoft/cms/cmsset_default.sh ]; then
        echo "sourcing environment: source $OSG_APP/cmssoft/cms/cmsset_default.sh"
        source "$OSG_APP"/cmssoft/cms/cmsset_default.sh
    elif [ -r /cvmfs/cms.cern.ch/cmsset_default.sh ]; then
        echo "sourcing environment: source /cvmfs/cms.cern.ch/cmsset_default.sh"
        source /cvmfs/cms.cern.ch/cmsset_default.sh
    else
        echo "ERROR! Couldn't find $OSGVO_CMSSW_Path/cmsset_default.sh or /cvmfs/cms.cern.ch/cmsset_default.sh or $OSG_APP/cmssoft/cms/cmsset_default.sh"
        exit 1
    fi
}

echo "Fragment file: $FRAGMENT"
export FRAG=${FRAGMENT#*/}

NEVENTS=$(getjobad param_nevents)
echo "Number of events: $NEVENTS"

echo "Unpacking shipped tarball"
tar -xf tarball.tar.gz
rm -f tarball.tar.gz

echo "ls -la after unpacking:"
ls -la

echo "Setting environment variables"
export BASEDIR=`pwd`
export MYPYTHIA=`pwd`/pythia8/
echo "BASEDIR: $BASEDIR"
echo "MYPYTHIA: $MYPYTHIA"

echo "Setting up CMSSW"
cd CMSSW_12_4_12/src/
echo "scram ProjectRename"
scramv1 b ProjectRename
echo "cmsenv"
cmsenv

echo "Replacing pythia path"
sed -i "/<environment name=\"PYTHIA8_BASE\" default=/c\\ \ \ \ <environment name=\"PYTHIA8_BASE\" default=\"$MYPYTHIA\"/>" $CMSSW_BASE/config/toolbox/slc7_amd64_gcc10/tools/selected/pythia8.xml

echo "Compile (maybe need to build clean?)"
scram b

echo "Running"
# This is only going to work now, need to adjust random seed etc
cmsRun darkshower_cfg.py

# mkdir -p Configuration/GenProduction/python/
# cp $FRAG Configuration/GenProduction/python/
# scram b
cmsDriver.py Configuration/GenProduction/python/$FRAG --customise_commands process.source.numberEventsInLuminosityBlock="cms.untracked.uint32($NEVENTS)" \
  --python_filename darkshower_cfg.py \
  --eventcontent RAWSIM \
  --customise Configuration/DataProcessing/Utils.addMonitoring \
  --datatier GEN-SIM \
  --fileout file: $OUTPUTNAME.root\
  --conditions 124X_mcRun3_2022_realistic_v12 \
  --beamspot Realistic25ns13p6TeVEarly2022Collision \
  --step GEN,SIM \
  --geometry DB:Extended \
  --era Run3 \
  --no_exec \
  --mc \
  -n $NEVENTS

echo "process.source.firstLuminosityBlock = cms.untracked.uint32($IFILE)" >> darkshower_cfg.py  # additional customize_commands

echo "ls -la after running"
ls -la


# FIXME stage out still missing
