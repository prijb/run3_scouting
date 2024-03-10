#!/bin/bash

ulimit -s unlimited
source /cvmfs/cms.cern.ch/cmsset_default.sh

startTime=`date +%s`

OUTPUTDIR=$1
OUTPUTNAME=$2
IFILE=$3
NEVENTS=$4
CMSSWVERSION=$5
SCRAMARCH=$6
GENCFG=$7
ERA=$8

# Make sure OUTPUTNAME doesn't have .root since we add it manually
OUTPUTNAME=$(echo $OUTPUTNAME | sed 's/\.root//')

SCRAM_ARCH=${SCRAMARCH}

function getjobad {
    grep -i "^$1" "$_CONDOR_JOB_AD" | cut -d= -f2- | xargs echo
}

echo -e "\n--- begin header output ---\n" #                     <----- section division
echo "OUTPUTDIR: $OUTPUTDIR"
echo "OUTPUTNAME: $OUTPUTNAME"
echo "IFILE: $IFILE"
echo "CMSSWVERSION: $CMSSWVERSION"
echo "SCRAMARCH: $SCRAMARCH"
echo "GENCFG: $GENCFG"
echo "ERA: $ERA"

echo "GLIDEIN_CMSSite: $GLIDEIN_CMSSite"
echo "hostname: $(hostname)"
echo "uname -a: $(uname -a)"
echo "time: $(date +%s)"
echo "args: $@"
#echo "tag: $(getjobad tag)"
#echo "taskname: $(getjobad taskname)"

echo -e "\n--- end header output ---\n" #                       <----- section division

# NOTE removed chirp, should be added back in

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

function edit_psets {
    nevents=$1
    # gensim
    echo "process.source.firstLuminosityBlock = cms.untracked.uint32($IFILE)" >> $gensimcfg
    echo "process.RandomNumberGeneratorService.generator.initialSeed = $IFILE" >> $gensimcfg
    echo "process.maxEvents.input = $nevents" >> $gensimcfg
    echo "process.options.numberOfThreads=cms.untracked.uint32(1)" >> $gensimcfg
    echo "process.options.numberOfStreams=cms.untracked.uint32(0)" >> $gensimcfg
 
    # rawsim
    echo "process.maxEvents.input = $nevents" >> $rawsimcfg

    # aodsim
    echo "process.maxEvents.input = $nevents" >> $aodsimcfg

    ## slimmer
    ## NOTE: no slimmer defined yet
    #echo "process.maxEvents.input = $nevents" >> $slimmercfg

}

echo "time at start: $(date +%s)"

GENCFG=${GENCFG#*/}
gensimcfg=$GENCFG
rawsimcfg="rawsim_${ERA}_cfg.py"
aodsimcfg="aodsim_${ERA}_cfg.py"

echo "RAW-sim file: $rawsimcfg"
echo "AOD-sim file: $aodsimcfg"

echo "Number of events: $NEVENTS"
#echo "Fragment file: $FRAGMENT"
#FRAG=${FRAGMENT#*/}

echo "Unpacking shipped tarball"
tar -xf package.tar.gz
rm -f package.tar.gz

BASEDIR=`pwd`

echo "ls -la after unpacking:"
ls -la

echo "Setting up CMSSW"
cd /cvmfs/cms.cern.ch/$SCRAMARCH/cms/cmssw/$CMSSWVERSION/src ; eval `scramv1 runtime -sh`; eval `cmsenv` ; cd -

echo "[DIR STAT] ls -la after cmsrel"
ls -la

echo "Preparing to run"
# This is only going to work now, need to adjust random seed etc
# cmsRun darkshower_cfg.py  # for debugging only FIXME remove

#mkdir -p Configuration/GenProduction/python/
#cp ../../$FRAG Configuration/GenProduction/python/
#scram b
#
cp ${BASEDIR}/signal/$gensimcfg .
cp ${BASEDIR}/psets/$rawsimcfg .
cp ${BASEDIR}/psets/$aodsimcfg .

echo "Editing psets"
edit_psets $NEVENTS
#
## redoing this everytime because the fragment changes for different signal points.
## could think of a less resource intensive way in the future once we've really figured out
## the signal models...
#cmsDriver.py Configuration/GenProduction/python/$FRAG --customise_commands process.source.numberEventsInLuminosityBlock="cms.untracked.uint32($NEVENTS)" \
#  --python_filename gensim_cfg.py \
#  --eventcontent RAWSIM \
#  --customise Configuration/DataProcessing/Utils.addMonitoring \
#  --datatier GEN-SIM \
#  --fileout file:${OUTPUTNAME}_gensim.root\
#  --conditions 124X_mcRun3_2022_realistic_v12 \
#  --beamspot Realistic25ns13p6TeVEarly2022Collision \
#  --step GEN,SIM \
#  --geometry DB:Extended \
#  --era Run3 \
#  --no_exec \
#  --mc \
#  -n $NEVENTS

# setting the random seed for pythia only samples works like this apparently
# in contrast to setting the seed for the external LHE producer

echo "[DIR STAT] ls -la befor GENSIM running"
ls -la

echo "[RUN] gensim $gensimcfg"
cmsRun $gensimcfg

echo "[DIR STAT] ls -la after GENSIM running"
ls -la

echo "[RUN] rawsim $rawsimcfg"
cmsRun $rawsimcfg

echo "[DIR STAT] ls -la after RAWSIM running"
ls -la

echo "[RUN] aodsim $aodsimcfg"
cmsRun $aodsimcfg
CMSRUN_STATUS=$?

echo "[DIR STAT] ls -la after AODSIM running"
ls -la

if [[ $CMSRUN_STATUS != 0 ]]; then
    echo "Removing output file because cmsRun crashed with exit code $?"
    rm output_aod.root
    exit 1
fi

echo "Done with running, Stageout time"

echo "time before copy: $(date +%s)"

echo "Local output dir"
echo ${OUTPUTDIR}

REP="/store"
OUTPUTDIR="${OUTPUTDIR/\/ceph\/cms\/store/$REP}"

echo "Final output path for xrootd:"
echo ${OUTPUTDIR}

COPY_SRC="file://`pwd`/output_aod.root"
COPY_DEST=" davs://redirector.t2.ucsd.edu:1095/${OUTPUTDIR}/${OUTPUTNAME}_${IFILE}.root"
stageout $COPY_SRC $COPY_DEST


echo -e "\n--- end copying output ---\n" #                      <----- section division

echo "time at end: $(date +%s)"
