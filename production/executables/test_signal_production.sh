#!/bin/bash

# This needs to run in the appropriate singularity container!

source /cvmfs/cms.cern.ch/cmsset_default.sh

startTime=`date +%s`

OUTPUTDIR=$1
OUTPUTNAME=$2
INPUTFILENAMES=$3
IFILE=$4
CMSSWVERSION=$5
SCRAMARCH=$6
GENCFG=$7

# Make sure OUTPUTNAME doesn't have .root since we add it manually
OUTPUTNAME=$(echo $OUTPUTNAME | sed 's/\.root//')

SCRAM_ARCH=${SCRAMARCH}

function getjobad {
    grep -i "^$1" "$_CONDOR_JOB_AD" | cut -d= -f2- | xargs echo
}

echo -e "\n--- begin header output ---\n" #                     <----- section division
echo "OUTPUTDIR: $OUTPUTDIR"
echo "OUTPUTNAME: $OUTPUTNAME"
echo "INPUTFILENAMES: $INPUTFILENAMES"
echo "IFILE: $IFILE"
echo "CMSSWVERSION: $CMSSWVERSION"
echo "SCRAMARCH: $SCRAMARCH"
echo "GENCFG: $GENCFG"

echo "GLIDEIN_CMSSite: $GLIDEIN_CMSSite"
echo "hostname: $(hostname)"
echo "uname -a: $(uname -a)"
echo "time: $(date +%s)"
echo "args: $@"
echo "tag: $(getjobad tag)"
echo "taskname: $(getjobad taskname)"

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

function setup_cmssw {
    echo "Setting up CMSSW $1 using SCRAM_ARCH $2"
    CMSSW=$1
    export SCRAM_ARCH=$2
    cmsrel $CMSSW
    cd $CMSSW
    eval $(scramv1 runtime -sh)
    cd -
}

echo "time at start: $(date +%s)"

GENCFG=${GENCFG#*/}
gensimcfg=$GENCFG
rawsimcfg="rawsim_cfg.py"
aodsimcfg="aodsim_cfg.py"

NEVENTS=$(getjobad param_nevents)
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

#setup_cmssw $CMSSWVERSION $SCRAM_ARCH
cmsrel $CMSSWVERSION
cd $CMSSWVERSION
CMSSW_BASE=`pwd`
echo "CMSSW_BASE: $CMSSW_BASE"
eval $(scramv1 runtime -sh)
cmsenv
cd -

echo "[DIR STAT] ls -la after cmsrel"
ls -la

echo "Preparing pythia"
mkdir pythia8
mv pythia8309.tgz pythia8
cd pythia8/
MYPYTHIA=`pwd`
tar -xf pythia8309.tgz
cd pythia8309

./configure --prefix=$MYPYTHIA --enable-shared --with-hepmc2=/cvmfs/cms.cern.ch/el8_amd64_gcc10/external/hepmc/2.06.10-2c9cd2f87f463ffd02a99f594368a84f --with-lhapdf6=/cvmfs/cms.cern.ch/el8_amd64_gcc10/external/lhapdf/6.4.0-ae790c99b90d02ddcd723a0f776517df

cd include/Pythia8Plugins/
rm JetMatching.h
wget http://amcatnlo.web.cern.ch/amcatnlo/JetMatching.h
cd ../..

echo "[DIR STAT] ls -la after setting up pythia"
ls -la

echo "Installing pythia"
make -j64
make install

sed -i "/<environment name=\"PYTHIA8_BASE\" default=/c\\ \ \ \ <environment name=\"PYTHIA8_BASE\" default=\"$MYPYTHIA\"/>" $CMSSW_BASE/config/toolbox/el8_amd64_gcc10/tools/selected/pythia8.xml

echo "CMSSW_BASE: $CMSSW_BASE"
cd $CMSSW_BASE/src
cmsenv
scram setup pythia8

echo "Removing evtgen"
cmsenv
git config --global user.name 'Indara Suarez'
git config --global user.github indarasuarez
git config --global user.email 'isuarez@cern.ch'
git cms-addpkg GeneratorInterface/Pythia8Interface
scram tool remove evtgen

cp ${BASEDIR}/signal/remove_evtgen.patch .
git apply remove_evtgen.patch

echo "Scram time"
scram b

echo "Preparing to run"
# This is only going to work now, need to adjust random seed etc
# cmsRun darkshower_cfg.py  # for debugging only FIXME remove

#mkdir -p Configuration/GenProduction/python/
#cp ../../$FRAG Configuration/GenProduction/python/
#scram b
#
cp ${BASEDIR}/psets/$gensimcfg .
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
