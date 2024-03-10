#!/bin/bash

export X509_USER_PROXY=$(voms-proxy-info -path)

usage()
{
    echo "Usage:"
    echo ""
    echo "  sh condor/runScoutingLooper_onCondor.sh ['notar'] ['2023'] input_dir"
    echo ""
    echo "The output directory will be created in /ceph/cms/store/user/$USER/Run3ScoutingOutput/"
    echo "Control the jobs to be run by editing runScoutingLooper_onCondor.sub or the corresponding 2023 file"
    echo ""
    exit
}

if [ -z $3 ]; then usage; fi

export ERA=$2
export GEN=$3
export STARTDIR=$PWD
export GENFILE="${GEN}_${ERA}_cfg.py"

echo "Samples will be saved in: $MCOUTPUTDIR"
echo "Selected era: $ERA"
echo "Running for sample: $GEN"
echo "Configuration file: $GENFILE"


if [ $ERA == "2022" ]; then export CMSSWVER="CMSSW_12_4_16"; fi

mkdir -p condor/mc_logs
mkdir -p /ceph/cms/store/user/$USER/Run3ScoutingOutput/$1
export MCOUTPUTDIR="/ceph/cms/store/user/${USER}/Run3ScoutingOutput/${1}"


condor_submit condor/produce.sub

