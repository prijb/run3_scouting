#!/bin/bash

export X509_USER_PROXY=$(voms-proxy-info -path)

usage()
{
    echo "Usage:"
    echo ""
    echo "  sh condor/runScoutingLooper_onCondor.sh ['notar'] input_dir"
    echo ""
    echo "The output directory will be created in /ceph/cms/store/user/$USER/Run3ScoutingOutput/"
    echo "Control the jobs to be run by editing runScoutingLooper_onCondor.sub"
    echo ""
    exit
}

if [ -z $1 ]; then usage; fi

notar=0
indir=""

if [ $# -gt 1 ]
then
    if [ "$1" == "notar" ]
    then
	notar=1
	indir=$2
    else
	indir=$1
    fi
fi

export SCOUTINGOUTPUTDIR=${indir}
export STARTDIR=$PWD

mkdir -p condor/plotting_logs
mkdir -p /ceph/cms/store/user/$USER/Run3ScoutingOutput/$SCOUTINGOUTPUTDIR

if [ ${notar} -gt 0 ]
then
    condor_submit condor/runScoutingLooper_onCondor.sub
else
    sh condor/create_package.sh
    condor_submit condor/runScoutingLooper_onCondor.sub    
fi
