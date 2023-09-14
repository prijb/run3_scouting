#!/bin/bash

export X509_USER_PROXY=$(voms-proxy-info -path)

usage()
{
    echo "Usage:"
    echo ""
    echo "  sh condor/runScoutingHistos_onCondor.sh input_dir output_dir"
    echo ""
    echo "The output_dir will be created in /ceph/cms/store/user/$USER/Run3ScoutingOutput/"
    echo "Control the jobs to be run by editing runScoutingHistos_onCondor.sub"
    echo ""
    exit
}

if [ -z $1 ]; then usage; fi

str=" "
i=1;
for arg in "$@"
do
    if [ $i -gt 2 ]; then
        str="$str $arg";
    fi
    i=$((i + 1));
done

export SCOUTINGINPUTDIR=$1
export SCOUTINGOUTPUTDIR=$2
export SCOUTINGARGS=$str
export STARTDIR=$PWD

mkdir -p condor/plotting_logs
mkdir -p /ceph/cms/store/user/$USER/Run3ScoutingOutput/$SCOUTINGOUTPUTDIR

sh condor/create_package.sh

condor_submit condor/runScoutingHistos_onCondor.sub
