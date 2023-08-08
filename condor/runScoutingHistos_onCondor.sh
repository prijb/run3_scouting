#!/bin/bash

export X509_USER_PROXY=$(voms-proxy-info -path)

usage()
{
    echo "Usage:"
    echo ""
    echo "  sh condor/runScoutingHistos_onCondor.sh ['notar'] input_dir output_dir"
    echo ""
    echo "The output_dir will be created in /ceph/cms/store/user/$USER/Run3ScoutingOutput/"
    echo "Control the jobs to be run by editing runScoutingHistos_onCondor.sub"
    echo ""
    exit
}

if [ -z $2 ]; then usage; fi

notar=0
indir=""
outdir=""
extraflags=""

if [ $# -gt 2 ]
then
    if [ "$1" == "notar" ]
    then
	notar=1
	indir=$2
	outdir=$3
    else
	indir=$1
	outdir=$2
    fi
fi

offset=$(($notar+2))
i=1
for arg in "$@"
do
    if [ $i -gt ${offset} ]
    then
        extraflags="$extraflags $arg"
    fi
    i=$((i + 1))
done

export SCOUTINGINPUTDIR=${indir}
export SCOUTINGOUTPUTDIR=${outdir}
export SCOUTINGARGS=${extraflags}
export STARTDIR=$PWD

mkdir -p condor/plotting_logs
mkdir -p /ceph/cms/store/user/$USER/Run3ScoutingOutput/$SCOUTINGOUTPUTDIR

if [ ${notar} -gt 0 ]
then
    condor_submit condor/runScoutingHistos_onCondor.sub
else
    sh condor/create_package.sh
    condor_submit condor/runScoutingHistos_onCondor.sub    
fi
