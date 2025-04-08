#!/bin/bash

export X509_USER_PROXY=$(voms-proxy-info -path)

usage()
{
    echo "Usage:"
    echo ""
    echo "  sh utils/condor_limits/runLimits_onCondor.sh input_dir output_dir period"
    echo ""
    echo "The output_dir will be created in /ceph/cms/store/user/$USER/ZPrimeSnTOutput/"
    echo "Control the jobs to be run by editing runLimits_onCondor.sub"
    echo ""
    exit
}

if [ -z $1 ]; then usage; fi

export SCOUTINGSNTINPUTDIRLIM=$1
export SCOUTINGSNTOUTPUTDIRLIM=$2
export PERIOD=$3
export NAME=$4
export HOMEDIR=$PWD
export LABEL=$(basename $SCOUTINGSNTOUTPUTDIRLIM)

echo "Creating dirs..."
mkdir -p condor/limits/limits_logs
mkdir -p /ceph/cms/store/user/$USER/Run3ScoutingOutput/$SCOUTINGSNTOUTPUTDIRLIM

echo "Preparing to create package..."
sh condor/limits/create_package.sh $SCOUTINGSNTINPUTDIRLIM
mv package.tar.gz package_${LABEL}.tar.gz

## Submission files to try extract limits on different subsets of signal samples
#condor_submit condor/limits/runLimits_onCondor.sub
#condor_submit condor/limits/runLimits_HTo2ZdTo2mu2x_onCondor.sub   # Mass grid for HTo2ZdTo2mu2x
#condor_submit condor/limits/runLimits_HTo2ZdTo2mu2x_ultrafine_onCondor.sub   # Mass grid for HTo2ZdTo2mu2x (fine)
#condor_submit condor/limits/runLimits_HTo2ZdTo2mu2x_toys_onCondor.sub   # Toys (reduced)
#condor_submit condor/limits/runLimits_BToPhi_onCondor.sub          # Mass grid for BToPhi
#condor_submit condor/limits/runLimits_BToPhi_SignalMC_onCondor.sub
#condor_submit condor/limits/runLimits_ScenarioA_SignalMC_onCondor.sub
#condor_submit condor/limits/runLimits_ScenarioB1_SignalMC_onCondor.sub
#condor_submit condor/limits/runLimits_reweighting_limited_onCondor.sub
condor_submit condor/limits/runLimits_reweighting_onCondor.sub
#condor_submit condor/limits/runLimits_reweighting_resubmit_onCondor.sub
#condor_submit condor/limits/runLimits_reweightingValidation_onCondor.sub
#condor_submit condor/limits/runLimits_limited_onCondor.sub
