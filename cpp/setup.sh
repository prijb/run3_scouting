SCRAMARCH=slc7_amd64_gcc10
CMSSWVERSION=CMSSW_12_6_0
if [ $# -gt 0 ]
then
    if [ $1=="EL8" ]
    then
	SCRAMARCH="el8_amd64_gcc10"
	CMSSWVERSION=CMSSW_12_6_0
    fi
fi

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/$SCRAMARCH/cms/cmssw/$CMSSWVERSION/src ; eval `scramv1 runtime -sh` ; cd -

voms-proxy-init -voms cms --valid 168:00

alias testFiller="python3 fillHistosScouting.py --inDir /store/user/fernance/Run3ScoutingOutput/looperOutput_2022_Feb-05-2024 --inSample Signal_HTo2ZdTo2mu2x_MZd-8p0_ctau-10mm_2022 --outSuffix test"
