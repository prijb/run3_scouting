SCRAMARCH=slc7_amd64_gcc10
CMSSWVERSION=CMSSW_12_6_0
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/$SCRAMARCH/cms/cmssw/$CMSSWVERSION/src ; eval `scramv1 runtime -sh` ; cd -

voms-proxy-init -voms cms --valid 168:00
