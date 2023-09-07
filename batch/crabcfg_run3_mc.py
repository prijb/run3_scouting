from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
from CRABAPI.RawCommand import crabCommand

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
config = config()

era = "2022" 
ntuple_version = "4"

config.General.requestName = 'skim__{}_{}'.format(
        era,
        "",
        )

import os 
base = os.environ["CMSSW_BASE"]
config.General.workArea = base+'/..'

config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Scouting/NtupleMaker/test/producer_Run3.py'

config.JobType.pyCfgParams=["era={}".format(era),"data=False",]

import sys
basedir = sys.argv[1]

#prepare local file list / it works when launching from ucsd
files=os.listdir(basedir)

for i in range(len(files)):
  files[i]=basedir+"/"+files[i]
for i in range(len(files)):
  files[i]=files[i].replace("/ceph/cms", "")

config.Data.userInputFiles = files 

config.Data.outputPrimaryDataset =basedir.split("/")[-1]
print(config.Data.outputPrimaryDataset)

config.General.requestName+=config.Data.outputPrimaryDataset

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 #can be larger

#edit the area and user name
config.Data.outLFNDirBase = '/store/group/Run3Scouting/RAWScouting_'+ntuple_version # DB no
config.Data.publication = False
config.Site.storageSite = "T2_US_UCSD"

print(config)
crabCommand('submit', config = config, dryrun = False) ## dryrun = True for local test
