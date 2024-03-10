from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
from CRABAPI.RawCommand import crabCommand

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
#config = config()

import sys

era = sys.argv[1] # [2022, 2022postEE, 2023 or 2023BPix]
signal = sys.argv[2] # HTo2ZdTo2mu2x

year=0
if era=='2022':
    datatag='Run3Summer22'
elif era=='2022postEE':
    datatag='Run3Summer22EE'
else:
    quit()

# This is only MC, should not be used to run on data
data=False

# ntuple version defined now
ntuple_version = "s1p0"

# Setup working environment
import os
base = os.environ["CMSSW_BASE"]
#if not os.path.exists(base + '/src/centralTasks'):
#
#    os.mkdir(base + '/src/centralTasks')

if len(sys.argv) > 2:
    config = config()
    config.General.workArea = base+'/../'
    config.General.transferLogs = True
    config.JobType.pluginName = 'PrivateMC' # Analysis for other steps
    config.JobType.psetName = 'signal/{signal}_{era}_cfg.py'.format(signal=signal, era=era)
    config.Data.splitting = 'EventBased'
    config.Data.unitsPerJob = 1000
    config.Data.totalUnits = 1000000
    config.Data.outLFNDirBase = '/store/group/Run3Scouting/GENScouting_'+ntuple_version # DB no
    config.Data.publication = True
    config.Data.outputPrimaryDataset = '{signal}-pythia8'.format(signal=signal)
    config.Data.outputDatasetTag = 'private-{datatag}'.format(datatag=datatag)
    config.Site.storageSite = "T2_US_UCSD"
    config.General.requestName = 'crab_{signal}_{era}'.format(signal=signal, era=era)
    print(config)
    crabCommand('submit', config = config, dryrun = False) ## dryrun = True for local test
else:
    quit()



