from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
from CRABAPI.RawCommand import crabCommand

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
config = config()

import sys

era = "2023"+sys.argv[1] # ABCDEFG
extra = sys.argv[2]
ntuple_version = "2"

config.General.requestName = 'skim__{}{}_{}'.format(
        era,
        extra,
        ntuple_version,
        )

config.Data.inputDataset = '/ScoutingPFRun3/Run{}-v1/RAW'.format(era)

import os 
base = os.environ["CMSSW_BASE"]
config.General.workArea = base+'/..'

config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Scouting/NtupleMaker/test/producer_Run3.py'

if extra=="-triggerV10":
   era=era+extra

config.JobType.pyCfgParams=["era={}".format(era),"data=True",]

config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = int(10e6/2)

#something like this can be useful for limited disk availability
#config.Data.inputBlocks = [
#'/ScoutingPFRun3/Run2022A-v1/RAW#eb476276-22f6-47f0-8e0a-17245b460227',
#'/ScoutingPFRun3/Run2022A-v1/RAW#ff292f67-5ce4-4240-92da-949bda18ade9',
#'/ScoutingPFRun3/Run2022A-v1/RAW#ff2aeb80-96a6-4007-a019-ff75eed7a527',
#]

if(era=="2023C-triggerV10"): 
   config.Data.lumiMask = "data/Cert_Collisions2023_eraC_367095_368823_Golden_1.json"
elif(era=="2023C"): 
   config.Data.lumiMask = "data/Cert_Collisions2023_eraC_367095_368823_Golden_2.json"
elif(era=="2023B"): 
   config.Data.lumiMask = "data/Cert_Collisions2023_eraB_366403_367079_Golden.json"

#edit the area and user name
config.Data.outLFNDirBase = '/store/user/legianni/testRAWScouting_'+ntuple_version # DB no
config.Data.publication = False
config.Site.storageSite = "T2_US_UCSD"

print(config)
crabCommand('submit', config = config, dryrun = False) ## dryrun = True for local test
