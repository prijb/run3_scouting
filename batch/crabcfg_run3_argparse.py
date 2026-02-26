# crabcfg with argparse (fewer edge cases)
import argparse

parser = argparse.ArgumentParser(description="Generate and submit crab config")
parser.add_argument("--era", type=str, required=True, help="Era used (example: 2022C, 2022X, 2023C)")
parser.add_argument("--dataset", type=str, help="Dataset to use (overrides era based defintion)")
parser.add_argument("--data", action="store_true", help="Is data")
parser.add_argument("--PFMonitor", action="store_true", help="Dataset is PFMonitor (only for data)")
parser.add_argument("--triggerV10", action="store_true", help="triggerV10 is used (only for 2023)")
parser.add_argument("--syst", action="store_true", help="Introduce orthogonal bits and avoid skimming (for both data and MC systematics)")
parser.add_argument("--filebased", action="store_true", help="Change job submission type to be file based (debug)")
parser.add_argument("--dryrun", action="store_true", help="Don't submit")
args = parser.parse_args()
ntuple_version = "6"

from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
from CRABAPI.RawCommand import crabCommand

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
config = config()

import sys

era = args.era
data = args.data

year=0
if ("2022") in era:
    year=2022
elif ("2023") in era:
    year=2023
else:
    print(f"Era unspecified! Quitting...")
    quit()

# Make default config params
config.JobType.pyCfgParams=["era={}".format(era),"data={}".format(data),]
config.Data.inputDataset = '/ScoutingPFRun3/Run{}-v1/RAW'.format(era)
requestName = "skim4__{}".format(era)

## Modifications
# Replace by PFMonitor if called
if args.PFMonitor:
    config.Data.inputDataset = '/ScoutingPFMonitor/Run{}-v1/RAW'.format(era)
    config.JobType.pyCfgParams.append("monitor=True")
    requestName += "_PFMonitor"

# Add triggerv10 to era
if args.triggerV10:
    config.JobType.pyCfgParams[0] += "-triggerV10"
    requestName += "_triggerV10"

# Manually override dataset
if args.dataset is not None:
    config.Data.inputDataset = args.dataset
    #requestName += "_{}".format(args.dataset.split("/")[1])
    # BToJPsi dataset name's too long for the CRAB 100 char limit for request name :p
    requestName += "_{}".format(args.dataset.split("/")[1].split("_pythia8")[0]) 

# Add systematics option
if args.syst:
    config.JobType.pyCfgParams.append("syst=True")
    requestName += "_syst"

## Finish making the config
config.General.requestName = "{}_{}".format(
    requestName,
    ntuple_version,
)

import os
base = os.environ["CMSSW_BASE"]
config.General.workArea = base+'/..'

config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Scouting/NtupleMaker/test/producer_Run3.py'

config.Data.splitting = 'EventAwareLumiBased'

if (args.data):
    config.Data.unitsPerJob = int(10e6/3)
else:
    config.Data.unitsPerJob = int(10e4)

if (data and year==2023):
    config.Data.lumiMask = "data/Cert_Collisions2023_366442_370790_Golden.json"
    if(era=="2023C" and args.triggerV10):
        config.Data.lumiMask = "data/Cert_Collisions2023_eraC_367095_368823_Golden_1.json"
    elif(era=="2023C"):
        config.Data.lumiMask = "data/Cert_Collisions2023_eraC_367095_368823_Golden_2.json"

if (data and year==2022):
   config.Data.lumiMask = "data/Cert_Collisions2022_355100_362760_Golden.json"

# Changing to filebased
if args.filebased:
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = int(1)


#edit the area and user name
config.Data.outLFNDirBase = '/store/user/ppradeep/Run3Scouting/RAWScouting_'+ntuple_version # DB no
config.Data.publication = False
config.Site.storageSite = "T2_US_UCSD"

print(config)
if not args.dryrun:
    crabCommand('submit', config = config, dryrun = False) ## dryrun = True for local test