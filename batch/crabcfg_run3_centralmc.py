from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
from CRABAPI.RawCommand import crabCommand

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
#config = config()

import sys

era = sys.argv[1] # [2022, 2022postEE, 2023 or 2023BPix]
signal = sys.argv[2] # HTo2ZdTo2mu2x

year=0
if ("2022") in era:
    year=2022
elif ("2023") in era:
    year=2023
else:
    quit()

# This is only MC, should not be used to run on data
data=False

# ntuple version defined now
ntuple_version = "5c"

# Setup working environment
import os
base = os.environ["CMSSW_BASE"]
#if not os.path.exists(base + '/src/centralTasks'):
#
#    os.mkdir(base + '/src/centralTasks')

# List of datasets
config_list = []
extras = []
if (len(sys.argv)>2):
    mass_points = []
    if "HTo2ZdTo2mu2x" in sys.argv[2]:
        config = config()
        config.General.workArea = base+'/../'
        config.General.transferLogs = True
        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = 'Scouting/NtupleMaker/test/producer_Run3.py'
        config.Data.splitting = 'EventAwareLumiBased'
        config.Data.unitsPerJob = int(10e4)
        config.Data.outLFNDirBase = '/store/group/Run3Scouting/RAWScouting_'+ntuple_version # DB no
        config.Data.publication = False
        config.Site.storageSite = "T2_US_UCSD"
        # Set the points to produce
        mass_points.append(['0p5', '1'])
        mass_points.append(['0p5', '10'])
        mass_points.append(['0p7', '1'])
        mass_points.append(['0p7', '10'])
        mass_points.append(['1p5', '1'])
        mass_points.append(['1p5', '10'])
        mass_points.append(['1p5', '100'])
        mass_points.append(['2p0', '1'])
        mass_points.append(['2p0', '10'])
        mass_points.append(['2p0', '100'])
        mass_points.append(['2p5', '1'])
        mass_points.append(['2p5', '10'])
        mass_points.append(['2p5', '100'])
        mass_points.append(['3p0', '1'])
        mass_points.append(['3p0', '10'])
        mass_points.append(['3p0', '100'])
        mass_points.append(['4p0', '1'])
        mass_points.append(['4p0', '10'])
        mass_points.append(['4p0', '100'])
        mass_points.append(['5p0', '1'])
        mass_points.append(['5p0', '10'])
        mass_points.append(['5p0', '100'])
        mass_points.append(['6p0', '1'])
        mass_points.append(['6p0', '10'])
        mass_points.append(['6p0', '100'])
        mass_points.append(['7p0', '1'])
        mass_points.append(['7p0', '10'])
        mass_points.append(['7p0', '100'])
        mass_points.append(['8p0', '1'])
        mass_points.append(['8p0', '10'])
        mass_points.append(['8p0', '100'])
        mass_points.append(['10p0', '1'])
        mass_points.append(['10p0', '10'])
        mass_points.append(['10p0', '100'])
        mass_points.append(['12p0', '1'])
        mass_points.append(['12p0', '10'])
        mass_points.append(['12p0', '100'])
        mass_points.append(['14p0', '1'])
        mass_points.append(['14p0', '10'])
        mass_points.append(['14p0', '100'])
        mass_points.append(['16p0', '1'])
        mass_points.append(['16p0', '10'])
        mass_points.append(['16p0', '100'])
        mass_points.append(['20p0', '1'])
        mass_points.append(['20p0', '10'])
        mass_points.append(['20p0', '100'])
        mass_points.append(['22p0', '1'])
        mass_points.append(['22p0', '10'])
        mass_points.append(['22p0', '100'])
        mass_points.append(['24p0', '1'])
        mass_points.append(['24p0', '10'])
        mass_points.append(['24p0', '100'])
        mass_points.append(['30p0', '1'])
        mass_points.append(['30p0', '10'])
        mass_points.append(['30p0', '100'])
        mass_points.append(['30p0', '1000'])
        mass_points.append(['34p0', '1'])
        mass_points.append(['34p0', '10'])
        mass_points.append(['34p0', '100'])
        mass_points.append(['34p0', '1000'])
        mass_points.append(['40p0', '1'])
        mass_points.append(['40p0', '10'])
        mass_points.append(['40p0', '100'])
        mass_points.append(['40p0', '1000'])
        mass_points.append(['44p0', '1'])
        mass_points.append(['44p0', '10'])
        mass_points.append(['44p0', '100'])
        mass_points.append(['44p0', '1000'])
        mass_points.append(['50p0', '1'])
        mass_points.append(['50p0', '10'])
        mass_points.append(['50p0', '100'])
        mass_points.append(['50p0', '1000'])
        for [m,t] in mass_points:
            config_list.append(config)
            config_list[-1].JobType.pyCfgParams=["era={}".format(era),"data=False",]
            if era=="2022":
                dataset_name = '/HTo2ZdTo2mu2x_MZd-{}_ctau-{}mm_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM'.format(m, t)
                config_list[-1].Data.inputDataset = dataset_name
            elif era=="2022postEE":
                dataset_name = '/HTo2ZdTo2mu2x_MZd-{}_ctau-{}mm_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM'.format(m, t)
                config_list[-1].Data.inputDataset = dataset_name
            elif era=="2023":
                dataset_name = '/HTo2ZdTo2mu2x_MZd-{}_ctau-{}mm_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23DRPremix-130X_mcRun3_2023_realistic_v14-v1/AODSIM'.format(m, t)
                config_list[-1].Data.inputDataset = dataset_name
            elif era=="2023BPix":
                dataset_name = '/HTo2ZdTo2mu2x_MZd-{}_ctau-{}mm_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v1/AODSIM'.format(m, t)
                config_list[-1].Data.inputDataset = dataset_name
            config_list[-1].General.requestName = 'centralSkim__{}_{}_m-{}_ctau-{}mm_{}'.format(signal, era, m, t, ntuple_version)
            print(config)
            crabCommand('submit', config = config, dryrun = False) ## dryrun = True for local test
    #elif "[signal]" in sys.argv[2]: (<--- Add additional signals here)
    else:
        quit()
    
else:
    quit()



