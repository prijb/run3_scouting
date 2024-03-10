from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
from CRABAPI.RawCommand import crabCommand

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
#config = config()

import sys

era = sys.argv[1] # [2022, 2022postEE, 2023 or 2023BPix]

year=0
if era=='2022':
    datatag='Run3Summer22-RAWSIM'
    points = []
    points.append(['BToPhi_MPhi-2p0_ctau-1mm', '/BToPhi_MPhi-2p0_ctau-1mm-pythia8/Run3Scouting-private-Run3Summer22-2bf55cfd546ea52b41ec0931f75ff588/USER'])
    points.append(['BToPhi_MPhi-2p0_ctau-10mm', '/BToPhi_MPhi-2p0_ctau-10mm-pythia8/Run3Scouting-private-Run3Summer22-214ff9e364fb8e0196c805684eeed6e1/USER'])
    points.append(['BToPhi_MPhi-2p0_ctau-100mm', '/BToPhi_MPhi-2p0_ctau-100mm-pythia8/Run3Scouting-private-Run3Summer22-4ae7615fe22e489e3d317c7c95439c15/USER'])
elif era=='2022postEE':
    datatag='Run3Summer22EE-RAWSIM'
    points = []
    points.append(['BToPhi_MPhi-2p0_ctau-1mm', '/BToPhi_MPhi-2p0_ctau-1mm-pythia8/Run3Scouting-private-Run3Summer22EE-ef24b12ae71cb5e29170ffcadda91cc9/USER'])
    points.append(['BToPhi_MPhi-2p0_ctau-10mm', '/BToPhi_MPhi-2p0_ctau-10mm-pythia8/Run3Scouting-private-Run3Summer22EE-8a8b8da91722241b3460388233bab23e/USER'])
    points.append(['BToPhi_MPhi-2p0_ctau-100mm', '/BToPhi_MPhi-2p0_ctau-100mm-pythia8/Run3Scouting-private-Run3Summer22EE-6a64a9c756367b13bfe7182069417686/USER'])
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

if len(sys.argv) > 1:
    configs = []
    for point in points:
        signal = point[0]
        dataset = point[1]
        configs.append(config())
        #config = config()
        configs[-1].General.workArea = base+'/../'
        configs[-1].General.transferLogs = True
        configs[-1].JobType.pluginName = 'Analysis' # Analysis for other steps
        configs[-1].JobType.psetName = 'psets/rawsim_{era}_cfg.py'.format(signal=signal, era=era)
        configs[-1].JobType.maxMemoryMB = 4000
        configs[-1].Data.inputDataset = dataset
        configs[-1].Data.inputDBS = 'phys03'
        configs[-1].Data.splitting = 'FileBased'
        configs[-1].Data.unitsPerJob = 1
        configs[-1].Data.outLFNDirBase = '/store/group/Run3Scouting/GENScouting_'+ntuple_version # DB no
        configs[-1].Data.publication = True
        configs[-1].Data.outputDatasetTag = 'private-{datatag}'.format(datatag=datatag)
        configs[-1].Site.storageSite = "T2_US_UCSD"
        configs[-1].Site.whitelist = ['T2_US_UCSD', 'T2_ES_CIEMAT', 'T2_CH_CERN', 'T2_IT_Bari', 'T2_ES_IFCA', 'T2_US_Nebraska']
        configs[-1].General.requestName = 'rawsim_{signal}_{era}'.format(signal=signal, era=era)
        print(configs[-1])
        crabCommand('submit', config = configs[-1], dryrun = False) ## dryrun = True for local test
else:
    quit()



