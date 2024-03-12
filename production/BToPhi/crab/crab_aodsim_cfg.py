from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
from CRABAPI.RawCommand import crabCommand

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
#config = config()

import sys

era = sys.argv[1] # [2022, 2022postEE, 2023 or 2023BPix]

year=0
if era=='2022':
    datatag='Run3Summer22-AODSIM'
    points = []
    points.append(['BToPhi_MPhi-2p0_ctau-1mm', '/BToPhi_MPhi-2p0_ctau-1mm-pythia8/Run3Scouting-private-Run3Summer22-RAWSIM-905680e847f46ad1a422d47384efadca/USER'])
    points.append(['BToPhi_MPhi-2p0_ctau-10mm', '/BToPhi_MPhi-2p0_ctau-10mm-pythia8/Run3Scouting-private-Run3Summer22-RAWSIM-905680e847f46ad1a422d47384efadca/USER'])
    points.append(['BToPhi_MPhi-2p0_ctau-100mm', '/BToPhi_MPhi-2p0_ctau-100mm-pythia8/Run3Scouting-private-Run3Summer22-RAWSIM-905680e847f46ad1a422d47384efadca/USER'])
elif era=='2022postEE':
    datatag='Run3Summer22EE-AODSIM'
    points = []
    points.append(['BToPhi_MPhi-2p0_ctau-1mm', '/BToPhi_MPhi-2p0_ctau-1mm-pythia8/Run3Scouting-private-Run3Summer22EE-RAWSIM-87e245638d0e8549d13398f07ce44cec/USER'])
    points.append(['BToPhi_MPhi-2p0_ctau-10mm', '/BToPhi_MPhi-2p0_ctau-10mm-pythia8/Run3Scouting-private-Run3Summer22EE-RAWSIM-87e245638d0e8549d13398f07ce44cec/USER'])
    points.append(['BToPhi_MPhi-2p0_ctau-100mm', '/BToPhi_MPhi-2p0_ctau-100mm-pythia8/Run3Scouting-private-Run3Summer22EE-RAWSIM-87e245638d0e8549d13398f07ce44cec/USER'])
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
        configs[-1].JobType.psetName = 'psets/aodsim_{era}_cfg.py'.format(signal=signal, era=era)
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
        configs[-1].General.requestName = 'aodsim_{signal}_{era}'.format(signal=signal, era=era)
        print(configs[-1])
        crabCommand('submit', config = configs[-1], dryrun = False) ## dryrun = True for local test
else:
    quit()



