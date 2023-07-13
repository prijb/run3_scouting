import FWCore.ParameterSet.VarParsing as VarParsing

opts     = VarParsing.VarParsing('python')
vpbool   = VarParsing.VarParsing.varType.bool
vpint    = VarParsing.VarParsing.varType.int
vpstring = VarParsing.VarParsing.varType.string

opts.register('data',    True ,         mytype = vpbool)
opts.register('era',     "2022D",       mytype = vpstring)
opts.register('output',  "output.root", mytype = vpstring)
opts.register('inputs',  "",            mytype = vpstring) # comma separated list of input files
opts.register('nevents', -1,            mytype = vpint)
opts.parseArguments()

def convert_fname(fname):
    if fname.startswith("/store") or fname.startswith("root://"): return fname
    return "file:" + fname.replace("file:","",1)

inputs = []
if opts.inputs:
    inputs = map(convert_fname,opts.inputs.split(","))

lin=list(inputs)

print(""" Running with options:data = {} era = {} output = {} inputs = {} nevents = {} """.format(opts.data,opts.era,opts.output,lin,opts.nevents))

import FWCore.ParameterSet.Config as cms
process = cms.Process('SLIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


gtag=""
# from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
if opts.data:
    if '2022' in opts.era:
        gtag="124X_dataRun3_Prompt_v10" # latest prompt RECO GT
        #gtag="124X_dataRun3_HLT_v7" # latest HLT GT
    else:
        gtag="124X_dataRun3_Prompt_v4" # latest prompt RECO GT (July 2023, CMSSW>=13_0_5_patch2)
        #gtag="130X_dataRun3_HLT_v2" # latest HLT GT (July 2023)
else:
    if '2022' in opts.era:
        gtag="124X_mcRun3_2022_realistic_v11" # latest MC GT
    else:
        gtag="124X_mcRun3_2022_realistic_v11" # latest MC GT (same as for 2022, as of July 2023)
process.GlobalTag.globaltag = gtag

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(opts.nevents))

if opts.data:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
else:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(True),
    fileNames = cms.untracked.vstring(),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop *_hltScoutingTrackPacker_*_*',
    ),
)
process.source.fileNames = lin

if not opts.data:
    process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# skim data, but keep all events for MC acceptance calculations
do_skim = opts.data

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:output.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('skimpath'),
        ),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_hltScoutingMuonPacker_*_*",
        "keep *_hltScoutingPFPacker_*_*",
        "keep *_hltScoutingPrimaryVertexPacker_*_*",
        "keep *_triggerMaker_*_*",
        "keep *_hitMaker_*_*",
        "keep *_beamSpotMaker_*_*",
        "keep *_genParticles_*_HLT",
        ),
     basketSize = cms.untracked.int32(128*1024), # 128kb basket size instead of ~30kb default
)

if len(opts.output):
    process.out.fileName = "file:" + str(opts.output).strip().replace("file:","")

process.outpath = cms.EndPath(process.out)

process.Timing = cms.Service("Timing",
        summaryOnly = cms.untracked.bool(True)
        )

process.countmu = cms.EDFilter("ScoutingMuonCountFilter",
    src = cms.InputTag("hltScoutingMuonPacker"),
    minNumber = cms.uint32(2)
)

process.countvtx = cms.EDFilter("ScoutingVertexCountFilter",
    src = cms.InputTag("hltScoutingMuonPacker","displacedVtx"),
    minNumber = cms.uint32(1)
)

#from Daniel's link - https://github.com/elfontan/Run3ScoutingAnalysisTools/blob/ScoutingPaper/tree.py#L38
L1Info  = []
HLTInfo = []
if '2022' in opts.era or '2023B' in opts.era or '2023C-triggerV10' in opts.era:
    L1Info = ["L1_DoubleMu_12_5","L1_DoubleMu_15_7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18","L1_DoubleMu4_SQ_OS_dR_Max1p2","L1_DoubleMu4p5_SQ_OS_dR_Max1p2"]
    HLTInfo=[ [ 'Run3_PFScouting', 'DST_Run3_PFScoutingPixelTracking_v*' ] ]
else:
    # for run >=367621 (during era Run2023C)
    L1Info = ["L1_DoubleMu_12_5","L1_DoubleMu_15_7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18","L1_DoubleMu4_SQ_OS_dR_Max1p2","L1_DoubleMu4p5_SQ_OS_dR_Max1p2","L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4","L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4","L1_DoubleMu8_SQ"]
    HLTInfo = [ [ 'Run3_DoubleMu3_PFScouting', 'DST_Run3_DoubleMu3_PFScoutingPixelTracking_v*' ] ]

do_trigger_objects = not opts.data

process.triggerMaker = cms.EDProducer("TriggerMaker",
        triggerAlias = cms.vstring(list(zip(*HLTInfo))[0]),
        triggerSelection = cms.vstring(list(zip(*HLTInfo))[1]),
        triggerConfiguration = cms.PSet(
            hltResults            = cms.InputTag('TriggerResults','','HLT'),
            l1tResults            = cms.InputTag(''),
            l1tIgnoreMaskAndPrescale = cms.bool(False),
            throw                 = cms.bool(True),
            usePathStatus = cms.bool(False),
            ),
        doL1 = cms.bool(True),
        doTriggerObjects = cms.bool(do_trigger_objects),
        AlgInputTag = cms.InputTag("gtStage2Digis"),
        l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
        l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
        ReadPrescalesFromFile = cms.bool(False),
        l1Seeds = cms.vstring(L1Info),
        )

process.hitMaker = cms.EDProducer("HitMaker",
        muonInputTag = cms.InputTag("hltScoutingMuonPacker"),
        dvInputTag = cms.InputTag("hltScoutingMuonPacker:displacedVtx"),
        measurementTrackerEventInputTag = cms.InputTag("MeasurementTrackerEvent"),
        )

process.beamSpotMaker = cms.EDProducer("BeamSpotMaker")

from RecoTracker.MeasurementDet.measurementTrackerEventDefault_cfi import measurementTrackerEventDefault as _measurementTrackerEventDefault
process.MeasurementTrackerEvent = _measurementTrackerEventDefault.clone()

process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )

process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")
#process.patTrigger.triggerResults = cms.InputTag("TriggerResults","","HLT")
#process.patTrigger.triggerEvent = cms.InputTag("hltTriggerSummaryAOD","","HLT")
process.patTrigger.stageL1Trigger = cms.uint32(2)

#hitMaker not needed (Mario)
if (opts.data): process.skimpath = cms.Path(process.countmu+process.countvtx+process.gtStage2Digis+process.triggerMaker+process.offlineBeamSpot+process.beamSpotMaker)#+process.MeasurementTrackerEvent+process.hitMaker)
else: process.skimpath = cms.Path(process.gtStage2Digis+process.patTrigger+process.triggerMaker+process.offlineBeamSpot+process.beamSpotMaker)#+process.MeasurementTrackerEvent+process.hitMaker)

