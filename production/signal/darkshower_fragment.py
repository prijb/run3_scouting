
import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunesRun3ECM13p6TeV.PythiaCP5Settings_cfi import *

generator = cms.EDFilter("Pythia8ConcurrentGeneratorFilter",
                         maxEventsToPrint = cms.untracked.int32(1),
                         pythiaPylistVerbosity = cms.untracked.int32(1),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         comEnergy = cms.double(13600.0),
                         PythiaParameters = cms.PSet(
                             pythia8CommonSettingsBlock,
                             pythia8CP5SettingsBlock,
                             processParameters = cms.vstring(
                                 'HiggsSM:gg2H = on',
                                 '25:m0 =125',
                                 '25:addChannel = 1 0.5 102 4900101 -4900101',
                                 '25:addChannel = 1 0.5 102 4900102 -4900102',
                                 '25:0:onMode=0',
                                 '25:1:onMode=0',
                                 '25:2:onMode=0',
                                 '25:3:onMode=0',
                                 '25:4:onMode=0',
                                 '25:5:onMode=0',
                                 '25:6:onMode=0',
                                 '25:7:onMode=0',
                                 '25:8:onMode=0',
                                 '25:9:onMode=0',
                                 '25:10:onMode=0',
                                 '25:11:onMode=0',
                                 '25:12:onMode=0',
                                 '25:13:onMode=0',
                                 'HiddenValley:Ngauge = 3',
                                 'HiddenValley:nFlav = 2',
                                 'HiddenValley:fragment = on',
                                 'HiddenValley:FSR = on',
                                 'HiddenValley:alphaOrder = 1',
                                 'HiddenValley:setLambda = on',  # known to only work in 8.309 or higher
                                 'HiddenValley:Lambda = 4.0',
                                 'HiddenValley:pTminFSR = 4.4',
                                 'HiddenValley:spinFv=0',
                                 'HiddenValley:probVector=0.75',
                                 'HiddenValley:separateFlav = on',  # unclear why not accepted
                                 'HiddenValley:probKeepEta1 = 1.0',  # unclear why not accepted
                                 '4900101:m0 = 4.162',
                                 '4900102:m0 = 4.198',  # this also fails
                                 '4900111:m0 = 1.2',
                                 '4900211:m0 = 1.2',
                                 '4900221:m0 = 5.0',
                                 '4900113:m0 = 4.0',
                                 '4900213:m0 = 4.0',
                                 '4900113:addChannel = 1 1.00 91 4900211 -4900211',
                                 '4900213:addChannel = 1 1.00 91 4900211 4900111',
                                 '999999:all = GeneralResonance void 1 0 0 0.5 0.001 0. 0. 0.',
                                 '999999:addChannel = 1 0.387 91 11 -11',
                                 '999999:addChannel = 1 0.383 91 13 -13',
                                 '999999:addChannel = 1 0.221 91 211 -211',
                                 '999999:addChannel = 1 0.00405 91 211 -211 111',
                                 '999999:tau0 = 0.000252',
                                 '4900221:addChannel = 1 0.5 91 4900211 4900211 4900111',
                                 '4900221:addChannel = 1 0.5 91 -4900211 -4900211 4900111',
                                 '4900221:tau0 = 0.0',
                                 '4900111:addChannel = 1 1.0 91 999999 999999',
                                 '4900211:onMode = 0',
                             ),
                             parameterSets = cms.vstring('pythia8CommonSettings',
                                                         'pythia8CP5Settings',
                                                         'processParameters',
                                                         )
                         )
                         )
