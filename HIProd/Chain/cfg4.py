import FWCore.ParameterSet.Config as cms

process = cms.Process("RECO")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Simulation_cff")
process.load("Configuration.StandardSequences.L1Emulator_cff")
process.load("Configuration.StandardSequences.L1TriggerDefaultMenu_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FakeConditions_cff")

process.load("HeavyIonsAnalysis.Configuration.Reconstruction_HI_cff")
from HLTrigger.Timer.timer_cfi import *
process.TimerService = cms.Service("TimerService",
                                   useCPUtime = cms.untracked.bool(True)
                                   )
process.recoTimer=myTimer

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:__INPUT__')
                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO'),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(-1)
        ),
        WARNING = cms.untracked.PSet(
            limit = cms.untracked.int32(-1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(-1)
        )
    ),
    destinations = cms.untracked.vstring('cout')
)

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
    oncePerEventMode = cms.untracked.bool(False),
    ignoreTotal = cms.untracked.int32(0)
)

process.options = cms.untracked.PSet(
    makeTriggerResults = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('__OUTPUT__')
)

process.p = cms.Path(process.RawToDigi*process.reconstruct_PbPb)
process.outpath = cms.EndPath(process.recoTimer+process.output)




