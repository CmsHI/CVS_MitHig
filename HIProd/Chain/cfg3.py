import FWCore.ParameterSet.Config as cms

process = cms.Process("DIGI")
process.load("Configuration.StandardSequences.Services_cff") 
process.load("Configuration.StandardSequences.VtxSmearedFlat_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Generator_cff")
process.load("Configuration.StandardSequences.Simulation_cff")
process.load("Configuration.StandardSequences.MixingNoPileUp_cff")
process.load("Configuration.StandardSequences.L1Emulator_cff")
process.load("Configuration.StandardSequences.L1TriggerDefaultMenu_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FakeConditions_cff")
from HLTrigger.Timer.timer_cfi import *
process.TimerService = cms.Service("TimerService",
                                   useCPUtime = cms.untracked.bool(True)
                                   )
process.digiTimer=myTimer

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

process.FEVT = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('__OUTPUT__')
)

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   moduleSeeds = cms.PSet(
    simMuonRPCDigis = cms.untracked.uint32(__RANDOM__),
    simEcalUnsuppressedDigis = cms.untracked.uint32(__RANDOM__),
    simSiStripDigis = cms.untracked.uint32(__RANDOM__),
    mix = cms.untracked.uint32(__RANDOM__),
    simHcalUnsuppressedDigis = cms.untracked.uint32(__RANDOM__),
    simMuonCSCDigis = cms.untracked.uint32(__RANDOM__),
    VtxSmeared = cms.untracked.uint32(__RANDOM__),
    g4SimHits = cms.untracked.uint32(__RANDOM__),
    simMuonDTDigis = cms.untracked.uint32(__RANDOM__),
    simHcalDigis = cms.untracked.uint32(__RANDOM__),
    simSiPixelDigis = cms.untracked.uint32(__RANDOM__)
    ),
                                                   sourceSeed = cms.untracked.uint32(__RANDOM__)
                                                   )


process.p2 = cms.Path(process.doAllDigi*process.L1Emulator*process.DigiToRaw)
#process.p2 = cms.Path(process.pdigi*process.L1Emulator*process.DigiToRaw)
process.outpath = cms.EndPath(process.FEVT)
#process.outpath = cms.EndPath(process.digiTimer+process.FEVT)



