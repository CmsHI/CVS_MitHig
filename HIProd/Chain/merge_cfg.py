import FWCore.ParameterSet.Config as cms

process = cms.Process("MERGE")

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(__INPUT__)
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
    oncePerEventMode = cms.untracked.bool(True),
    ignoreTotal = cms.untracked.int32(0)
)

myOutputCommands = cms.untracked.vstring('keep *',
                                         'drop *Digi*_*_*_*',
                                         'keep *_*sim*Digi*_*_*',
                                         'drop *CrossingFrame*_*_*_*',
                                         'drop *Raw*_*_*_*',
                                         'drop *_*_*_RECO'
                                         )


process.output = cms.OutputModule("PoolOutputModule",
                                  outputCommands = cms.untracked.vstring(),
                                  fileName = cms.untracked.string('__OUTPUT__')
                                  )

process.output.outputCommands = myOutputCommands

process.outpath = cms.EndPath(process.output)




