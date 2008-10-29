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

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('__OUTPUT__')
)

process.outpath = cms.EndPath(process.output)




