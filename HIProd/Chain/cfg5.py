import FWCore.ParameterSet.Config as cms

process = cms.Process("TIME")
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:__INPUT__')
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('__OUTPUT__.hist')
)

process.ana = cms.EDFilter("RecoTimeAnalyzer",
                           doMixed = cms.untracked.bool(False)
                           )

process.p1 = cms.Path(process.ana)


