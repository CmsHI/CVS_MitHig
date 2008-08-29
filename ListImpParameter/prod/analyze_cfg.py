import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:__INPUT__.root')
                            )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('__OUTPUT__.root')
                                   )

process.ana = cms.EDFilter("ListImpParameter",
                           )

process.p1 = cms.Path(process.ana)

