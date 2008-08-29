import FWCore.ParameterSet.Config as cms

process = cms.Process("GEN")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Generator_cff")
process.load("GeneratorInterface.Pythia6Interface.PythiaSourceMinBias_cfi")
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   sourceSeed = cms.untracked.uint32(__RANDOM__)
                                                   )

process.PythiaSource.comEnergy = cms.untracked.double(10000)

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(100)
        )

process.FEVT = cms.OutputModule("PoolOutputModule",
                                fileName = cms.untracked.string('__BASE__.root')
                                )
process.outpath = cms.EndPath(process.FEVT)

