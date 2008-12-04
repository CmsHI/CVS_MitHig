import FWCore.ParameterSet.Config as cms

process = cms.Process("GEN")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Generator_cff")

process.load("RecoJets.Configuration.RecoGenJets_cff")
process.load("RecoJets.JetProducers.IconeJetParameters_cfi")
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#process.load("GeneratorInterface.PyquenInterface.pyquenSourceDefault_cfi")
#process.PyquenSource.doQuench = False

process.load("GeneratorInterface.HydjetInterface.hydjetSourceDefault_cfi")
process.HydjetSource.allowEmptyEvents = cms.bool(False)
process.HydjetSource.hydjetMode = 'kQJetsOnly'
#process.HydjetSource.hydjetMode = 'kJetsOnly'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
                                        ignoreTotal=cms.untracked.int32(0),
                                        oncePerEventMode = cms.untracked.bool(False)
                                        )

process.iterativeCone5HiGenJets = cms.EDProducer("IterativeConeHiGenJetProducer",
                                                 process.GenJetParameters,  
                                                 process.IconeJetParameters, 
                                                 jetType = cms.untracked.string('GenJet'), 
                                                 alias = cms.untracked.string('IC5HiGenJet'),
                                                 coneRadius = cms.double(0.5)
                                                 )


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   sourceSeed = cms.untracked.uint32(1)
                                                   )

process.subevent = cms.EDAnalyzer('HeavyIonJetAnalyzer',
                                                               jetSrc = cms.vstring('iterativeCone5HiGenJets')
                                                               )

process.allevent =  cms.EDAnalyzer('HeavyIonJetAnalyzer',
                                                                      jetSrc = cms.vstring('iterativeCone5GenJets'),
                                                                      doParticles = cms.untracked.bool(False)
                                                                      )

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('genjets.root')
                                   )

process.output = cms.OutputModule("PoolOutputModule",
    compressionLevel = cms.untracked.int32(2),
    commitInterval = cms.untracked.uint32(1),
    fileName = cms.untracked.string('PbPb_edm.root')
)

process.runjets = cms.Sequence(process.GeneInfo*process.genJetParticles*process.iterativeCone5GenJets*process.iterativeCone5HiGenJets)
process.p = cms.Path(process.runjets*process.subevent*process.allevent)
process.outpath = cms.EndPath(process.output)



