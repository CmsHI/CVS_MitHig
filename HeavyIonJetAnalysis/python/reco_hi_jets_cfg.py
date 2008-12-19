import FWCore.ParameterSet.Config as cms

process = cms.Process("JET")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Generator_cff")
process.load("Configuration.StandardSequences.MixingNoPileUp_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Simulation_cff")
process.load("Configuration.StandardSequences.L1Emulator_cff")
process.load("Configuration.StandardSequences.L1TriggerDefaultMenu_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FakeConditions_cff")

process.load("HeavyIonsAnalysis.Configuration.Reconstruction_HI_cff")

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('dcache:__INPUT__'
                                                              )
                            )


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
                                                 alias = cms.untracked.string('IC5HiGenJet'),
                                                 coneRadius = cms.double(0.5)
                                                 )


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   moduleSeeds = cms.PSet(simMuonRPCDigis = cms.untracked.uint32(__RANDOM__),
                                                                          simEcalUnsuppressedDigis = cms.untracked.uint32(__RANDOM__),
                                                                          simSiStripDigis = cms.untracked.uint32(__RANDOM__),
                                                                          mix = cms.untracked.uint32(__RANDOM__),
                                                                          simHcalUnsuppressedDigis = cms.untracked.uint32(__RANDOM__),
                                                                          simMuonCSCDigis = cms.untracked.uint32(__RANDOM__),
                                                                          VtxSmeared = cms.untracked.uint32(__RANDOM__),
                                                                          g4SimHits = cms.untracked.uint32(__RANDOM__),
                                                                          simMuonDTDigis = cms.untracked.uint32(__RANDOM__),
                                                                          simSiPixelDigis = cms.untracked.uint32(__RANDOM__)
                                                                          ),
                                                   sourceSeed = cms.untracked.uint32(__RANDOM__)
                                                   )


process.subevent = cms.EDAnalyzer('HeavyIonJetAnalyzer',
                                                               jetSrc = cms.vstring('iterativeCone5HiGenJets')
                                                               )

process.allevent =  cms.EDAnalyzer('HeavyIonJetAnalyzer',
                                                                      jetSrc = cms.vstring('iterativeCone5GenJets'),
                                                                      doParticles = cms.untracked.bool(False)
                                                                      )


process.recoevent = cms.EDAnalyzer('HeavyIonJetAnalyzer',
                                   jetSrc = cms.vstring('iterativeConePu5CaloJets'),
                                   doParticles = cms.untracked.bool(False)
                                   )


process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('__OUTPUT__.hist')
                                   )

process.output = cms.OutputModule("PoolOutputModule",
    compressionLevel = cms.untracked.int32(2),
    commitInterval = cms.untracked.uint32(1),
    fileName = cms.untracked.string('__OUTPUT__')
)

process.myjets = cms.Sequence(process.GeneInfo*process.runjets*process.iterativeCone5HiGenJets)
process.pre = cms.Path(process.mix*process.doAllDigi*process.L1Emulator*process.DigiToRaw*process.RawToDigi*process.ecalloc*process.hcalLocalRecoSequence*process.hiCentrality)
process.p = cms.Path(process.myjets*process.subevent*process.allevent*process.recoevent)
process.outpath = cms.EndPath(process.output)



