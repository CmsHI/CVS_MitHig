import FWCore.ParameterSet.Config as cms

process = cms.Process("RECO")

from Configuration.Generator.PyquenDefaultSettings_cff import *

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RecoHI.Configuration.Reconstruction_HI_cff")
#process.load("CmsHi.EgammaAnalysis.EcalGenTrigger_cff")

process.MessageLogger.debugModules = cms.untracked.vstring("mix")
                             
process.source = cms.Source('PoolSource',
                            fileNames = cms.untracked.vstring('dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/mc/Summer09/Hydjet_MinBias_4TeV/GEN-SIM-RAW/MC_31X_V2-GaussianVtx_311_ver1/0000/8C53B062-9673-DE11-94C6-001EC94BF0EF.root'),
                            inputCommands = cms.untracked.vstring('keep *',
                                                                  'drop *_*rawData*_*_*',
                                                                  'drop *_*Digis_*_*',
                                                                  'drop *_genParticles_*_*'
                                                                  )
                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

from GeneratorInterface.PyquenInterface.pyquenDefault_cfi import *
process.signal = generator.clone()
process.signal.doQuench = False
process.signal.embeddingMode = True
process.signal.PythiaParameters.parameterSets = cms.vstring('pythiaDefault','pythiaJets','kinematics')
process.signal.PythiaParameters.kinematics = cms.vstring('CKIN(3) = 120','CKIN(4) = 170')

process.RandomNumberGeneratorService.signal = cms.PSet(
    initialSeed = cms.untracked.uint32(1)
    )

process.RandomNumberGeneratorService.signalSIM = cms.PSet(process.RandomNumberGeneratorService.g4SimHits)

process.RandomNumberGeneratorService.signalSIM.initialSeed = 1

from CmsHi.Utilities.EventEmbedding_cff import *
process.mix=cms.EDProducer('HiEventEmbedder',
                           simEventEmbeddingMixParameters,
                           signalTag = cms.vstring("signal","signalSIM")
                           )

process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
                                        ignoreTotal=cms.untracked.int32(0),
                                        oncePerEventMode = cms.untracked.bool(False)
                                        )

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Generator_cff")
process.load("Configuration.StandardSequences.Digi_cff")
process.load("Configuration.StandardSequences.L1Emulator_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_31X_V3::All'

process.load("CmsHi.Utilities.HiGenParticles_cfi")
process.load("CmsHi.Utilities.HiAnalysisEventContent_cff")
process.load("CmsHi.Utilities.MatchVtx_cfi")

from Configuration.StandardSequences.Simulation_cff import *

process.signalSIM = g4SimHits
process.signalSIM.Generator.HepMCProductLabel = 'signal'
process.genParticles.src = cms.InputTag("signal")

process.output = cms.OutputModule("PoolOutputModule",
                                  process.HITrackAnalysisObjects,
                                  compressionLevel = cms.untracked.int32(2),
                                  commitInterval = cms.untracked.uint32(1),
                                  fileName = cms.untracked.string('pyquen_mixed_into_hydjet.root')
                                  )

process.load('Configuration.EventContent.EventContent_cff')
process.HITrackAnalysisObjects.outputCommands.extend(process.RAWEventContent.outputCommands)

process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.mergedtruth.HepMCDataLabels = ['signal']                      # by default: 'VtxSmeared', 'PythiaSource', 'source' (and 'generator' in 3_1_x)

# My Filter
#process.filter = cms.Sequence(process.partontrig100*process.ecaltrig100)

# Paths
process.sim = cms.Path(process.signal*process.matchVtx*process.signalSIM*process.mix)
process.gen = cms.Path(process.hiGenParticles+process.genParticles)
process.digi = cms.Path(process.doAllDigi*process.trackingParticles*process.L1Emulator*process.DigiToRaw*process.RawToDigi)
process.reco = cms.Path(process.reconstruct_PbPb_CaloOnly)
# End Path
process.end = cms.EndPath(process.output)

# Schedule
process.schedule = cms.Schedule(process.sim,process.gen,process.digi,process.reco,process.end)
