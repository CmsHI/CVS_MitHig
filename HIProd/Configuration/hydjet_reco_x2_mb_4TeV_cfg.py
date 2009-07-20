import FWCore.ParameterSet.Config as cms

process = cms.Process("RECO")

from GeneratorInterface.PyquenInterface.pyquenPythiaDefault_cff import *


process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RecoHI.Configuration.Reconstruction_HI_cff")

process.source = cms.Source('PoolSource',
                            fileNames = cms.untracked.vstring(__INPUT__),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
                                        ignoreTotal=cms.untracked.int32(0),
                                        oncePerEventMode = cms.untracked.bool(False)
                                        )

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MixingNoPileUp_cff")
process.load("Configuration.StandardSequences.Digi_cff")
process.load("Configuration.StandardSequences.L1Emulator_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_31X_V3::All'

process.load("CmsHi.Utilities.HiAnalysisEventContent_cff")
process.load("CmsHi.Utilities.HiGenParticles_cfi")

process.hiGenParticles.src = ['generator']

process.output = cms.OutputModule("PoolOutputModule",
                                  process.HITrackAnalysisObjects,
                                  compressionLevel = cms.untracked.int32(2),
                                  commitInterval = cms.untracked.uint32(1),
                                  fileName = cms.untracked.string('__OUTPUT__'),
                                  )

process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.mergedtruth.HepMCDataLabels = ['generator']                      # by default: 'VtxSmeared', 'PythiaSource', 'source' (and 'generator' in 3_1_x)

# Paths
process.gen = cms.Path(process.hiGenParticles)
process.digi = cms.Path(process.RawToDigi)
process.reco = cms.Path(process.reconstruct_PbPb_CaloOnly)
# End Path
process.end = cms.EndPath(process.output)

# Schedule
process.schedule = cms.Schedule(process.gen,process.digi,process.reco,process.end)
