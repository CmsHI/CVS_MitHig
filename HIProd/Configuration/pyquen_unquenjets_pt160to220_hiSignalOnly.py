import FWCore.ParameterSet.Config as cms

process = cms.Process("RECO")

from Configuration.Generator.PyquenDefaultSettings_cff import *

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RecoHI.Configuration.Reconstruction_hiSignal_cff")

process.MessageLogger.debugModules = cms.untracked.vstring("mix")
                             
process.source = cms.Source('PoolSource',
                            fileNames = cms.untracked.vstring('dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/mc/Summer09/Hydjet_MinBias_4TeV/GEN-SIM-RAW/MC_31X_V2-GaussianVtx_311_ver1/0000/8C53B062-9673-DE11-94C6-001EC94BF0EF.root'),

                            inputCommands = cms.untracked.vstring('keep *',
                                                                  'drop *_*rawData*_*_*',
                                                                  'drop *_*Digis_*_*',
                                                                  'drop *_genParticles_*_*',
                                                                  'drop *_OTHER_LOCAL_RECO'
                                                                  )
                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
                                        ignoreTotal=cms.untracked.int32(0),
                                        oncePerEventMode = cms.untracked.bool(False)
                                        )

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Digi_cff")
process.load("Configuration.StandardSequences.L1Emulator_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_31X_V3::All'

process.load("CmsHi.Utilities.HiGenParticles_cfi")
process.load("CmsHi.Utilities.HiAnalysisEventContent_cff")

from CmsHi.Utilities.EventEmbedding_cff import *
process.mix = noMix.clone()

from CmsHi.Utilities.HiGenParticles_cfi import *
process.hiSignalGenParticles = hiGenParticles.clone()
process.hiSignalGenParticles.src = cms.InputTag("signal") #"hiSignal"

process.output = cms.OutputModule("PoolOutputModule",
                                  process.HITrackAnalysisObjects,
                                  compressionLevel = cms.untracked.int32(2),
                                  commitInterval = cms.untracked.uint32(1),
                                  fileName = cms.untracked.string('pyquen_unquenjets_pt160to220_hydjet_quen_b0_4TeV.root')
                                  )

process.load('Configuration.EventContent.EventContent_cff')
process.HITrackAnalysisObjects.outputCommands.extend(process.RAWEventContent.outputCommands)

process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.mergedtruth.HepMCDataLabels = ['signal','hiSignal']  # Label "signal" to be replaced with "hiSignal" in the next production 

# My Filter
#process.filter = cms.Sequence(process.partontrig100*process.ecaltrig100)

# Paths
process.restart = cms.Path(process.mix)
process.gen = cms.Path(process.hiSignalGenParticles)
process.digi = cms.Path(process.doAllDigi*process.trackingParticles*process.L1Emulator*process.DigiToRaw*process.RawToDigi)
process.reco = cms.Path(process.reconstruct_hiSignalOnly)

# End Path
process.end = cms.EndPath(process.output)

# Schedule
process.schedule = cms.Schedule(process.restart,process.gen,process.digi,process.reco,process.end)
