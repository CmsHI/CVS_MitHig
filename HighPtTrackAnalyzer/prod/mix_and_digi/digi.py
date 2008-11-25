import FWCore.ParameterSet.Config as cms

process = cms.Process("DIGI")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Simulation_cff")
process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")

###############################################################################
# Timing service
process.Timing = cms.Service("Timing")

###############################################################################
# Memory Check
process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
	oncePerEventMode = cms.untracked.bool(False),
	ignoreTotal = cms.untracked.int32(0)
)

###############################################################################
# Message logger
process.MessageLogger = cms.Service("MessageLogger",
    cerr = cms.untracked.PSet( threshold = cms.untracked.string('INFO') ),
    destinations = cms.untracked.vstring('cerr'),
	#categories = cms.untracked.vstring('MinBiasTracking'),
	#debugModules = cms.untracked.vstring('pixel3ProtoTracks','pixelVertices','pixel3PrimTracks','primSeeds')
)

###############################################################################
# Source
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(__SKIPEVENT__),
    fileNames  = cms.untracked.vstring(
		'dcache:/pnfs/cmsaf.mit.edu/hibat/cms/users/yetkin/sim/__DATATAG__/__DATATAG___r__RUN__.root'
	)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
		
process.load("Configuration.StandardSequences.MixingNoPileUp_cff")			
					
#from HeavyIonsAnalysis.Configuration.EventEmbedding_cff import *
#process.mix=mixSim
#process.mix.input.fileNames = cms.untracked.vstring('dcache:/pnfs/cmsaf.mit.edu/hibat/cms/users/yetkin/sim/__MIXTAG__/__MIXTAG___r__MIXRUN__.root')			
				
###############################################################################
# Random Number Seeds
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

##############################################################################
# RawToDigi

process.rawDataCollector.currentProcessOnly = True

##############################################################################
# Paths

process.p = cms.Path(process.mix
					 * process.trackingParticles
					 * process.doAllDigi	  
					 * process.DigiToRaw
					 * process.RawToDigi
					 )
					 
##############################################################################
# Output
process.output = cms.OutputModule("PoolOutputModule",
	fileName = cms.untracked.string('__OUTPUT__.root'),
	outputCommands = cms.untracked.vstring('drop *',
											'keep edmHepMCProduct_*_*_*',
											'keep *_mergedtruth_*_*',
											'keep PSimHitCrossingFrame_mix_g4SimHitsTracker*_DIGI',
											'keep *_simSiPixelDigis_*_DIGI',
											'keep *_simSiStripDigis_*_DIGI',
											'keep *_siStripDigis_*_DIGI',
											'keep *_siPixelDigis_*_DIGI'
    )
)
process.save = cms.EndPath(process.output)	
							   
