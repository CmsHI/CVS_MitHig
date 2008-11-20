import FWCore.ParameterSet.Config as cms

process = cms.Process("RECO")

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Simulation_cff")
process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("HeavyIonsAnalysis.Configuration.HeavyIonTracking_cff")

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
    skipEvents = cms.untracked.uint32(0),
    fileNames  = cms.untracked.vstring(
		'dcache:/pnfs/cmsaf.mit.edu/hibat/cms/users/edwenger/mix_and_digi/__MIXTAG__/__INPUT__.root'
	)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

##############################################################################
# Track Analyzer

process.testHighPtGlobalTracks = cms.EDAnalyzer("HighPtTrackAnalyzer",
						trackCollection = cms.vstring("globalPrimTracks"),
						resultFile      = cms.string("__OUTPUT__.root"),
						plotEvent = cms.bool(False),
						zipFiles  = cms.bool(False)
						)

##############################################################################
# Paths

process.p = cms.Path(process.hiTrackingWithOfflineBeamSpot * process.testHighPtGlobalTracks)