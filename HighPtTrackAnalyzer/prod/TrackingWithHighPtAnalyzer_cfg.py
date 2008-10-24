import FWCore.ParameterSet.Config as cms

process = cms.Process("RECO")

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Simulation_cff")
process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Digi_cff")

process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.load("RecoLocalTracker.Configuration.RecoLocalTracker_SimData_cff")
#process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")

process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("RecoVertex.BeamSpotProducer.BeamSpotFakeConditionsEarlyCollision_cff")
#process.load("HeavyIonsAnalysis.Configuration.HeavyIonTracking_cff")
process.load("HeavyIonsAnalysis.Configuration.HighPtTracking_PbPb_cff")

###############################################################################
# Message logger
process.MessageLogger = cms.Service("MessageLogger",
    cerr = cms.untracked.PSet( threshold = cms.untracked.string('INFO') ),
    destinations = cms.untracked.vstring('cerr'),
	#categories = cms.untracked.vstring('MinBiasTracking'),
	debugModules = cms.untracked.vstring()
)

###############################################################################
# Source
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames  = cms.untracked.vstring(
        #'dcache:/pnfs/cmsaf.mit.edu/hibat/cms/users/yetkin/sim/pythia_dijet_pt100to9999_d20081021/pythia_dijet_pt100to9999_d20081021_r000001.root'
		#'dcache:/pnfs/cmsaf.mit.edu/hibat/cms/users/yetkin/sim/pythia_mb_5TeV_d20080823/pythia_mb_5TeV_d20080823_r000002.root'
		#'dcache:/pnfs/cmsaf.mit.edu/hibat/cms/users/yetkin/sim/hydjet_sim_x2_mb_d20080819/hydjet_sim_x2_mb_d20080819_r000211.root'
		#'dcache:/pnfs/cmsaf.mit.edu/hibat/cms/users/yetkin/sim/hydjet_sim_x2_c1_d20080819/hydjet_sim_x2_c1_d20080819_r000002.root'
		'dcache:/pnfs/cmsaf.mit.edu/hibat/cms/users/yetkin/sim/hydjet_x2_b0_d20081020/hydjet_x2_b0_d20081020_r000116.root'
		)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
		
process.load("Configuration.StandardSequences.MixingNoPileUp_cff")			
					
#from HeavyIonsAnalysis.Configuration.EventEmbedding_cff import *
#process.mix=mixSim
#process.mix.input.fileNames = cms.untracked.vstring(
			#'dcache:/pnfs/cmsaf.mit.edu/hibat/cms/users/yetkin/sim/hydjet_sim_x2_mb_d20080819/hydjet_sim_x2_mb_d20080819_r000211.root'
			#'dcache:/pnfs/cmsaf.mit.edu/hibat/cms/users/yetkin/sim/pythia_dijet_pt100to9999_d20080731/pythia_dijet_pt100to9999_d20080731_r000037.root'
			#'dcache:/pnfs/cmsaf.mit.edu/hibat/cms/users/yetkin/sim/hydjet_x2_b0_d20081020/hydjet_x2_b0_d20081020_r000116.root'
#)			
				
###############################################################################
# Random Number Seeds
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                  moduleSeeds = cms.PSet(simMuonRPCDigis = cms.untracked.uint32(6),
										 simEcalUnsuppressedDigis = cms.untracked.uint32(8),
										 simSiStripDigis = cms.untracked.uint32(7),
										 mix = cms.untracked.uint32(4),
										 simHcalUnsuppressedDigis = cms.untracked.uint32(9),
										 simMuonCSCDigis = cms.untracked.uint32(6),
										 VtxSmeared = cms.untracked.uint32(2),
										 g4SimHits = cms.untracked.uint32(3),
										 simMuonDTDigis = cms.untracked.uint32(6),
										 simSiPixelDigis = cms.untracked.uint32(7)
										 ),
				  sourceSeed = cms.untracked.uint32(1)
                  )
				  
##############################################################################
# Track Analyzer

process.testHighPtGlobalTracks = cms.EDAnalyzer("HighPtTrackAnalyzer",
						trackCollection = cms.vstring("globalPrimTracks"),
						resultFile      = cms.string("RecoStudyOutput.root"),
						plotEvent = cms.bool(False),
						zipFiles  = cms.bool(False)
						)

##############################################################################
# Paths

process.p = cms.Path(process.mix
					 * process.trackingParticles
					 * process.offlineBeamSpot
					 * process.trDigi	  
					 * process.trackerlocalreco
					 * process.heavyIonTracking
					 * process.testHighPtGlobalTracks
					 )
					 
##############################################################################
# Output
#process.output = cms.OutputModule("PoolOutputModule",
#						fileName = cms.untracked.string("output_tracking_test.root")
#						)
#process.outpath = cms.EndPath(process.output)	
							   
