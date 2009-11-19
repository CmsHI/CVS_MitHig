import FWCore.ParameterSet.Config as cms

process = cms.Process("myRECO")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

# Timing service
process.Timing = cms.Service("Timing") 

process.GlobalTag.globaltag = 'STARTUP31X_V4::All'

process.pixelVertexFromClusters = cms.EDProducer('PixelVertexProducerClusters')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:RelValMinBias_314_STARTUP31X_V2-v1-Reco.root'
    )
)


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   moduleSeeds = cms.PSet(simMuonRPCDigis = cms.untracked.uint32(49835),
                                                                          simEcalUnsuppressedDigis = cms.untracked.uint32(49835),
                                                                          simSiStripDigis = cms.untracked.uint32(49835),
                                                                          mix = cms.untracked.uint32(49835),
                                                                          simHcalUnsuppressedDigis = cms.untracked.uint32(49835),
                                                                          simMuonCSCDigis = cms.untracked.uint32(49835),
                                                                          VtxSmeared = cms.untracked.uint32(49835),
                                                                          g4SimHits = cms.untracked.uint32(49835),
                                                                          simMuonDTDigis = cms.untracked.uint32(49835),
                                                                          simSiPixelDigis = cms.untracked.uint32(49835)
                                                                          ),
                                                   sourceSeed = cms.untracked.uint32(49835)
                                                   )


process.ana = cms.EDAnalyzer('PixelHitAnalyzer',
                             vertexSrc = cms.vstring('pixelVertices'),
                             trackSrc = cms.InputTag('pixelTracks'),
                             doTracking = cms.untracked.bool(False)
                             )

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('PixelTree-4.root')
                                   )


process.analyze = cms.Path(process.siPixelRecHits*process.ana)


