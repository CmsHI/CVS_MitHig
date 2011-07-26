import FWCore.ParameterSet.Config as cms

anaTrack = cms.EDAnalyzer('TrackAnalyzer',
                          trackPtMin = cms.untracked.double(0.4),
                          simTrackPtMin = cms.untracked.double(0.4),
                          vertexSrc = cms.vstring('hiSelectedVertex'),
                          trackSrc = cms.InputTag('hiGoodTightTracks'),
                          pfCandSrc = cms.InputTag('particleFlow'),
                          doPFMatching = cms.untracked.bool(False),
                          doSimTrack = cms.untracked.bool(False)
                          )
