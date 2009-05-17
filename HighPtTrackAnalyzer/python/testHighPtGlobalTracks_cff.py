import FWCore.ParameterSet.Config as cms

# track associator settings
from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *	 # sim to reco associator
TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')      # matching fraction: shared/nrechit! default='sim', shared/nsimhit
TrackAssociatorByHits.UseGrouped = cms.bool(False)                   # grouping hits on overlap layers? default=True

# high pt track analyzer settings
testHighPtGlobalTracks = cms.EDAnalyzer("HighPtTrackAnalyzer",
										trackCollection = cms.vstring("globalPrimTracks"),
										resultFile      = cms.string("TrkStudyOutput.root"),
										useAbsoluteNumberOfHits = cms.untracked.bool(False),  # in SimToReco associator
										keepLowPtSimTracks = cms.untracked.bool(False),  # keep pt<2 GeV/c simtracks
										infoHiEventTopology = cms.untracked.bool(True)   # store heavy ion event info (e.g. npart)
										)
