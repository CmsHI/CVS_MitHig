import FWCore.ParameterSet.Config as cms

#
# Tracker Local Reco
# Initialize magnetic field
#
from RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitConverter_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitMatcher_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.StripCPEfromTrackAngle_cfi import *
from RecoLocalTracker.SiStripZeroSuppression.SiStripZeroSuppression_SimData_cfi import * #sim
from RecoLocalTracker.SiStripClusterizer.SiStripClusterizer_SimData2_cfi import * #sim
from RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi import *
siPixelClusters.src = 'simSiPixelDigis' #sim
from RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi import *
pixeltrackerlocalrecoSim = cms.Sequence(siPixelClusters*siPixelRecHits)
striptrackerlocalrecoSim = cms.Sequence(siStripZeroSuppression*siStripClusters*siStripMatchedRecHits)
trackerlocalrecoSim = cms.Sequence(pixeltrackerlocalrecoSim*striptrackerlocalrecoSim)


