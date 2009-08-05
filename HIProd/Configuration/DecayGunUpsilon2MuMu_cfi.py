import FWCore.ParameterSet.Config as cms

from Configuration.Generator.PyquenDefaultSettings_cff import *

generator = cms.EDProducer("Pythia6PtGun",
                           maxEventsToPrint = cms.untracked.int32(0),
                           pythiaHepMCVerbosity = cms.untracked.bool(False),
                           pythiaPylistVerbosity = cms.untracked.int32(0),
                           
                           PGunParameters = cms.PSet(ParticleID = cms.vint32(553),
                                                     MinEta = cms.double(-2.5),
                                                     MaxEta = cms.double(2.5),
                                                     MinPhi = cms.double(-3.14159265359),
                                                     MaxPhi = cms.double(3.14159265359),
                                                     MinPt = cms.double(0),
                                                     MaxPt = cms.double(20),
                                                     AddAntiParticle = cms.bool(False)
                                                     ),

                           PythiaParameters = cms.PSet(pyquenPythiaDefaultBlock,
                                                       parameterSets = cms.vstring('pythiaDefault','pythiaUpsilontoMuons')
                                                       )
                           )


