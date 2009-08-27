import FWCore.ParameterSet.Config as cms

process = cms.Process("RECO")

from Configuration.Generator.PyquenDefaultSettings_cff import *

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RecoHI.Configuration.Reconstruction_HI_cff")

process.MessageLogger.debugModules = cms.untracked.vstring("mix")
                             
process.source = cms.Source('EmptySource')

process.load("GeneratorInterface.PyquenInterface.pyquenDefault_cfi")
process.generator.doQuench = False
process.generator.embeddingMode = False
process.generator.doIsospin = False
process.generator.comEnergy = 10000
process.generator.PythiaParameters = cms.PSet(pythiaUESettings = cms.vstring('MSTJ(11)=3     ! Choice of the fragmentation function',
                                                                             'MSTJ(22)=2     ! Decay those unstable particles',
                                                                             'PARJ(71)=10 .  ! for which ctau  10 mm',
                                                                             'MSTP(2)=1      ! which order running alphaS',
                                                                             'MSTP(33)=0     ! no K factors in hard cross sections',
                                                                             'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)',
                                                                             'MSTP(52)=2     ! work with LHAPDF',
                                                                             'MSTP(81)=1     ! multiple parton interactions 1 is Pythia default',
                                                                             'MSTP(82)=4     ! Defines the multi-parton model',
                                                                             'MSTU(21)=1     ! Check on possible errors during program execution',
                                                                             'PARP(82)=1.8387   ! pt cutoff for multiparton interactions',
                                                                             'PARP(89)=1960. ! sqrts for which PARP82 is set',
                                                                             'PARP(83)=0.5   ! Multiple interactions: matter distrbn parameter',
                                                                             'PARP(84)=0.4   ! Multiple interactions: matter distribution parameter',
                                                                             'PARP(90)=0.16  ! Multiple interactions: rescaling power',
                                                                             'PARP(67)=2.5    ! amount of initial-state radiation',
                                                                             'PARP(85)=1.0  ! gluon prod. mechanism in MI',
                                                                             'PARP(86)=1.0  ! gluon prod. mechanism in MI',
                                                                             'PARP(62)=1.25   ! ',
                                                                             'PARP(64)=0.2    ! ',
                                                                             'MSTP(91)=1      !',
                                                                             'PARP(91)=2.1   ! kt distribution',
                                                                             'PARP(93)=15.0  ! '),
                                              processParameters = cms.vstring('MSEL=1               ! QCD hight pT processes',
                                                                              'CKIN(3)=230.          ! minimum pt hat for hard interactions',
                                                                              'CKIN(4)=300.         ! maximum pt hat for hard interactions'),
                                              parameterSets = cms.vstring('pythiaUESettings',
                                                                          'processParameters')
                                              )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
                                        ignoreTotal=cms.untracked.int32(0),
                                        oncePerEventMode = cms.untracked.bool(False)
                                        )

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.Generator_cff")
process.load("Configuration.StandardSequences.VtxSmearedGauss_cff")
process.load("Configuration.StandardSequences.Simulation_cff")
process.load("Configuration.StandardSequences.MixingNoPileUp_cff")
process.load("Configuration.StandardSequences.Digi_cff")
process.load("Configuration.StandardSequences.L1Emulator_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_31X_V3::All'

process.load("CmsHi.Utilities.HiGenParticles_cfi")
process.load("RecoHI.Configuration.RecoHI_EventContent_cff")

process.output = cms.OutputModule("PoolOutputModule",
                                  process.FEVTEventContent,
                                  compressionLevel = cms.untracked.int32(2),
                                  commitInterval = cms.untracked.uint32(1),
                                  fileName = cms.untracked.string('pyquen_dijet_pt230to300.root')
                                  )

# Paths
process.gen = cms.Path(process.generator*process.VtxSmeared*process.hiGenParticles)
process.sim = cms.Path(process.psim)
process.digi = cms.Path(process.mix*process.doAllDigi*process.trackingParticles*process.L1Emulator*process.DigiToRaw*process.RawToDigi)
process.reco = cms.Path(process.reconstruct_PbPb)
# End Path
process.end = cms.EndPath(process.output)

# Schedule
process.schedule = cms.Schedule(process.sim,process.gen,process.digi,process.reco,process.end)
