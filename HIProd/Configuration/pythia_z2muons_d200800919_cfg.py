#
# SCRAM_DIR: the directory where the intended version of CMSSW is
# JOB_HEPMC, JOB_DIGI: the input and output directory of the data
# JOB_BASE, JOB_TAG: place where condor job is based
# JOB_VERSION: a string indicating the version of the job
# JOB_GRID: (default) gatekeeper to use when submitting to the grid
#

#CFGSTART
#
# SCRAM_ARCH=slc4_ia32_gcc345
# SCRAM_DIR=/net/hisrv0001/home/yetkin/cms217
#
# JOB_INPUT=
# JOB_DIGI=
#
# JOB_TAG=pythia_z2muons_d200800919
# JOB_BASE=/net/pstore01/d00/scratch/yetkin/prod/pp/${JOB_TAG}
#
# JOB_VERSION=217
#
# JOB_RUN_HIROOT=
# JOB_OUTPUT=/pnfs/cmsaf.mit.edu/hibat/cms/users/yetkin/sim/${JOB_TAG}
#
#
# JOB_USERANDOM=1
#
#CFGEND

#VARSTART
#
# BASE=__BASE__
#
#VAREND

import FWCore.ParameterSet.Config as cms

process = cms.Process("SIM")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.VtxSmearedFlat_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Generator_cff")
process.load("Configuration.StandardSequences.Simulation_cff")
process.load("Configuration.StandardSequences.MixingNoPileUp_cff")
process.load("Configuration.StandardSequences.L1Emulator_cff")
process.load("Configuration.StandardSequences.L1TriggerDefaultMenu_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FakeConditions_cff")


process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(100)
            )

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
                                                                          simHcalDigis = cms.untracked.uint32(__RANDOM__),
                                                                          simSiPixelDigis = cms.untracked.uint32(__RANDOM__)
                                                                          ),
                                                   sourceSeed = cms.untracked.uint32(__RANDOM__)
                                                   )

process.source = cms.Source("PythiaSource",
                            pythiaHepMCVerbosity = cms.untracked.bool(False),
                            pythiaPylistVerbosity = cms.untracked.int32(0),
                            PythiaParameters = cms.PSet(myParameters = cms.vstring(),
                                                        parameterSets = cms.vstring('default',
                                                                                    'kinematics',
                                                                                    'process'),
                                                        default = cms.vstring('MSTU(21)=1'),
                                                        kinematics = cms.vstring('CKIN(3) = 20',
                                                                                'CKIN(7)=-1',
                                                                                'CKIN(8)=1'),
                                                        process = cms.vstring("MSEL=0         ! User defined processes", 
                                                                                    "MSUB(15)=1",             # !qqbar->g+Z,
                                                                                    "MSUB(30)=1",             # !qg->q+Z,
                                                                                    "MSTP(43)=2",             # !Only Z0,
                                                                                    "MDME(174,1)=0",          # !Z decay into d dbar,
                                                                                    "MDME(175,1)=0",          # !Z decay into u ubar,
                                                                                    "MDME(176,1)=0",          # !Z decay into s sbar,
                                                                                    "MDME(177,1)=0",          # !Z decay into c cbar,
                                                                                    "MDME(178,1)=0",          # !Z decay into b bbar,
                                                                                    "MDME(179,1)=0",          # !Z decay into t tbar,
                                                                                    "MDME(182,1)=0",          # !Z decay into e- e+,
                                                                                    "MDME(183,1)=0",          # !Z decay into nu_e nu_ebar,
                                                                                    "MDME(184,1)=1",          # !Z decay into mu- mu+,
                                                                                    "MDME(185,1)=0",          # !Z decay into nu_mu nu_mubar,
                                                                                    "MDME(186,1)=0",          # !Z decay into tau- tau+,
                                                                                    "MDME(187,1)=0"           # !Z decay into nu_tau nu_taubar
                                                                              )
                                                        )
                            )

process.MessageLogger = cms.Service("MessageLogger",
                                        debugModules = cms.untracked.vstring('PythiaSource'),
                                        cout = cms.untracked.PSet(
            threshold = cms.untracked.string('INFO'),
                    default = cms.untracked.PSet(
                limit = cms.untracked.int32(-1)
                        ),
                    WARNING = cms.untracked.PSet(
                limit = cms.untracked.int32(-1)
                        ),
                    INFO = cms.untracked.PSet(
                limit = cms.untracked.int32(-1)
                        )
                ),
                                        destinations = cms.untracked.vstring('cout')
                                    )

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
                                            oncePerEventMode = cms.untracked.bool(True),
                                            ignoreTotal = cms.untracked.int32(0)
                                        )

process.options = cms.untracked.PSet(
        makeTriggerResults = cms.untracked.bool(True),
            wantSummary = cms.untracked.bool(True)
        )

process.FEVT = cms.OutputModule("PoolOutputModule",
                                fileName = cms.untracked.string('__BASE__.root')
                                )

process.pgen0 = cms.Sequence(process.VertexSmearing+process.GeneInfo)
process.p0 = cms.Path(process.pgen0)
process.p1 = cms.Path(process.psim)
process.p2 = cms.Path(process.mix*process.doAllDigi*process.L1Emulator*process.DigiToRaw*process.RawToDigi)
process.outpath = cms.EndPath(process.FEVT)

process.VtxSmeared.MinX = 0.0001
process.VtxSmeared.MaxX = 0.0001
process.VtxSmeared.MinY = 0.0002
process.VtxSmeared.MaxY = 0.0002
process.VtxSmeared.MinZ = 2.
process.VtxSmeared.MaxZ = 2.



