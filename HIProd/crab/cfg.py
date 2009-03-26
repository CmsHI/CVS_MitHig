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
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load("GeneratorInterface.PyquenInterface.pyquenSourceDefault_cfi")

process.PyquenSource.comEnergy = 4000
process.PyquenSource.doQuench = False

process.PyquenSource.PythiaParameters.parameterSets = cms.vstring('pythiaDefault','csa08Settings','pythiaJets','pythiaPromptPhotons','mySettings')
process.PyquenSource.PythiaParameters.mySettings = cms.vstring(
                                                               'CKIN(3)=0.',    #PtHatMin Lowered
                                                               'MSTP(81)=1'     #Pythia Underlying Event ON
                                                               )

process.GlobalTag.globaltag = 'IDEAL_V11::All'

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(100)
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
                                fileName = cms.untracked.string('pyquen_off_sim_mb_4TeV_vtxFixed_d20090326.root')
                                )

process.pgen0 = cms.Sequence(process.VertexSmearing+process.GeneInfo)
process.p0 = cms.Path(process.pgen0)
process.p1 = cms.Path(process.psim)
#process.p2 = cms.Path(process.mix*process.doAllDigi*process.L1Emulator*process.DigiToRaw*process.RawToDigi)
process.outpath = cms.EndPath(process.FEVT)

process.VtxSmeared.MinX = 0.0001
process.VtxSmeared.MaxX = 0.0001
process.VtxSmeared.MinY = 0.0002
process.VtxSmeared.MaxY = 0.0002
process.VtxSmeared.MinZ = 2.
process.VtxSmeared.MaxZ = 2.



