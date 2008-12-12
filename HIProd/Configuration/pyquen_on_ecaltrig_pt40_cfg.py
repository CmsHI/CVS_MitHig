# Auto generated configuration file
# using: 
# Revision: 1.77.2.1 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/ecaltrig_cff.py -s GEN:ProductionFilterSequence,SIM --eventcontent RAWSIM --datatier GEN-SIM -n 10 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff') 
process.load('Configuration/StandardSequences/Geometry_cff') # Ecal Pre-Shower Included
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedFlat_cff') # For Fixed Vertex Position (For mixing)
process.load('Configuration/StandardSequences/Sim_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.6 $'),
    annotation = cms.untracked.string('GeneratorInterface/PyquenInterface/python/pyquenSourceDefault_cfi.py nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

# Old physics list due to Memory constraints
process.g4SimHits.Physics.type = "SimG4Core/Physics/QGSP_EMV"

# Input source
process.source = cms.Source("PyquenSource",

    doQuench = cms.bool(True),

    qgpInitialTemperature = cms.double(1.0),
    cFlag = cms.int32(0),
    doCollisionalEnLoss = cms.bool(True),
    bFixed = cms.double(0.0),
    angularSpectrumSelector = cms.int32(0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    PythiaParameters = cms.PSet(
        myParameters = cms.vstring(),
        pythiaDefault = cms.vstring('MSEL=1', 
            'MSTU(21)=1', 
            'PARU(14)=1.', 
            'MSTP(81)=0', 
            'PMAS(5,1)=4.8', 
            'PMAS(6,1)=175.0', 
            'CKIN(3)=7.', 
            'MSTJ(22)=2', 
            'PARJ(71)=10.'),
        csa08Settings = cms.vstring('MSEL=0', 
            'PARP(67)=1.', 
            'PARP(82)=1.9', 
            'PARP(85)=0.33', 
            'PARP(86)=0.66', 
            'PARP(89)=1000.', 
            'PARP(91)=1.0', 
            'MSTJ(11)=3', 
            'MSTJ(22)=2'),
        pythiaMuonCandidates = cms.vstring('CKIN(3)=20', 
            'MSTJ(22)=2', 
            'PARJ(71)=40.'),
        pythiaJets = cms.vstring('MSUB(11)=1', 
            'MSUB(12)=1', 
            'MSUB(13)=1', 
            'MSUB(28)=1', 
            'MSUB(53)=1', 
            'MSUB(68)=1'),
        pythiaPromptPhotons = cms.vstring('MSUB(14)=1', 
            'MSUB(18)=1', 
            'MSUB(29)=1', 
            'MSUB(114)=1', 
            'MSUB(115)=1'),
        parameterSets = cms.vstring('pythiaDefault', 
            'csa08Settings', 
            'pythiaJets', 
            'pythiaPromptPhotons', 
            'kinematics', 
            'multipleScatter'),
        kinematics = cms.vstring('CKIN(3)=40', 
            'CKIN(4)=9999', 
            'CKIN(7)=-4.', 
            'CKIN(8)=4.'),
        multipleScatter = cms.vstring('MSTP(81)=0')
    ),
    qgpProperTimeFormation = cms.double(0.1),
    comEnergy = cms.double(5500),
    numQuarkFlavor = cms.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    maxEventsToPrint = cms.untracked.int32(0),
    aBeamTarget = cms.double(207.0),
    doRadiativeEnLoss = cms.bool(True)
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('ecaltrig_cff_py_GEN_SIM.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'IDEAL_V9::All'
process.ecaltrig = cms.EDFilter("MCSingleParticleFilter",
    Status = cms.untracked.vint32(2, 2, 2, 2, 2, 
        2, 1, 1, 2, 2, 
        1, 1, 1),
    MaxEta = cms.untracked.vdouble(3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 
        3, 3, 3),
    MinEta = cms.untracked.vdouble(-3, -3, -3, -3, -3, 
        -3, -3, -3, -3, -3, 
        -3, -3, -3),
    MinPt = cms.untracked.vdouble(17.5, 17.5, # eta
                                  16., 16., # eta'
                                  18., 18., # omega
                                  14.5, 14.5, # pi
                                  17.5,       # pi0
                                  15,       # K0
                                  17.5, 17.5,   # e
                                  38.5         # gamma           
                                  ),

    ParticleID = cms.untracked.vint32(221, -221, # eta
                                      331, -331, # eta'
                                      223, -223, # omega
                                      211, -211, # pi
                                      111,       # pi0
                                      311,       # K0
                                      11, -11,   # e
                                      22         # gamma
                                      )
)

process.partontrig = cms.EDFilter("MCSingleParticleFilter",
                                  Status = cms.untracked.vint32(2,2,2,2,2,2,2,1),
                                  MaxEta = cms.untracked.vdouble(3,3,3,3,3,3,3,3),
                                  MinEta = cms.untracked.vdouble(-3,-3,-3,-3,-3,-3,-3,-3),
                                  MinPt = cms.untracked.vdouble(38.5,38.5,38.5,38.5,38.5,38.5,38.5,38.5),
                                  ParticleID = cms.untracked.vint32(1,2,3,4,5,6,21,22)
                                  )

process.ProductionFilterSequence = cms.Sequence(process.partontrig*process.ecaltrig)

process.VtxSmeared.MinX = 0.0001
process.VtxSmeared.MaxX = 0.0001
process.VtxSmeared.MinY = 0.0002
process.VtxSmeared.MaxY = 0.0002
process.VtxSmeared.MinZ = 2.
process.VtxSmeared.MaxZ = 2.

process.pgenhi = cms.Sequence(cms.SequencePlaceholder("randomEngineStateProducer")+process.VertexSmearing+process.GeneInfo)

# Path and EndPath definitions
process.generation_step = cms.Path(process.ProductionFilterSequence*process.pgenhi)
process.simulation_step = cms.Path(process.ProductionFilterSequence*process.psim)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.simulation_step,process.out_step)
