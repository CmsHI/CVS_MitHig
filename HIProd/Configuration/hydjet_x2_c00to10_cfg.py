# Auto generated configuration file
# using: 
# Revision: 1.77.2.1 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: hydjetSourceDefault_cfi.py -s GEN,SIM --eventcontent RAWSIM --datatier GEN-SIM --conditions FrontierConditions_GlobalTag,IDEAL_V9::All -n 10 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/Geometry_cff') # Ecal Pre-Shower In
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedFlat_cff') # For Mixing
process.load('Configuration/StandardSequences/Sim_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('Hydjet, 10% Most Central'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/UserCode/MitHig/HIProd/Configuration/hydjet_x2_c00to10_cfg.py,v $')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("HydjetSource",
    shadowingSwitch = cms.int32(0),
    maxTransverseRapidity = cms.double(1.0),
    comEnergy = cms.double(5500.0),
    sigmaInelNN = cms.double(58),
    allowEmptyEvents = cms.bool(False),
    doCollisionalEnLoss = cms.bool(True),
    qgpInitialTemperature = cms.double(1.0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    bMin = cms.double(0),
    bMax = cms.double(5.02707),
    cFlag = cms.int32(1),
    firstEvent = cms.untracked.uint32(1),
    hydjetMode = cms.string('kHydroQJets'),
    hadronFreezoutTemperature = cms.double(0.14),
    nMultiplicity = cms.int32(26000),
    qgpNumQuarkFlavor = cms.int32(0),
    doRadiativeEnLoss = cms.bool(True),
    bFixed = cms.double(0),
    maxLongitudinalRapidity = cms.double(3.75),
    fracSoftMultiplicity = cms.double(1.0),
    maxEventsToPrint = cms.untracked.int32(0),
    aBeamTarget = cms.double(208.0),
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
            'pythiaPromptPhotons')
    ),
    firstRun = cms.untracked.uint32(1),
    qgpProperTimeFormation = cms.double(0.1)
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('hydjetSourceDefault_cfi_py_GEN_SIM.root'),
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

# Fixed Vertex
process.VtxSmeared.MinX = 0.0001
process.VtxSmeared.MaxX = 0.0001
process.VtxSmeared.MinY = 0.0002
process.VtxSmeared.MaxY = 0.0002
process.VtxSmeared.MinZ = 2.
process.VtxSmeared.MaxZ = 2.

# pgen re-defined to exclude all genJets
process.pgenhi = cms.Sequence(cms.SequencePlaceholder("randomEngineStateProducer")+process.VertexSmearing+process.GeneInfo)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgenhi)
process.simulation_step = cms.Path(process.psim)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.simulation_step,process.out_step)
