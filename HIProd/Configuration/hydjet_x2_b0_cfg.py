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
process.load('Configuration/StandardSequences/VtxSmearedGauss_cff') # For Heavy-Ion conditions
process.load('Configuration/StandardSequences/Sim_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('Hydjet, Most Central (b = 0 fm)'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/UserCode/MitHig/HIProd/Configuration/hydjet_x2_b0_cfg.py,v $')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
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
    bMax = cms.double(0),
    cFlag = cms.int32(0),
    firstEvent = cms.untracked.uint32(1),
    hydjetMode = cms.string('kHydroQJets'),
    hadronFreezoutTemperature = cms.double(0.14),
    nMultiplicity = cms.int32(26000),
    qgpNumQuarkFlavor = cms.int32(0),
    doRadiativeEnLoss = cms.bool(True),
    bFixed = cms.double(0),
    maxLongitudinalRapidity = cms.double(3.75),
    bMin = cms.double(0),
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
	pythiaBoson = cms.vstring(
            'MSUB(1)=1',
            'MSUB(2)=1'),
        pythiaCharmoniumNRQCD = cms.vstring(                              
            'MSUB(421) = 1',
            'MSUB(422) = 1',
            'MSUB(423) = 1',
            'MSUB(424) = 1',
            'MSUB(425) = 1',
            'MSUB(426) = 1',
            'MSUB(427) = 1',
            'MSUB(428) = 1',
            'MSUB(429) = 1',
            'MSUB(430) = 1',
            'MSUB(431) = 1',
            'MSUB(432) = 1',
            'MSUB(433) = 1',
            'MSUB(434) = 1',
            'MSUB(435) = 1',
            'MSUB(436) = 1',
            'MSUB(437) = 1',
            'MSUB(438) = 1',
            'MSUB(439) = 1',
#            'MSTP(142) = 2',
            'PARJ(13)=0.60',
            'PARJ(14)=0.162',
            'PARJ(15)=0.018',
            'PARJ(16)=0.054',
            'MSTP(145)=0',
            'MSTP(146)=0',
            'MSTP(147)=0',
            'MSTP(148)=1',
            'MSTP(149)=1',
            'PARP(141)=1.16',
            'PARP(142)=0.0119',
            'PARP(143)=0.01',
            'PARP(144)=0.01',
            'PARP(145)=0.05',
            'PARP(146)=9.28',
            'PARP(147)=0.15',
            'PARP(148)=0.02',
            'PARP(149)=0.02',
            'PARP(150)=0.085',
            'BRAT(861)=0.202',
            'BRAT(862)=0.798',
            'BRAT(1501)=0.013',
            'BRAT(1502)=0.987',
            'BRAT(1555)=0.356',
            'BRAT(1556)=0.644'),
        pythiaBottomoniumNRQCD = cms.vstring(
            'MSUB(461) = 1',
            'MSUB(462) = 1',
            'MSUB(463) = 1',
            'MSUB(464) = 1',
            'MSUB(465) = 1',
            'MSUB(466) = 1',
            'MSUB(467) = 1',
            'MSUB(468) = 1',
            'MSUB(469) = 1',
            'MSUB(470) = 1',
            'MSUB(471) = 1',
            'MSUB(472) = 1',
            'MSUB(473) = 1',
            'MSUB(474) = 1',
            'MSUB(475) = 1',
            'MSUB(476) = 1',
            'MSUB(477) = 1',
            'MSUB(478) = 1',
            'MSUB(479) = 1',
 #           'MSTP(142)=2',    
            'PARJ(13)=0.750',  
            'PARJ(14)=0.162',  
            'PARJ(15)=0.018',  
            'PARJ(16)=0.054',
            'MSTP(145)=0',    
            'MSTP(146)=0',     
            'MSTP(147)=0',   
            'MSTP(148)=1',    
            'MSTP(149)=1',     
            'PARP(141)=1.16', 
            'PARP(142)=0.0119',
            'PARP(143)=0.01',  
            'PARP(144)=0.01',  
            'PARP(145)=0.05',  
            'PARP(146)=9.28', 
            'PARP(147)=0.15',  
            'PARP(148)=0.02',  
            'PARP(149)=0.02',  
            'PARP(150)=0.085'),                 
#        decaytiming = cms.vstring('PARJ(71)=40'),                # only for gen-level
        parameterSets = cms.vstring('pythiaDefault', 
            'csa08Settings', 
            'pythiaJets', 
            'pythiaPromptPhotons',
            'pythiaBoson',
            'pythiaCharmoniumNRQCD',
            'pythiaBottomoniumNRQCD')
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
process.GlobalTag.globaltag = 'IDEAL_V12::All'

# pgen re-defined to exclude all genJets
process.pgenhi = cms.Sequence(cms.SequencePlaceholder("randomEngineStateProducer")+process.VertexSmearing+process.GeneInfo)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgenhi)
process.simulation_step = cms.Path(process.psim)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.simulation_step,process.out_step)

