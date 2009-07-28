# Auto generated configuration file
# using: 
# Revision: 1.123 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: HYDJET_MinBias_4TeV_cfi.py -s GEN:ProductionFilterSequence,SIM,DIGI,L1,DIGI2RAW --conditions FrontierConditions_GlobalTag,MC_31X_V3::All --datatier GEN-SIM-RAW --eventcontent FEVTDEBUGHLT -n 10 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('DIGI2RAW')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedEarly10TeVCollision_cff')
process.load('Configuration/StandardSequences/Sim_cff')
process.load('Configuration/StandardSequences/Digi_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load('Configuration/StandardSequences/DigiToRaw_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.3 $'),
    annotation = cms.untracked.string('Hydjet-B0 at 4TeV'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/CMSSW/GeneratorInterface/HydjetInterface/python/hydjetDefault_cfi.py,v $')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("EmptySource")

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('HYDJET_MinBias_4TeV_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RAW'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'MC_31X_V3::All'
process.generator = cms.EDFilter("HydjetGeneratorFilter",
    nMultiplicity = cms.int32(26000),
    qgpInitialTemperature = cms.double(1.0),
    doCollisionalEnLoss = cms.bool(True),
    maxLongitudinalRapidity = cms.double(3.75),
    shadowingSwitch = cms.int32(0),
    maxTransverseRapidity = cms.double(1.0),
    PythiaParameters = cms.PSet(
        pythiaUpsilonToMuons = cms.vstring('BRAT(1034) = 0 ', 
            'BRAT(1035) = 1 ', 
            'BRAT(1036) = 0 ', 
            'BRAT(1037) = 0 ', 
            'BRAT(1038) = 0 ', 
            'BRAT(1039) = 0 ', 
            'BRAT(1040) = 0 ', 
            'BRAT(1041) = 0 ', 
            'BRAT(1042) = 0 ', 
            'MDME(1034,1) = 0 ', 
            'MDME(1035,1) = 1 ', 
            'MDME(1036,1) = 0 ', 
            'MDME(1037,1) = 0 ', 
            'MDME(1038,1) = 0 ', 
            'MDME(1039,1) = 0 ', 
            'MDME(1040,1) = 0 ', 
            'MDME(1041,1) = 0 ', 
            'MDME(1042,1) = 0 '),
        myParameters = cms.vstring(),
        pythiaDefault = cms.vstring('MSEL=0', 
            'MSTU(21)=1', 
            'PARU(14)=1.', 
            'MSTP(81)=0', 
            'PMAS(5,1)=4.8', 
            'PMAS(6,1)=175.0', 
            'CKIN(3)=7.', 
            'MSTJ(22)=2', 
            'PARJ(71)=10.', 
            'PARP(67)=1.', 
            'PARP(82)=1.9', 
            'PARP(85)=0.33', 
            'PARP(86)=0.66', 
            'PARP(89)=1000.', 
            'PARP(91)=1.0', 
            'MSTJ(11)=3', 
            'MSTJ(22)=2'),
        pythiaBottomoniumNRQCD = cms.vstring('MSUB(461) = 1', 
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
            'MSUB(479) = 1'),
        pythiaCharmoniumNRQCD = cms.vstring('MSUB(421) = 1', 
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
            'MSUB(439) = 1'),
        pythiaMuonCandidates = cms.vstring('CKIN(3)=20', 
            'MSTJ(22)=2', 
            'PARJ(71)=40.'),
        pythiaQuarkoniaSettings = cms.vstring('PARP(141)=1.16', 
            'PARP(142)=0.0119', 
            'PARP(143)=0.01', 
            'PARP(144)=0.01', 
            'PARP(145)=0.05', 
            'PARP(146)=9.28', 
            'PARP(147)=0.15', 
            'PARP(148)=0.02', 
            'PARP(149)=0.02', 
            'PARP(150)=0.085', 
            'PARJ(13)=0.60', 
            'PARJ(14)=0.162', 
            'PARJ(15)=0.018', 
            'PARJ(16)=0.054', 
            'MSTP(145)=0', 
            'MSTP(146)=0', 
            'MSTP(147)=0', 
            'MSTP(148)=1', 
            'MSTP(149)=1', 
            'BRAT(861)=0.202', 
            'BRAT(862)=0.798', 
            'BRAT(1501)=0.013', 
            'BRAT(1502)=0.987', 
            'BRAT(1555)=0.356', 
            'BRAT(1556)=0.644'),
        pythiaWeakBosons = cms.vstring('MSUB(1)=1', 
            'MSUB(2)=1'),
        pythiaJets = cms.vstring('MSUB(11)=1', 
            'MSUB(12)=1', 
            'MSUB(13)=1', 
            'MSUB(28)=1', 
            'MSUB(53)=1', 
            'MSUB(68)=1'),
        pythiaZtoMuons = cms.vstring('MDME(174,1)=0', 
            'MDME(175,1)=0', 
            'MDME(176,1)=0', 
            'MDME(177,1)=0', 
            'MDME(178,1)=0', 
            'MDME(179,1)=0', 
            'MDME(182,1)=0', 
            'MDME(183,1)=0', 
            'MDME(184,1)=1', 
            'MDME(185,1)=0', 
            'MDME(186,1)=0', 
            'MDME(187,1)=0'),
        pythiaPromptPhotons = cms.vstring('MSUB(14)=1', 
            'MSUB(18)=1', 
            'MSUB(29)=1', 
            'MSUB(114)=1', 
            'MSUB(115)=1'),
        pythiaJpsiToMuons = cms.vstring('BRAT(858) = 0 ', 
            'BRAT(859) = 1 ', 
            'BRAT(860) = 0 ', 
            'MDME(858,1) = 0 ', 
            'MDME(859,1) = 1 ', 
            'MDME(860,1) = 0 '),
        parameterSets = cms.vstring('pythiaDefault', 
            'pythiaJets', 
            'pythiaPromptPhotons')
    ),
    qgpNumQuarkFlavor = cms.int32(0),
    qgpProperTimeFormation = cms.double(0.1),
    doRadiativeEnLoss = cms.bool(True),
    hydjetMode = cms.string('kHydroQJets'),
    sigmaInelNN = cms.double(58),
    fracSoftMultiplicity = cms.double(1.0),
    hadronFreezoutTemperature = cms.double(0.14),
    aBeamTarget = cms.double(208.0),
    allowEmptyEvents = cms.bool(False),
    rotateEventPlane = cms.bool(True),
    cFlag = cms.int32(1),
    firstRun = cms.untracked.uint32(1),
    bFixed = cms.double(0),
    bMin = cms.double(0.0),
    firstEvent = cms.untracked.uint32(1),
    comEnergy = cms.double(4000.0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    maxEventsToPrint = cms.untracked.int32(0),
    bMax = cms.double(30.0)
)
process.ProductionFilterSequence = cms.Sequence(process.generator)

# pgen re-defined to exclude all genJets
process.pgenhi = cms.Sequence(cms.SequencePlaceholder("randomEngineStateProducer")+process.VertexSmearing+process.GeneInfo)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgenhi)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.endjob_step,process.out_step)

# special treatment in case of production filter sequence  
for path in process.paths: 
    getattr(process,path)._seq = process.ProductionFilterSequence*getattr(process,path)._seq
