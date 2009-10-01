# Auto generated configuration file
# using: 
# Revision: 1.2 
# Source: /cvs_server/repositories/CMSSW/UserCode/edwenger/Misc/ConfigBuilder.py,v 
# with command line options: hiReco -s RAW2DIGI,RECO --scenario HeavyIons --conditions FrontierConditions_GlobalTag,MC_31X_V8::All --datatier GEN-SIM-RECO --eventcontent=RECODEBUG --filein=/store/relval/CMSSW_3_3_0_pre3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_31X_V8-v1/0015/22EFEE98-4EA1-DE11-BEFC-001D09F2B30B.root --fileout=22EFEE98-4EA1-DE11-BEFC-001D09F2B30B.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/RawToDigi_cff')
process.load('Configuration/StandardSequences/ReconstructionHeavyIons_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContentHeavyIons_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.2 $'),
    annotation = cms.untracked.string('hiReco nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/relval/CMSSW_3_3_0_pre3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_31X_V8-v1/0015/22EFEE98-4EA1-DE11-BEFC-001D09F2B30B.root')
)

process.Timing = cms.Service("Timing")

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
                                        oncePerEventMode = cms.untracked.bool(False),
                                        ignoreTotal = cms.untracked.int32(0)
                                        )

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RECODEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('edmfile.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'MC_31X_V8::All'

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstructionHeavyIons)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.endjob_step,process.out_step)
