import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.FakeConditions_cff import *

siPixelCabling.connect = 'sqlite_file:/net/pstore01/d00/scratch/yetkin/cms/CondCore/SQLiteData/data/siPixelCabling200.db'
cscPackingCabling.connect = 'sqlite_file:/net/pstore01/d00/scratch/yetkin/cms/CondCore/SQLiteData/data/CSCChamberMapValues.db'
cscUnpackingCabling.connect = 'sqlite_file:/net/pstore01/d00/scratch/yetkin/cms/CondCore/SQLiteData/data/CSCCrateMapValues.db'
DTCabling.connect = 'sqlite_file:/net/pstore01/d00/scratch/yetkin/cms/CondCore/SQLiteData/data/DTFullMap_fix17X.db'
RPCCabling.connect = 'sqlite_file:/net/pstore01/d00/scratch/yetkin/cms/CondCore/SQLiteData/data/RPCEMap_181007.db'
trackProbabilityFakeCond.connect = 'sqlite_file:/net/pstore01/d00/scratch/yetkin/cms/CondCore/SQLiteData/data/btagTrackProbability200.db'
BTauMVAJetTagComputerRecord.connect = 'sqlite_file:/net/pstore01/d00/scratch/yetkin/cms/CondCore/SQLiteData/data/MVAJetTagsFakeConditions.db'


