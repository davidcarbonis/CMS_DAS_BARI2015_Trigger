import FWCore.ParameterSet.Config as cms

process = cms.Process('TRIGANA')

# import of standard configurations
process.load("TrackingTools.PatternTools.TSCBLBuilderNoMaterial_cfi")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source

process.source = cms.Source("PoolSource",
    # replace this file with the one you actually want to use
fileNames = cms.untracked.vstring(
  'file:///data/shared/Short_Exercise_HLT/data/CMSDAS_HLT_DoubleMuSkim_10k.root')
)


# Additional output definition

#process.ana = cms.EDAnalyzer('HLTDASexercise',
process.ana = cms.EDAnalyzer('TriggerMuMuAnalysis',

                             OutFileName=cms.untracked.string("TriggerProductionHistos.root"),
                             TriggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
                             TriggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                             OfflineMuonsTag = cms.untracked.InputTag("muons","","RECO"),
                             BeamSpotTag = cms.untracked.InputTag("offlineBeamSpot","","RECO"),

                             PathName = cms.untracked.string("HLT_Mu17_Mu8_v"),
                             L3MuonFilter = cms.untracked.InputTag("hltDiMuonGlb17Glb8DzFiltered0p2","","HLT"),
                             DiMuonVertexTag = cms.untracked.InputTag("hltDisplacedmumuFilterDoubleMu4Jpsi","","HLT"),

                             maxDRforRecomatching = cms.untracked.double(0.1),

                             debug = cms.untracked.bool(False),
)

# Other statements
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_hlt_GRun', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:hltonline_8E33v2', '')

# Path and EndPath definitions
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("HLTexercise2012D.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

process.ana_step = cms.Path(process.ana)
