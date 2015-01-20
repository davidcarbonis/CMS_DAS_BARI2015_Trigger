import FWCore.ParameterSet.Config as cms

process = cms.Process("TestTriggerReport")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['TriggerReport.txt']
process.MessageLogger.categories.append('HLTrigReport')
process.MessageLogger.categories.append('L1GtTrigReport')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace this file with the one you actually want to use
    fileNames = cms.untracked.vstring(
        'file:///data/shared/Short_Exercise_HLT/data/CMSDAS_HLT_DoubleMuSkim_10k.root'
    )
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_hlt_GRun', '')

process.load("HLTrigger.HLTanalyzers.hltTrigReport_cfi")
process.hltTrigReport.HLTriggerResults = cms.InputTag( 'TriggerResults','','HLT') 

process.load("L1Trigger.GlobalTriggerAnalyzer.l1GtTrigReport_cfi")
process.l1GtTrigReport.L1GtRecordInputTag = cms.InputTag( 'gtDigis' )

process.HLTAnalyzerEndpath = cms.EndPath( process.l1GtTrigReport + process.hltTrigReport )
