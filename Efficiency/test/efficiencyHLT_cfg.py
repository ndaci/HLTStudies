import FWCore.ParameterSet.Config as cms

process = cms.Process("EFF")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:hlt_100_1_knd.root'
    )
)

process.load("HLTStudies.Efficiency.efficiencyHLT_cfi")
process.p = cms.Path(process.eff)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('file:outeff.root')
)
