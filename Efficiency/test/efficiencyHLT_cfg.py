import FWCore.ParameterSet.Config as cms

process = cms.Process("EFF")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:hlt_100_1_knd.root'
        'file:hltSimMonojet_Phys14_14e33_100_1_SOi.root'
    )
)

#from list_MonojetM1AV_r731_i731HLT1_V13_7e33_V2_cff import *
from list_MonojetM1AV_r731_i731HLT1_V13_14e33_cff import *
process.source = source

process.load("HLTStudies.Efficiency.efficiencyHLT_cfi")
process.eff.namePaths = cms.vstring("")
process.eff.hltProcessName   = cms.string("TEST")
process.p = cms.Path(process.eff)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('file:outeff.root')
)
