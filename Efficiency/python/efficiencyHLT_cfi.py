import FWCore.ParameterSet.Config as cms

eff = cms.EDAnalyzer(
    'Efficiency',
    hltProcessName   = cms.string("HLT"),
    genParticleLabel = cms.string("genParticles"),
    pfjetCollection    = cms.InputTag("ak5PFJets"),
    pfmetCollection    = cms.InputTag("pfMet"),
    muCollection     = cms.InputTag("muons")
)
