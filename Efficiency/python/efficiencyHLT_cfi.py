import FWCore.ParameterSet.Config as cms

eff = cms.EDAnalyzer(
    'Efficiency',
    hltProcessName   = cms.string("HLT"),
    genParticleLabel = cms.InputTag("genParticles"),
    jetCollection    = cms.InputTag("ak5PFJets")
)
