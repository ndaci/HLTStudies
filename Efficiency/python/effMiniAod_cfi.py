import FWCore.ParameterSet.Config as cms

eff = cms.EDAnalyzer(
    'MiniAodEff',
    verbose        = cms.int32(0),
    namePaths      = cms.vstring(""),
    metfilter      = cms.InputTag("TriggerResults","","PAT"),
    bits           = cms.InputTag("TriggerResults","","HLT"),
    prescales      = cms.InputTag("patTrigger"),
    objects        = cms.InputTag("selectedPatTrigger"),
    vertices       = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons          = cms.InputTag("slimmedMuons"),
    electrons      = cms.InputTag("slimmedElectrons"),
    taus           = cms.InputTag("slimmedTaus"),
    photons        = cms.InputTag("slimmedPhotons"),
    jets           = cms.InputTag("slimmedJets"),
    genjets        = cms.InputTag("slimmedGenJets"),
    met            = cms.InputTag("slimmedMETs"),
    usePtMHT       = cms.bool(True) ,
    minPtJetHt     = cms.double(20.0) , 
    maxEtaJetHt    = cms.double(5.2) , 
    minPtJetMht    = cms.double(20.0) , 
    maxEtaJetMht   = cms.double(5.2)
)
