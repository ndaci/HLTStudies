import FWCore.ParameterSet.Config as cms

eff = cms.EDAnalyzer(
    'Efficiency',
    namePaths        = cms.vstring(""),
    hltProcessName   = cms.string("HLT"),
    genParticleLabel = cms.string("genParticles"),
    pfjetCollection  = cms.InputTag("ak4PFJets"),
    pfmetCollection  = cms.InputTag("pfMet"),
    muCollection     = cms.InputTag("muons"),
    #vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    usePtMHT     = cms.bool(True) ,
    minPtJetHt   = cms.double(20.0) , 
    maxEtaJetHt  = cms.double(5.2) , 
    minPtJetMht  = cms.double(20.0) , 
    maxEtaJetMht = cms.double(5.2)
)
