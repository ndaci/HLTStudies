import FWCore.ParameterSet.Config as cms

triggerStudy = cms.EDAnalyzer('TriggerStudy',
  triggerResults = cms.InputTag('TriggerResults', '', 'HLT'),
  filterResults = cms.InputTag('TriggerResults', '', 'RECO'),
  prescales = cms.InputTag('patTrigger', '', 'RECO'),
  triggerObjects = cms.InputTag('slimmedPatTrigger', '', 'RECO'),
  triggerL1EG = cms.InputTag('caloStage2Digis', 'EGamma', 'RECO'),
  triggerL1Tau = cms.InputTag('caloStage2Digis', 'Tau', 'RECO'),
  triggerL1Jet = cms.InputTag('caloStage2Digis', 'Jet', 'RECO'),
  triggerL1Sums = cms.InputTag('caloStage2Digis', 'EtSum', 'RECO'),
  triggerL1Mu = cms.InputTag('gmtStage2Digis', 'Muon', 'RECO'),
  triggerL1algos = cms.InputTag('gtStage2Digis', '', 'RECO'),
  muons = cms.InputTag('slimmedMuons'),
  electrons = cms.InputTag('slimmedElectrons'),
  photons = cms.InputTag('slimmedPhotons'),
  taus = cms.InputTag('slimmedTaus'),
  jets = cms.InputTag('slimmedJets'),
  t1met = cms.InputTag('slimmedMETs')
)
