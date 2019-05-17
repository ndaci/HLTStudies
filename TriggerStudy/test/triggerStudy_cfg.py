import FWCore.ParameterSet.Config as cms

process = cms.Process("TRG")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'root://cms-xrd-global.cern.ch//store/data/Run2017B/SingleElectron/MINIAOD/17Nov2017-v1/70000/029A4179-3DE6-E711-AC8B-0025905A60B0.root'#,
        'root://cms-xrd-global.cern.ch//store/data/Run2018A/EGamma/MINIAOD/17Sep2018-v2/20000/306C7738-A329-7E45-8710-3724D419CFD8.root'
    )
)

from HLTStudies.TriggerStudy.triggerStudy_cfi import triggerStudy
process.triggerStudy = triggerStudy

process.p = cms.Path(process.triggerStudy)
