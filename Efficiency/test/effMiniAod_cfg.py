import os,sys
import FWCore.ParameterSet.Config as cms

process = cms.Process("EFFMINIAOD")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

tag=""
for i in range(0,len(sys.argv)):
    if str(sys.argv[i])=="_input" and len(sys.argv)>i+1:
        tag = str(sys.argv[i+1])

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:hlt_100_1_knd.root'
        #'file:hltSimMonojet_Phys14_14e33_100_1_SOi.root'
        #'file:/user/ndaci/Data/DarkMonojet/TestHCAL/hltSim10000_Method3.root'
        #'file:/user/ndaci/Data/MiniAOD/miniAOD-prod_PAT.root'
        'file:/user/ndaci/Data/MiniAOD/store/relval/CMSSW_7_4_4/RelValADDMonoJet_d3MD3_13/MINIAODSIM/MCRUN2_74_V9_38Tbis-v1/00000/3A4505E8-2F09-E511-8459-003048FFD732.root'
    )
)

process.load("HLTStudies.Efficiency.effMiniAod_cfi")
process.eff.namePaths = cms.vstring("")
process.eff.hltProcessName   = cms.string("HLT")

process.p = cms.Path(process.eff)

# Output
if tag!="":
    tag = "_"+tag

process.TFileService = cms.Service('TFileService',
    fileName = cms.string('outeff'+tag+'.root')
)
