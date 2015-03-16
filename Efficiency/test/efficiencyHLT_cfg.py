import os,sys
import FWCore.ParameterSet.Config as cms

process = cms.Process("EFF")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(9500) )

tag="test"
for i in range(0,len(sys.argv)):
    if str(sys.argv[i])=="_input" and len(sys.argv)>i+1:
        tag = str(sys.argv[i+1])

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:hlt_100_1_knd.root'
        'file:hltSimMonojet_Phys14_14e33_100_1_SOi.root'
    )
)

from list_MonojetM1AV_r731_i731HLT1_V13b_14e33_cff    import *
from list_MonojetM1AV_r731_i731HLT1_V13c_14e33_cff    import *
from list_MonojetM1AV_r731_i731HLT1_V13d_14e33_v1_cff import *
from list_MonojetM1AV_r731_i731HLT1_V13b_7e33_V2_cff  import *
from list_MonojetM1AV_r731_i731HLT1_V13b_5e33_V3_cff  import *
from list_MonojetM1AV_r731_i731HLT1_V13_7e33_V2_cff   import *
from list_files_DYMuMu_EDM_r731_i731HLT1_V13b_5e33_v5_cff import *
from list_files_DYMuMu_EDM_r731_i731HLT1_V13b_7e33_v5_cff import *
from list_files_DYMuMu_EDM_r731_i731HLT1_V13c_14e33_cff import *
from list_MonojetM10AV_r740p8_i740p7HLT1_V8_14e33_v5_cff import *
from list_MonojetM10AV_r740p8_i740p7HLT1_V8_7e33_v5_cff import *
from list_MonojetM1AV_r740p8_i740p7HLT1_V8_5e33_v5_cff import *

if   tag=="M10_74X_V8_14e33_v5":
    process.source = source_M10_74X_V8_14e33_v5
elif tag=="M10_74X_V8_7e33_v5":
    process.source = source_M10_74X_V8_7e33_v5
elif tag=="M1_74X_V8_5e33_v5":
    process.source = source_M1_74X_V8_5e33_v5
elif tag=="V13b_7e33_V2":
    process.source = source_V13b_7e33_V2
elif tag=="V13_7e33_V2":
    process.source = source_V13_7e33_V2
elif tag=="V13b_5e33_V3":
    process.source = source_V13b_5e33_V3
elif tag=="V13b_14e33":
    process.source = source_V13b_14e33
elif tag=="V13c_14e33":
    process.source = source_V13c_14e33
elif tag=="V13d_14e33":
    process.source = source_V13d_14e33
elif tag=="DYMuMu_V13b_5e33_v5":
    process.source = source_DYMuMu_V13b_5e33_v5
elif tag=="DYMuMu_V13b_7e33_v5":
    process.source = source_DYMuMu_V13b_7e33_v5
elif tag=="DYMuMu_V13c_14e33":
    process.source = source_DYMuMu_V13c_14e33


process.load("HLTStudies.Efficiency.efficiencyHLT_cfi")
process.eff.namePaths = cms.vstring("")
process.eff.hltProcessName   = cms.string("TEST")

# because in this GEN-SIM-RECO there are no ak4PFJets, this is scandalous
if tag=="DYMuMu_V13c_14e33": 
    process.eff.pfjetCollection  = cms.InputTag("ak5PFJets")

process.p = cms.Path(process.eff)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('file:outeff_'+tag+'.root')
)
