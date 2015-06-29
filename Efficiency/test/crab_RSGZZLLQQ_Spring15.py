from WMCore.Configuration import Configuration
config = Configuration()

name = 'MiniEff_Spring15_V3'
proc = 'RSGZZLLQQ_50ns'
dataset = '/RSGravToZZToLLQQ_kMpl01_M-1800_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/MINIAODSIM'

# GENERAL
config.section_("General")
config.General.requestName = name+"_"+proc
config.General.workArea    = '/user/ndaci/CRABBY/StudyHLT/Monojet/'+name+'/'+proc+'/'
config.General.transferLogs    = True
#config.General.transferOutput = True

# JOB TYPE
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'effMiniAod_cfg.py'
#config.JobType.pyCfgParams = ['reEmulation=True','reEmulMuons=True','reEmulCalos=True','patchNtuple=True','force2012Config=True','customDTTF=True','dttfLutsFile=sqlite:src/L1TriggerDPG/L1Menu/data/dttf_config.db','useUct2015=True','globalTag=POSTLS162_V2::All','runOnMC=True','runOnPostLS1=True','whichPU=40']
#config.JobType.inputFiles = '../../data/dttf_config.db'
config.JobType.allowUndistributedCMSSW = True

# INPUT DATA
config.section_("Data")
config.Data.inputDataset = dataset
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits  = 6
config.Data.publication = False
config.Data.publishDBS  = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.publishDataName = name+"_"+proc
config.Data.ignoreLocality = True # allows to process inputs on CE != site hosting inputs
#config.Data.lumiMask = 
#config.Data.runRange = 

#A custom string to insert in the output file name inside the CRAB-created directory path to allow organizing groups of tasks.
#config.Data.prefix =  

# USER
config.section_("User")
#config.User.email = 'nadir.daci@cern.ch'
#config.User.voRole = 
config.User.voGroup = 'becms'

# GRID
config.section_("Site")
config.Site.storageSite = 'T2_BE_IIHE'
#config.Site.whitelist = 
config.Site.blacklist = ['T1_US_FNAL']
