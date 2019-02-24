from WMCore.Configuration import Configuration
config = Configuration()
from CRABClient.UserUtilities import getUsernameFromSiteDB
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.section_('JobType')
config.JobType.outputFiles = ['check.root']
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/home/sanders/pPb_2016/CheckFlattening/checkflattening_cfg.py'
config.JobType.maxJobRuntimeMin = 1315
config.section_('Data')
config.Data.unitsPerJob = 20
config.Data.publication = False
config.Data.splitting = 'LumiBased'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_US_Vanderbilt'

config.General.requestName = 'pPb2016_pPb_MB1_Check'
config.Data.outLFNDirBase = '/store/user/ssanders/pPb2016_pPb_MB1_Check'
config.Data.lumiMask = 'Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt'
config.JobType.inputFiles = ['Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt','HeavyIonRPRcd_pPb2016_MB_offline.db','offset_pPb2016_MB_1_285832.root']
config.JobType.pyCfgParams = ['noprint','ntrkMin=0','lumifile=Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt','offset=offset_pPb2016_MB_1_285832.root']
config.Data.inputDataset = '/PAMinimumBias1/PARun2016C-PromptReco-v1/AOD'
config.Data.runRange = '285479-285832'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_projects'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print hte.headers

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################
submit(config)

config.General.requestName = 'pPb2016_Pbp_MB1_Check'
config.Data.outLFNDirBase = '/store/user/ssanders/pPb2016_Pbp_MB1_Check'
config.Data.lumiMask = 'Cert_285952-286009_HI8TeV_PromptReco_Pbp_Collisions16_JSON_NoL1T.txt'
config.JobType.inputFiles = ['Cert_285952-286009_HI8TeV_PromptReco_Pbp_Collisions16_JSON_NoL1T.txt','HeavyIonRPRcd_pPb2016_MB_offline.db','offset_pPb2016_MB_285833_286009.root']
config.JobType.pyCfgParams = ['noprint','ntrkMin=0','lumifile=Cert_285952-286009_HI8TeV_PromptReco_Pbp_Collisions16_JSON_NoL1T.txt','offset=offset_pPb2016_MB_285833_286009.root']
config.Data.inputDataset = '/PAMinimumBias1/PARun2016C-PromptReco-v1/AOD'
config.Data.runRange = '285952-286009'
submit(config)

