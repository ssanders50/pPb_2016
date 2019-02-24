#!/bin/bash
rm Pbp_MB.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB1/* > Pbp_MB.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB2/* > Pbp_MB.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB3/* > Pbp_MB.lis
mkdir RescorTables_Pbp_MB
root -l -b -q 'EPCalib/EPCalib.C+(285833,286009,"Pbp_MB")'
rm offsets/offset_pPb2016_Pbp_MB.root
mv foff_Pbp_MB.root offsets/offset_pPb2016_Pbp_MB.root
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_pPb2016_Pbp_MB.db outputTag=HeavyIonRPRcd begin=285833 end=286009 infile="/rfs/sanders/EP_Pbp_MB.root" rescor="RescorTables_Pbp_MB"
rm /rfs/sanders/tmp_Pbp_MB
rm save/EP_pPb2016_Pbp_MB.root
mv /rfs/sanders/EP_Pbp_MB.root save/EP_pPb2016_Pbp_MB.root
rm -rf RescorSave/RescorTables_pPb2016_Pbp_MB
mv RescorTables_Pbp_MB RescorSave/RescorTables_pPb2016_Pbp_MB
