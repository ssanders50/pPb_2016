#!/bin/bash
rm Pbp_HM120.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM120_1/* > Pbp_HM120.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM120_2/* >> Pbp_HM120.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM120_3/* >> Pbp_HM120.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM120_4/* >> Pbp_HM120.lis
mkdir RescorTables_Pbp_HM120
root -l -b -q 'EPCalib/EPCalib.C+(285833,286009,"Pbp_HM120")'
rm offsets/offset_pPb2016_Pbp_HM120.root
mv foff_Pbp_HM120.root offsets/offset_pPb2016_Pbp_HM120.root
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_pPb2016_Pbp_HM120.db outputTag=HeavyIonRPRcd begin=285833 end=286009 infile="/rfs/sanders/EP_Pbp_HM120.root" rescor="RescorTables_Pbp_HM120"
rm /rfs/sanders/tmp_Pbp_HM120
rm save/EP_pPb2016_Pbp_HM120.root
mv /rfs/sanders/EP_Pbp_HM120.root save/EP_pPb2016_Pbp_HM120.root
rm -rf RescorSave/RescorTables_pPb2016_Pbp_HM120
mv RescorTables_Pbp_HM120 RescorSave/RescorTables_pPb2016_Pbp_HM120
