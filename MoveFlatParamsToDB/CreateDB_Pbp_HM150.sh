#!/bin/bash
rm Pbp_HM150.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM150_1/* > Pbp_HM150.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM150_2/* >> Pbp_HM150.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM150_3/* >> Pbp_HM150.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM150_4/* > Pbp_HM150.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM150_5/* >> Pbp_HM150.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM150_6/* >> Pbp_HM150.lis
mkdir RescorTables_Pbp_HM150
root -l -b -q 'EPCalib/EPCalib.C+(285833,286009,"Pbp_HM150")'
rm offsets/offset_pPb2016_Pbp_HM150.root
mv foff_Pbp_HM150.root offsets/offset_pPb2016_Pbp_HM150.root
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_pPb2016_Pbp_HM150.db outputTag=HeavyIonRPRcd begin=285833 end=286009 infile="/rfs/sanders/EP_Pbp_HM150.root" rescor="RescorTables_Pbp_HM150"
rm /rfs/sanders/tmp_Pbp_HM150
rm save/EP_pPb2016_Pbp_HM150.root
mv /rfs/sanders/EP_Pbp_HM150.root save/EP_pPb2016_Pbp_HM150.root
rm -rf RescorSave/RescorTables_pPb2016_Pbp_HM150
mv RescorTables_Pbp_HM150 RescorSave/RescorTables_pPb2016_Pbp_HM150
