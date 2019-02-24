#!/bin/bash
rm Pbp_HM250.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM250_1/* > Pbp_HM250.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM250_2/* >> Pbp_HM250.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM250_3/* >> Pbp_HM250.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM250_4/* >> Pbp_HM250.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM250_5/* >> Pbp_HM250.lis
mkdir RescorTables_Pbp_HM250
root -l -b -q 'EPCalib/EPCalib.C+(285833,286009,"Pbp_HM250")'
rm offsets/offset_pPb2016_Pbp_HM250.root
mv foff_Pbp_HM250.root offsets/offset_pPb2016_Pbp_HM250.root
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_pPb2016_Pbp_HM250.db outputTag=HeavyIonRPRcd begin=285833 end=286009 infile="/rfs/sanders/EP_Pbp_HM250.root" rescor="RescorTables_Pbp_HM250"
rm /rfs/sanders/tmp_Pbp_HM250
rm save/EP_pPb2016_Pbp_HM250.root
mv /rfs/sanders/EP_Pbp_HM250.root save/EP_pPb2016_Pbp_HM250.root
rm -rf RescorSave/RescorTables_pPb2016_Pbp_HM250
mv RescorTables_Pbp_HM250 RescorSave/RescorTables_pPb2016_Pbp_HM250
