#!/bin/bash
rm Pbp_HM185.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM185_1/* > Pbp_HM185.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM185_2/* >> Pbp_HM185.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM185_3/* >> Pbp_HM185.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM185_4/* >> Pbp_HM185.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM185_5/* >> Pbp_HM185.lis
mkdir RescorTables_Pbp_HM185
root -l -b -q 'EPCalib/EPCalib.C+(285833,286009,"Pbp_HM185")'
rm offsets/offset_pPb2016_Pbp_HM185.root
mv foff_Pbp_HM185.root offsets/offset_pPb2016_Pbp_HM185.root
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_pPb2016_Pbp_HM185.db outputTag=HeavyIonRPRcd begin=285833 end=286009 infile="/rfs/sanders/EP_Pbp_HM185.root" rescor="RescorTables_Pbp_HM185"
rm /rfs/sanders/tmp_Pbp_HM185
rm save/EP_pPb2016_Pbp_HM185.root
mv /rfs/sanders/EP_Pbp_HM185.root save/EP_pPb2016_Pbp_HM185.root
rm -rf RescorSave/RescorTables_pPb2016_Pbp_HM185
mv RescorTables_Pbp_HM185 RescorSave/RescorTables_pPb2016_Pbp_HM185
