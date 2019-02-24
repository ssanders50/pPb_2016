#!/bin/bash
rm pPb_HM185.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM185_1/* > pPb_HM185.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM185_2/* >> pPb_HM185.lis
mkdir RescorTables_pPb_HM185
root -l -b -q 'EPCalib/EPCalib.C+(1,285832,"pPb_HM185")'
rm offsets/offset_pPb2016_pPb_HM185.root
mv foff_pPb_HM185.root offsets/offset_pPb2016_pPb_HM185.root
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_pPb2016_pPb_HM185.db outputTag=HeavyIonRPRcd begin=1 end=285832 infile="/rfs/sanders/EP_pPb_HM185.root" rescor="RescorTables_pPb_HM185"
rm /rfs/sanders/tmp_pPb_HM185
rm save/EP_pPb2016_pPb_HM185.root
mv /rfs/sanders/EP_pPb_HM185.root save/EP_pPb2016_pPb_HM185.root
rm -rf RescorSave/RescorTables_pPb2016_pPb_HM185
mv RescorTables_pPb_HM185 RescorSave/RescorTables_pPb2016_pPb_HM185

