#!/bin/bash
rm pPb_HM120.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM120_1/* > pPb_HM120.lis
#ls -1  /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM120_3/* >> pPb_HM120.lis
mkdir RescorTables_pPb_HM120
root -l -b -q 'EPCalib/EPCalib.C+(1,285832,"pPb_HM120")'
rm offsets/offset_pPb2016_pPb_HM120.root
mv foff_pPb_HM120.root offsets/offset_pPb2016_pPb_HM120.root
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_pPb2016_pPb_HM120.db outputTag=HeavyIonRPRcd begin=1 end=285832 infile="/rfs/sanders/EP_pPb_HM120.root" rescor="RescorTables_pPb_HM120"
rm /rfs/sanders/tmp_pPb_HM120
rm save/EP_pPb2016_pPb_HM120.root
mv /rfs/sanders/EP_pPb_HM120.root save/EP_pPb2016_pPb_HM120.root
rm -rf RescorSave/RescorTables_pPb2016_pPb_HM120
mv RescorTables_pPb_HM120 RescorSave/RescorTables_pPb2016_pPb_HM120

