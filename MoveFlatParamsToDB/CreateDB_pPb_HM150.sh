#!/bin/bash
rm pPb_HM150.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM150_1/* > pPb_HM150.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM150_2/* >> pPb_HM150.lis
mkdir RescorTables_pPb_HM150
root -l -b -q 'EPCalib/EPCalib.C+(1,285832,"pPb_HM150")'
rm offsets/offset_pPb2016_pPb_HM150.root
mv foff_pPb_HM150.root offsets/offset_pPb2016_pPb_HM150.root
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_pPb2016_pPb_HM150.db outputTag=HeavyIonRPRcd begin=1 end=285832 infile="/rfs/sanders/EP_pPb_HM150.root" rescor="RescorTables_pPb_HM150"
rm /rfs/sanders/tmp_pPb_HM150
rm save/EP_pPb2016_pPb_HM150.root
mv /rfs/sanders/EP_pPb_HM150.root save/EP_pPb2016_pPb_HM150.root
rm -rf RescorSave/RescorTables_pPb2016_pPb_HM150
mv RescorTables_pPb_HM150 RescorSave/RescorTables_pPb2016_pPb_HM150

