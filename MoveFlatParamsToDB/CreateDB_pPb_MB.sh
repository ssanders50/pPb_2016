#!/bin/bash
#cd EPCalib
#rm *.so
#rm *.d
#rm *.pcm
#rm -rf HiEvtPlaneList.h
#rm -rf HiEvtPlaneFlatten.h
#ln -s $CMSSW_BASE/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h
#ln -s $CMSSW_BASE/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h
#cd ..
#rm -rf data/*.root
#rm -rf RescorTables
#rm *.db
#rm /rfs/sanders/tmp_pPb_MB
#rm /rfs/sanders/EP_pPb_MB.root
rm pPb_MB.lis
ls -1  /panfs/crab_projects/crab_pPb2016_pPb_MB*/* > pPb_MB.lis
mkdir RescorTables_pPb_MB
root -l -b -q 'EPCalib/EPCalib.C+(1,285832,"pPb_MB")'
rm offsets/offset_pPb2016_pPb_MB.root
mv foff_pPb_MB.root offsets/offset_pPb2016_pPb_MB.root
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_pPb2016_pPb_MB.db outputTag=HeavyIonRPRcd begin=1 end=285832 infile="/panfs/EP_pPb_MB.root" rescor="RescorTables_pPb_MB"
rm /panfs/tmp_pPb_MB
rm save/EP_pPb2016_pPb_MB.root
mv /panfs/EP_pPb_MB.root save/EP_pPb2016_pPb_MB.root
rm -rf RescorSave/RescorTables_pPb2016_pPb_MB
mv RescorTables_pPb_MB RescorSave/RescorTables_pPb2016_pPb_MB

