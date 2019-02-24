#!/bin/bash
cd EPCalib
rm *.so
rm *.d
rm *.pcm
rm -rf HiEvtPlaneList.h
rm -rf HiEvtPlaneFlatten.h
ln -s $CMSSW_BASE/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h
ln -s $CMSSW_BASE/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h
cd ..
rm -rf data/*.root
rm -rf RescorTables
rm *.db
rm /rfs/sanders/tmpsave
rm /rfs/sanders/EP.root
rm tmp.lis
#ls -1 ../CalibTree/*.root > tmp.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM120_1/* > tmp.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM120_2/* > tmp.lis
mkdir RescorTables
root -l -b -q "EPCalib/EPCalib.C+(1,285832)"
cp foff.root offsets/offset_pPb2016_HM120_1_285832.root
cd data
rm rpflat_combined.root
ln -s /rfs/sanders/EP.root rpflat_combined.root
cd ..
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_pPb2016_HM120_1_285832.db outputTag=HeavyIonRPRcd begin=1 end=285832
rm /rfs/sanders/tmpsave
rm save/EP_pPb2016_HM120_1_285832.root
mv /rfs/sanders/EP.root save/EP_pPb2016_HM120_1_285832.root
rm -rf RescorSave/RescorTables_pPb2016_HM120_1_285832
mv RescorTables RescorSave/RescorTables_pPb2016_HM120_1_285832

rm tmp.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM120_1/* > tmp.lis
ls -1  /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM120_2/* > tmp.lis
mkdir RescorTables
root -l -b -q "EPCalib/EPCalib.C+(285833,286009)"
cp foff.root offsets/offset_pPb2016_HM120_285833_286009.root
cd data
rm rpflat_combined.root
ln -s /rfs/sanders/EP.root rpflat_combined.root
cd ..
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_pPb2016_HM120_285833_286009.db outputTag=HeavyIonRPRcd begin=285833 end=286009
rm /rfs/sanders/tmpsave
rm save/EP_pPb2016_HM120_285833_286009.root
mv /rfs/sanders/EP.root save/EP_pPb2016_HM120_285833_286009.root
rm -rf RescorSave/RescorTables_pPb2016_HM120_285833_286009
mv RescorTables RescorSave/RescorTables_pPb2016_HM120_285833_286009

conddb_import -f sqlite_file:HeavyIonRPRcd_pPb2016_HM120_1_285832.db -c sqlite_file:HeavyIonRPRcd_pPb2016_HM120_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_pPb2016_MB_offline -b 1 -e 285832
conddb_import -f sqlite_file:HeavyIonRPRcd_pPb2016_HM120_285833_286009.db -c sqlite_file:HeavyIonRPRcd_pPb2016_HM120_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_pPb2016_HM120_offline -b 285833 -e 286009
