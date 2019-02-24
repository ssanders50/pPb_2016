#!/bin/bash
conddb_import -f sqlite_file:HeavyIonRPRcd_pPb2016_pPb_MB.db -c sqlite_file:HeavyIonRPRcd_pPb2016_MB_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_pPb2016_MB_offline -b 1 -e 285832
conddb_import -f sqlite_file:HeavyIonRPRcd_pPb2016_Pbp_MB.db -c sqlite_file:HeavyIonRPRcd_pPb2016_MB_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_pPb2016_MB_offline -b 285833 -e 286496

conddb_import -f sqlite_file:HeavyIonRPRcd_pPb2016_pPb_HM120.db -c sqlite_file:HeavyIonRPRcd_pPb2016_HM120_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_pPb2016_HM120_offline -b 1 -e 285832
conddb_import -f sqlite_file:HeavyIonRPRcd_pPb2016_Pbp_HM120.db -c sqlite_file:HeavyIonRPRcd_pPb2016_HM120_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_pPb2016_HM120_offline -b 285833 -e 286496

conddb_import -f sqlite_file:HeavyIonRPRcd_pPb2016_pPb_HM150.db -c sqlite_file:HeavyIonRPRcd_pPb2016_HM150_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_pPb2016_HM150_offline -b 1 -e 285832
conddb_import -f sqlite_file:HeavyIonRPRcd_pPb2016_Pbp_HM150.db -c sqlite_file:HeavyIonRPRcd_pPb2016_HM150_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_pPb2016_HM150_offline -b 285833 -e 286496

conddb_import -f sqlite_file:HeavyIonRPRcd_pPb2016_pPb_HM185.db -c sqlite_file:HeavyIonRPRcd_pPb2016_HM185_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_pPb2016_HM185_offline -b 1 -e 285832
conddb_import -f sqlite_file:HeavyIonRPRcd_pPb2016_Pbp_HM185.db -c sqlite_file:HeavyIonRPRcd_pPb2016_HM185_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_pPb2016_HM185_offline -b 285833 -e 286496

conddb_import -f sqlite_file:HeavyIonRPRcd_pPb2016_pPb_HM250.db -c sqlite_file:HeavyIonRPRcd_pPb2016_HM250_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_pPb2016_HM250_offline -b 1 -e 285832
conddb_import -f sqlite_file:HeavyIonRPRcd_pPb2016_Pbp_HM250.db -c sqlite_file:HeavyIonRPRcd_pPb2016_HM250_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_pPb2016_HM250_offline -b 285833 -e 286496
