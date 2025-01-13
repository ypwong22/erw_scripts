#!/bin/sh

RES=ELMMOS_USRDAT
COMPSET=IELM
MACH=pm-cpu
COMPILER=gnu
PROJECT=m3780

FORCING=gfdl-esm4
SCENARIO=ssp585

SRC_DIR=/global/homes/d/donghui/e3sm_surface_water
CASE_DIR=${SRC_DIR}/cime/scripts

cd ${SRC_DIR}

GIT_HASH=`git log -n 1 --format=%h`
CASE_NAME=GLOBE_Surface_Water_Projection_${FORCING}_${SCENARIO}_${GIT_HASH}.`date "+%Y-%m-%d-%H%M%S"`

cd ${SRC_DIR}/cime/scripts

./create_newcase -case ${CASE_DIR}/${CASE_NAME} \
-res ${RES} -mach ${MACH} -compiler ${COMPILER} -compset ${COMPSET} --project ${PROJECT}


cd ${CASE_DIR}/${CASE_NAME}

./xmlchange -file env_run.xml -id DOUT_S             -val FALSE
./xmlchange -file env_run.xml -id INFO_DBUG          -val 2

./xmlchange DATM_MODE=CLMGSWP3v1
./xmlchange CLM_USRDAT_NAME=test_r05_r05
./xmlchange LND_DOMAIN_FILE=domain_lnd_GLOBE_1d.nc
./xmlchange ATM_DOMAIN_FILE=domain_lnd_GLOBE_1d.nc
./xmlchange LND_DOMAIN_PATH=/global/cfs/projectdirs/m3780/donghui/Global-Surface-Water-Dynamics-Projection/inputdata
./xmlchange ATM_DOMAIN_PATH=/global/cfs/projectdirs/m3780/donghui/Global-Surface-Water-Dynamics-Projection/inputdata
./xmlchange CIME_OUTPUT_ROOT=/global/cfs/projectdirs/m3780/donghui/Global-Surface-Water-Dynamics-Projection/outputs

./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_END -val 2100
./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_START -val 2015
./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_ALIGN -val 1

./xmlchange DATM_CO2_TSERIES=SSP5-8.5
./xmlchange ELM_CO2_TYPE=diagnostic
./xmlchange CCSM_BGC=CO2A

./xmlchange PIO_BUFFER_SIZE_LIMIT=67108864
./xmlchange STOP_N=43,STOP_OPTION=nyears
./xmlchange JOB_QUEUE=regular
./xmlchange RESUBMIT=1
./xmlchange NTASKS=384
./xmlchange JOB_WALLCLOCK_TIME=24:00:00

./preview_namelists

cat >> user_nl_mosart << EOF
frivinp_rtm = '/global/cfs/projectdirs/m3780/donghui/Global-Surface-Water-Dynamics-Projection/inputdata/MOSART_GLOBE_1d_c230915.nc'
finidat_rtm = '/global/cfs/projectdirs/m3780/donghui/Global-Surface-Water-Dynamics-Projection/outputs/GLOBE_Surface_Water_Projection_gfdl-esm4_historical_0b7cdf4760.2024-02-22-134555.mosart.r.0065-01-01-00000.nc'
rtmhist_fincl1='RIVER_DISCHARGE_OVER_LAND_LIQ','RIVER_DISCHARGE_TO_OCEAN_LIQ','FLOODED_FRACTION','FLOODPLAIN_FRACTION'
inundflag = .true.
opt_elevprof = 1
EOF

cat >> user_nl_elm << EOF
fsurdat = '/global/cfs/projectdirs/m3780/donghui/Global-Surface-Water-Dynamics-Projection/inputdata/surfdata_GLOBE_1d_calibrated.nc'
finidat = '/global/cfs/projectdirs/m3780/donghui/Global-Surface-Water-Dynamics-Projection/outputs/GLOBE_Surface_Water_Projection_gfdl-esm4_historical_0b7cdf4760.2024-02-22-134555.elm.r.0065-01-01-00000.nc'
flanduse_timeseries = '/global/cfs/projectdirs/m3780/donghui/Global-Surface-Water-Dynamics-Projection/inputdata/landuse.timeseries_0.5x0.5_SSP5_RCP85_simyr2015-2100_GLOBAL_1d_c240209.nc'
check_dynpft_consistency = .false.
use_modified_infil = .true.
hist_empty_htapes = .true.
hist_fincl1 = 'QOVER', 'QDRAI', 'QH2OSFC', 'QRUNOFF', 'QINFL', 'QSNOMELT',    \
'FH2OSFC', 'EFLX_LH_TOT', 'RAIN', 'ZWT', 'ZWT_PERCH','FROST_TABLE','TSA',     \
'FSNO','FSAT','TWS','FSDS','FSA','FSR','FIRE','QSNOMELT','FPSN', 'FCTR', 'SNOW', \
'SOILICE','H2OSOI','H2OSFC','SOILWATER_10CM'
EOF

cat >> user_nl_datm << EOF
tintalgo = 'coszen', 'nearest', 'linear', 'linear', 'lower'
dtlimit=2.0e0,2.0e0,2.0e0,2.0e0,2.0e0
EOF

./case.setup

cp ${CASE_DIR}/${CASE_NAME}/CaseDocs/datm.streams.txt.co2tseries.SSP5-8.5 ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.co2tseries.SSP5-8.5 
chmod +rw ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.co2tseries.SSP5-8.5
perl -w -i -p -e "s@/global/cfs/cdirs/e3sm/inputdata/atm/datm7/CO2@/global/cfs/projectdirs/m3780/donghui/inputdata/CO2@" ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.co2tseries.SSP5-8.5 

# ---------------------------------------------------------------------------- #
# **************************************************************************** #
# ---------------------------------------------------------------------------- #
files1=""
for i in {2015..2100}
do
   for j in {1..12}
   do
      if [ $j -lt 10 ]
      then
         files1="${files1}clmforc.${FORCING}.${SCENARIO}.c2107.0.5x0.5.Prec.$i-0$j.nc\n"
      else
         if [ $i == 2100 ] && [ $j == 12 ]
         then
             files1="${files1}clmforc.${FORCING}.${SCENARIO}.c2107.0.5x0.5.Prec.$i-$j.nc"
         else
             files1="${files1}clmforc.${FORCING}.${SCENARIO}.c2107.0.5x0.5.Prec.$i-$j.nc\n"
         fi
      fi
   done
done

cp ${CASE_DIR}/${CASE_NAME}/CaseDocs/datm.streams.txt.CLMGSWP3v1.Precip ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Precip2
chmod +rw ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Precip2
sed '30,1060d' ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Precip2 > ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Precip
rm ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Precip2
chmod +rw ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Precip

perl -w -i -p -e "s@/global/cfs/cdirs/e3sm/inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v1.c170516/Precip@/global/cfs/projectdirs/m3780/donghui/data/forcings/${FORCING}/${SCENARIO}/Prec@" ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Precip
perl -w -i -p -e "s@/global/cfs/cdirs/e3sm/inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v1.c170516@/global/cfs/projectdirs/m3780/donghui/Runoff_Projection_Uncertainty/inputdata@" ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Precip
perl -w -i -p -e "s@domain.lnd.360x720_gswp3.0v1.c170606.nc@domain.lnd.360x720_isimip.3b.c211109.nc@" ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Precip
perl -w -i -p -e "s@clmforc.GSWP3.c2011.0.5x0.5.Prec.2015-01.nc@${files1}@" ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Precip
sed -i '/ZBOT/d' ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Precip
# ---------------------------------------------------------------------------- #
# **************************************************************************** #
# ---------------------------------------------------------------------------- #
files2=""
for i in {2015..2100}
do
   for j in {1..12}
   do
      if [ $j -lt 10 ]
      then
         files2="${files2}clmforc.${FORCING}.${SCENARIO}.c2107.0.5x0.5.Solr.$i-0$j.nc\n"
      else
         if [ $i == 2100 ] && [ $j == 12 ]
         then
             files2="${files2}clmforc.${FORCING}.${SCENARIO}.c2107.0.5x0.5.Solr.$i-$j.nc"
         else
             files2="${files2}clmforc.${FORCING}.${SCENARIO}.c2107.0.5x0.5.Solr.$i-$j.nc\n"
         fi
      fi
   done
done


cp ${CASE_DIR}/${CASE_NAME}/CaseDocs/datm.streams.txt.CLMGSWP3v1.Solar ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Solar2
chmod +rw ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Solar2
sed '30,1060d' ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Solar2 > ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Solar
rm ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Solar2
chmod +rw ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Solar

perl -w -i -p -e "s@/global/cfs/cdirs/e3sm/inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v1.c170516/Solar@/global/cfs/projectdirs/m3780/donghui/data/forcings/${FORCING}/${SCENARIO}/Solr@" ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Solar
perl -w -i -p -e "s@/global/cfs/cdirs/e3sm/inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v1.c170516@/global/cfs/projectdirs/m3780/donghui/Runoff_Projection_Uncertainty/inputdata@" ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Solar
perl -w -i -p -e "s@domain.lnd.360x720_gswp3.0v1.c170606.nc@domain.lnd.360x720_isimip.3b.c211109.nc@" ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Solar
perl -w -i -p -e "s@clmforc.GSWP3.c2011.0.5x0.5.Solr.2015-01.nc@${files2}@" ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Solar
sed -i '/ZBOT/d' ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.Solar
# ---------------------------------------------------------------------------- #
# **************************************************************************** #
# ---------------------------------------------------------------------------- #
files3=""
for i in {2015..2100}
do
   for j in {1..12}
   do
      if [ $j -lt 10 ]
      then
         files3="${files3}clmforc.${FORCING}.${SCENARIO}.c2107.0.5x0.5.TPQWL.$i-0$j.nc\n"
      else
         if [ $i == 2100 ] && [ $j == 12 ]
         then
             files3="${files3}clmforc.${FORCING}.${SCENARIO}.c2107.0.5x0.5.TPQWL.$i-$j.nc"
         else
             files3="${files3}clmforc.${FORCING}.${SCENARIO}.c2107.0.5x0.5.TPQWL.$i-$j.nc\n"
         fi
      fi
   done
done

cp ${CASE_DIR}/${CASE_NAME}/CaseDocs/datm.streams.txt.CLMGSWP3v1.TPQW ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.TPQW2
chmod +rw ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.TPQW2
sed '27d' ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.TPQW2 > ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.TPQW3
rm ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.TPQW2
chmod +rw ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.TPQW3
sed '33,1063d' ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.TPQW3 > ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.TPQW
rm ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.TPQW3
chmod +rw ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.TPQW

perl -w -i -p -e "s@/global/cfs/cdirs/e3sm/inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v1.c170516/TPHWL@/global/cfs/projectdirs/m3780/donghui/data/forcings/${FORCING}/${SCENARIO}/TPQWL@" ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.TPQW
perl -w -i -p -e "s@/global/cfs/cdirs/e3sm/inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v1.c170516@/global/cfs/projectdirs/m3780/donghui/Runoff_Projection_Uncertainty/inputdata@" ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.TPQW
perl -w -i -p -e "s@domain.lnd.360x720_gswp3.0v1.c170606.nc@domain.lnd.360x720_isimip.3b.c211109.nc@" ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.TPQW
perl -w -i -p -e "s@clmforc.GSWP3.c2011.0.5x0.5.TPQWL.2015-01.nc@${files3}@" ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.TPQW
sed -i '/ZBOT/d' ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLMGSWP3v1.TPQW
# ---------------------------------------------------------------------------- #
# **************************************************************************** #
# ---------------------------------------------------------------------------- #


./case.setup

./case.build

./case.submit
