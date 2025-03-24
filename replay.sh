#!/bin/bash

set -e

# Created 2024-10-03 18:13:33
# Modified 2024-10-15 18:47:33

CASEDIR="/gpfs/wolf2/cades/cli185/proj-shared/ywo/E3SM/case_dirs/20241003_conus_ICB20TRCNPRDCTCBC_erw_pft1"

/autofs/nccsopen-svm1_home/ywo/models/E3SM_ERW/cime/scripts/create_newcase --case "${CASEDIR}" --mach cades-baseline --compset ICB20TRCNPRDCTCBC --res hcru_hcru --walltime 24:0:00 --handle-preexisting-dirs u --project CLI185

cd "${CASEDIR}"

./xmlchange SAVE_TIMING=FALSE

./xmlchange EXEROOT=/gpfs/wolf2/cades/cli185/proj-shared/ywo/ERW/output/20241003_conus_ICB20TRCNPRDCTCBC_erw_pft1/bld

./xmlchange PIO_VERSION=2

./xmlchange MOSART_MODE=NULL

./xmlchange RUNDIR=/gpfs/wolf2/cades/cli185/proj-shared/ywo/ERW/output/20241003_conus_ICB20TRCNPRDCTCBC_erw_pft1/run

./xmlchange DIN_LOC_ROOT=/gpfs/wolf2/cades/cli185/proj-shared/ywo/E3SM/inputdata

./xmlchange DIN_LOC_ROOT_CLMFORC=/gpfs/wolf2/cades/cli185/proj-shared/ywo/E3SM/inputdata/atm/datm7/

./xmlchange RUN_STARTDATE=1-01-01

./xmlchange DOUT_S=FALSE

./xmlchange ATM_NCPL=24

./xmlchange RUN_REFDATE=0401-01-01

./xmlchange CCSM_BGC=CO2A

./xmlchange ELM_CO2_TYPE=diagnostic

./xmlchange NTASKS_ATM=512

./xmlchange NTHRDS_ATM=1

./xmlchange NTASKS_LND=512

./xmlchange NTHRDS_LND=1

./xmlchange NTASKS_ICE=512

./xmlchange NTHRDS_ICE=1

./xmlchange NTASKS_OCN=512

./xmlchange NTHRDS_OCN=1

./xmlchange NTASKS_CPL=512

./xmlchange NTHRDS_CPL=1

./xmlchange NTASKS_GLC=512

./xmlchange NTHRDS_GLC=1

./xmlchange NTASKS_ROF=512

./xmlchange NTHRDS_ROF=1

./xmlchange NTASKS_WAV=512

./xmlchange NTHRDS_WAV=1

./xmlchange NTASKS_ESP=512

./xmlchange NTHRDS_ESP=1

./xmlchange NTASKS_IAC=512

./xmlchange NTHRDS_IAC=1

./xmlchange STOP_OPTION=nyears

./xmlchange STOP_N=150

./xmlchange REST_N=150

./xmlchange REST_N=20

./xmlchange PIO_TYPENAME=netcdf

./case.setup

./xmlchange ATM_DOMAIN_PATH='${RUNDIR}'

./xmlchange LND_DOMAIN_PATH='${RUNDIR}'

./xmlchange ATM_DOMAIN_FILE=domain.nc

./xmlchange LND_DOMAIN_FILE=domain.nc

./xmlchange BUILD_COMPLETE=TRUE

./preview_namelists

./preview_namelists

