# Subset the CRU JRA forcing to 1980-2022 because the self calibration procedure needs a full cycle of the forcing data


cd /gpfs/wolf2/cades/cli185/proj-shared/ywo/E3SM/inputdata/atm/datm7/atm_forcing.CRUJRA_trendy_2023/cpl_bypass_full

#crujra.v2.4.5d_FLDS_1901-2022_z01.nc  crujra.v2.4.5d_PRECTmms_1901-2022_z01.nc  crujra.v2.4.5d_QBOT_1901-2022_z01.nc  crujra.v2.4.5d_UWIND_1901-2022_z01.nc
#crujra.v2.4.5d_FSDS_1901-2022_z01.nc  crujra.v2.4.5d_PSRF_1901-2022_z01.nc      crujra.v2.4.5d_TBOT_1901-2022_z01.nc  crujra.v2.4.5d_VWIND_1901-2022_z01.nc


# $ncdump -k crujra.v2.4.5d_FLDS_1901-2022_z01.nc
# 64-bit offset

# Use -6 for NetCDF3 64-bit offset
# add --ppc default=0 so NCO explicitly preserves packing/compression (short format)
ncks -O -6 \
  --map DTIME=double,LONGXY=short,LATIXY=short,FLDS=short \
  -d DTIME,28854.,44561. \
  crujra.v2.4.5d_FLDS_1901-2022_z01.nc \
  crujra.v2.4.5d_FLDS_1980-2022.nc


ncks -O -6 \
  --ppc LONGXY=0,LATIXY=0,FLDS=0 \
  -d DTIME,28854.,44561. \
  crujra.v2.4.5d_FLDS_1901-2022_z01.nc \
  crujra.v2.4.5d_FLDS_1980-2022.nc
