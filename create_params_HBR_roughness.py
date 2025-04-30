# Adjust the roughness length of vegetation 
import os
from netCDF4 import Dataset

os.setwd(os.path.join(os.environ['PROJDIR'], 'E3SM', 'inputdata', 'lnd', 'clm2', 'PTCLM', 'HBR'))
os.system('cp clm_params_20250427_HBR_ICB20TRCNPRDCTCBC_erw_03772_doubleZ0.nc_bak clm_params_20250427_HBR_ICB20TRCNPRDCTCBC_erw_03772_doubleZ0.nc')

nc = Dataset('clm_params_20250427_HBR_ICB20TRCNPRDCTCBC_erw_03772_doubleZ0.nc', 'r+')
# double all vegetation except grass and crop; halve crop
nc['z0mr'][:] = nc['z0mr'][:]*2
nc.sync()
nc.close()