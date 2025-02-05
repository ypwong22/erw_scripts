import os
import numpy as np
import pandas as pd
from netCDF4 import Dataset


path_root = os.path.join(os.environ['E3SM_ROOT'], 'inputdata', 'lnd', 'clm2','surfdata_map')

file_orig = os.path.join(path_root, 
                         'landuse.timeseries_0.5x0.5_combined_simyr1850-2100_c240508_newlon.nc')
file_dest = os.path.join(path_root, 
                         'landuse.conus_erw_off_combined_simyr1850-2100_c240508_newlon.nc')


nc = Dataset(file_orig, 'r')

lat = nc.variables['LATIXY'][:,0]
lon = nc.variables['LONGXY'][0,:]

lat_min_val, lat_max_val = 23.25, 54.25
lon_min_val, lon_max_val = 234.75, 293.25

lat_inds = np.where((lat >= lat_min_val) & (lat <= lat_max_val))[0]
lon_inds = np.where((lon >= lon_min_val) & (lon <= lon_max_val))[0]

lat_start, lat_end = lat_inds[0], lat_inds[-1]
lon_start, lon_end = lon_inds[0], lon_inds[-1]

print(f"Subsetting lat indices: {lat_start} to {lat_end}")
print(f"Subsetting lon indices: {lon_start} to {lon_end}")

os.system(
   f"ncks -d lsmlat,{lat_start},{lat_end} -d lsmlon,{lon_start},{lon_end} {file_orig} {file_dest}"
)

nc.close()