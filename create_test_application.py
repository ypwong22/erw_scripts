import xarray as xr
import os
import numpy as np

surfdata_path = os.path.join(
    os.environ['E3SM_ROOT'], 'inputdata', 'lnd', 'clm2', 'surfdata_map'
)

# edit
inputpath = os.path.join(surfdata_path, 'surfdata.pftdyn_firstday_application.nc')
outputpath = os.path.join(surfdata_path, 'surfdata.pftdyn_thirdyear_application.nc')
oldyear = 0
newyear = 3


hr = xr.open_dataset(inputpath)

hr['SOIL_AMENDMENTS_DOY'][newyear, :, :] = hr['SOIL_AMENDMENTS_DOY'][oldyear, :, :].values
hr['SOIL_AMENDMENTS_RATE'][newyear, :, :] = hr['SOIL_AMENDMENTS_RATE'][oldyear, :, :].values
hr['SOIL_AMENDMENTS_GRAINSIZE'][newyear, :] = hr['SOIL_AMENDMENTS_GRAINSIZE'][oldyear, :].values
hr['SOIL_AMENDMENTS_PCT'][newyear, :, :] = hr['SOIL_AMENDMENTS_PCT'][oldyear, :, :].values

hr['SOIL_AMENDMENTS_DOY'][oldyear, :, :] = np.nan
hr['SOIL_AMENDMENTS_RATE'][oldyear, :, :] = np.nan
hr['SOIL_AMENDMENTS_GRAINSIZE'][oldyear, :] = np.nan
hr['SOIL_AMENDMENTS_PCT'][oldyear, :, :] = np.nan

hr.to_netcdf(outputpath, format = 'NETCDF4_CLASSIC')
hr.close()