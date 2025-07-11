from netCDF4 import Dataset
import numpy as np
import os
from concurrent.futures import ProcessPoolExecutor
from utils import vert_interp


def vert_interp2(data):
    """ input_data: (levgrnd, lat, lon) """
    elm_nodes = np.array([0.71, 2.79, 6.23, 11.89, 21.22, 36.61, 61.98, 103.80])
    elm_interface = np.array([0, 1.75, 4.51, 9.06, 16.55, 28.91, 49.29, 82.89, 138.28])

    data = data[:len(elm_nodes), :, :]
    data = data.reshape((data.shape[0], data.shape[1] * data.shape[2])).T
    data_new = vert_interp(np.array([2.5, 5, 15, 50]), elm_nodes,
                          data, False,
                          np.array([0, 5, 10, 30, 100]), elm_interface)
    return data_new


path_root = "/gpfs/wolf2/cades/cli185/proj-shared/zdr"


def process_met():
    """ identical for all members """
    member_output = {}

    # meteorological driver
    # z01: [-180, -90], z02: [-90, 0], z03: [0, 90], z04: [90, 180]
    for var in ['TBOT','QBOT','WIND','FSDS','PRECTmms']:
        myfile = f"{path_root}/atm_forcing.ISIMIP.DonghuiXu.2024/cpl_bypass_full/gfdl-esm4.historical.c2107.0.5x0.5_{var}_1951-2100_z01.nc"
        mydata = Dataset(myfile, 'r')
        lon_1 = mydata['LONGXY'][:]
        lat_1 = mydata['LATIXY'][:]
        sect_1 = mydata[var][:, (365*(2025-1951)):(365*(2080-1951))].mean(axis = 1)
        mydata.close()

        myfile = f"{path_root}/atm_forcing.ISIMIP.DonghuiXu.2024/cpl_bypass_full/gfdl-esm4.historical.c2107.0.5x0.5_{var}_1951-2100_z02.nc"
        mydata = Dataset(myfile, 'r')
        lon_2 = mydata['LONGXY'][:]
        lat_2 = mydata['LATIXY'][:]
        sect_2 = mydata[var][:, (365*(2025-1951)):(365*(2080-1951))].mean(axis = 1)
        mydata.close()

        lon = np.concatenate([lon_1, lon_2]) + 360
        lat = np.concatenate([lat_1, lat_2])
        sect = np.concatenate([sect_1, sect_2])

        subset = (lon >= 234.75) & (lon <= 293.25) & (lat >= 23.25) & (lat <= 54.25)
        lon = lon[subset]
        lat = lat[subset]
        sect = sect[subset].reshape(len(lon[::len(np.unique(lat))]), -1)

        member_output[var] = sect.T # re-adjust to lat-lon order

    # reverse latitude
    lat = lat[::-1]
    for key in member_output:
        member_output[key] = member_output[key][::-1, :]

    return member_output, lat, lon


def process_surf():
    """ identical for all members """
    member_output = {}
    myfile = f"{path_root}/ERW/output/UQ/pft1/20250613_conus_ICB20TRCNPRDCTCBC_erw_pft1_ens0060/run/surfdata.nc"
    mydata = Dataset(myfile, 'r')
    for var in ['ORGANIC','PCT_CLAY','PCT_SAND','SOIL_PH','CEC_TOT','CEC_ACID',
                'CEC_EFF_1','CEC_EFF_2','CEC_EFF_3','CEC_EFF_4','CEC_EFF_5']:
        output = vert_interp2(mydata[var])
        for i, d in enumerate([5, 10, 50, 100]):
            member_output[f'{var}_{d}'] = output[:, i].reshape((nlat, nlon))
    member_output['SLOPE'] = mydata['SLOPE'][:,:]
    member_output['STD_ELEV'] = mydata['STD_ELEV'][:,:]
    lat = mydata['LATIXY'][:, 0]
    lon = mydata['LONGXY'][0, :]
    mydata.close()

    # fix range
    for var in ['SOIL_PH','CEC_TOT','CEC_ACID',
                'CEC_EFF_1','CEC_EFF_2','CEC_EFF_3','CEC_EFF_4','CEC_EFF_5']:
        for i, d in enumerate([5, 10, 50, 100]):
            member_output[f'{var}_{d}'] = np.where(member_output[f'{var}_{d}'] > 400, 1e36, 
                                                   member_output[f'{var}_{d}'])

    return member_output, lat, lon


nens = 2267
nlon = 118
nlat = 63
startyear = 2025
endyear = 2080

nt = endyear - startyear + 1


def process_member(n):
    member_output = {}

    est = str(10000 + n)[1:]
    date = 20251231
    mydir = '_conus_ICB20TRCNPRDCTCBC_erw_pft1_ens' + est
    while not os.path.isdir(f"{path_root}/ERW/output/UQ/pft1/{date}{mydir}") and date >= 20250101:
        date = date - 1

    # qrunoff, qdrai, soil moisture 5cm, 10cm, 30cm, 100cm
    member_output['QRUNOFF'] = np.full((nt, nlat, nlon), np.nan, dtype=float)
    member_output['QDRAI'] = np.full((nt, nlat, nlon), np.nan, dtype=float)
    for i, d in enumerate([5, 10, 50, 100]):
        member_output[f'SM_{d}'] = np.full((nt, nlat, nlon), np.nan, dtype=float)
    lat, lon = None, None
    for y in range(startyear, endyear+1):
        myfile = f"{path_root}/ERW/output/UQ/pft1/{date}{mydir}/run/{date}{mydir}.elm.h0.{y}-01-01-00000.nc"
        if os.path.isfile(myfile):
            mydata = Dataset(myfile, 'r')
            if lat is None:
                lat = mydata['lat'][:]
                lon = mydata['lon'][:]
            member_output['QRUNOFF'][y - startyear, :, :] = np.where(
                mydata['QRUNOFF'][0, :, :].mask, np.nan, mydata['QRUNOFF'][0, :, :].data)
            member_output['QDRAI'][y - startyear, :, :] = np.where(
                mydata['QDRAI'][0, :, :].mask, np.nan, mydata['QRUNOFF'][0, :, :].data)

            output = vert_interp2(np.where(
                mydata['H2OSOI'][0, :, :, :].mask, np.nan, mydata['H2OSOI'][0, :, :, :].data))
            for i, d in enumerate([5, 10, 50, 100]):
                member_output[f'SM_{d}'][y - startyear, :, :] = output[:,i].reshape((nlat, nlon))

            mydata.close()

    # Average over time
    for key,item in member_output.items():
        member_output[key] = np.nanmean(item, axis=0)
        member_output[key] = np.where(np.isnan(member_output[key]), 1e36, member_output[key])

    print(f"Finished member {n}")
    return n, member_output, lat, lon


# Parallel processing
met_result, lat1, lon1 = process_met()
print(lat1)
print(lon1)

surf_result, lat2, lon2 = process_surf()
print(lat2)
print(lon2)

# n, member_out, lat2, lon2 = process_member(1)

results = []
with ProcessPoolExecutor() as executor:
    for result in executor.map(process_member, range(nens)):
        results.append(result)

# Unpack results
output = {}
lat3, lon3 = None, None

for n, member_output, lat0, lon0 in results:
    if n == 0:
        for key in member_output:
            output[key] = np.full((nens, nlat, nlon), np.nan, dtype=float)
    output[key][n, :, :] = member_output[key]
    if lat3 is None and lat0 is not None:
        lat3 = lat0
        lon3 = lon0
print(lat3)
print(lon3)

# Mask fill values
for key in member_output:
    output[key] = np.ma.masked_values(output[key], 1e36)

# Create coordinate arrays
ensemble = np.arange(nens)

# Create a NetCDF file
ncfile = Dataset(f"{path_root}/ERW/output/UQ/pft1/predictors.nc", "w", format="NETCDF4")

# Define dimensions
ncfile.createDimension("ensemble", nens)
ncfile.createDimension("lat", nlat)
ncfile.createDimension("lon", nlon)

# Create coordinate variables
ensemble_var = ncfile.createVariable("ensemble", "i4", ("ensemble",))
lat_var = ncfile.createVariable("lat", "f4", ("lat",))
lon_var = ncfile.createVariable("lon", "f4", ("lon",))

# Assign coordinate values
ensemble_var[:] = ensemble
lat_var[:] = lat3
lon_var[:] = lon3

# Add attributes
lat_var.units = "degrees_north"
lon_var.units = "degrees_east"

# Create the meteorological variables
for met, metdata in met_result.items():
    f = ncfile.createVariable(met, "f4", ("lat", "lon"), fill_value = 1e36)
    if met == "TBOT":
        f.units = "K"
    elif met == "WIND":
        f.units = "m/s"
    elif met == "PRECTmms":
        f.units = "mm/s"
    else:
        f.units = ""
    f[:,:] = np.ma.filled(metdata, 1e36)


# Create the surfdata variables
for key, surfdata in surf_result.items():
    f = ncfile.createVariable(key, "f4", ("lat", "lon"), fill_value = 1e36)
    if key == "ORGANIC":
        f.units = "kg/m3"
    elif key in ["PCT_SAND","PCT_CLAY"]:
        f.units = "%"
    elif "CEC" in key:
        f.units = "meq 100g-1 dry soil"
    else:
        f.units = ""
    f[:,:] = np.ma.filled(surfdata, 1e36)


# Create the ensemble variables
for key, keydata in member_output.items():
    f = ncfile.createVariable(key, "f4", ("ensemble", "lat", "lon"), fill_value=1e36)
    if key in ["QRUNOFF", "QDRAI"]:
        f.units = "mm/s"
    else:
        f.units = "m3/m3"
    f[:,:,:] = np.ma.masked_values(keydata, 1e36)


# Close the file
ncfile.close()

print("NetCDF file 'predictors.nc' created successfully.")