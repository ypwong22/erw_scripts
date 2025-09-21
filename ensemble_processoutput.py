from netCDF4 import Dataset
import numpy as np
import os
from concurrent.futures import ProcessPoolExecutor
from glob import glob


nens = 2276
nlon = 118
nlat = 63
startyear = 2025
endyear = 2080
pft = 'pft15'

path_root = "/gpfs/wolf2/cades/cli185/proj-shared/zdr/ERW/output/UQ/" + pft

nt = endyear - startyear + 1

def process_member(n):
    est = str(10000 + n)[1:]
    date = 20251231
    mydir = '_conus_ICB20TRCNPRDCTCBC_erw_' + pft + '_ens' + est
    while not len(glob(f"{path_root}/{date}{mydir}/run/*.h0.*.nc")) == 57 and date >= 20250101:
        date = date - 1
    print(f"{date}{mydir}")
    member_output = np.full((nt, nlat, nlon), 1e20, dtype=float)
    member_n2o = np.full((nt, nlat, nlon), 1e20, dtype=float)
    member_primary_added = np.full((nt, nlat, nlon), 1e20, dtype=float)
    lat, lon = None, None
    for y in range(startyear, endyear+1):
        myfile = f"{path_root}/{date}{mydir}/run/{date}{mydir}.elm.h0.{y}-01-01-00000.nc"
        if os.path.isfile(myfile):
            mydata = Dataset(myfile, 'r')
            if lat is None:
                lat = mydata['lat'][:]
                lon = mydata['lon'][:]

            # reset the null values to - some are 9.99999962e+35, and primary_added got sumed to 0
            temp = mydata['r_sequestration'][0, :, :]
            temp = np.where(temp.mask, 1e20, temp.data)
            member_output[y - startyear, :, :] = temp

            temp = mydata['F_N2O_DENIT'][0, :, :] + mydata['F_N2O_NIT'][0, :, :]
            temp = np.where(temp.mask, 1e20, temp.data)
            member_n2o[y - startyear, :, :] = temp

            primary_added = mydata['primary_added'][:]  # (time, mineral, lat, lon)
            temp = np.sum(primary_added[0, :, :, :], axis=0)
            temp = np.where(temp.mask, 1e20, temp.data)

            member_primary_added[y - startyear, :, :] = temp
            mydata.close()

            print(f"Finished member {n}")
    return n, member_output, member_n2o, member_primary_added, lat, lon

# Parallel processing
results = []
with ProcessPoolExecutor() as executor:
    for result in executor.map(process_member, range(nens)):
        results.append(result)

# Unpack results
output = np.full((nens, nt, nlat, nlon), 1e20, dtype=float)
n2o_output = np.full((nens, nt, nlat, nlon), 1e20, dtype=float)
primary_added_output = np.full((nens, nt, nlat, nlon), 1e20, dtype=float)
lat, lon = None, None

for n, member_output, member_n2o, member_primary_added, lat0, lon0 in results:
    output[n] = member_output
    n2o_output[n] = member_n2o
    primary_added_output[n] = member_primary_added
    if lat is None and lat0 is not None:
        lat = lat0
        lon = lon0

# Mask fill values before averaging
output_masked = np.ma.masked_values(output, 1e20)
n2o_output_masked = np.ma.masked_values(n2o_output, 1e20)
primary_added_output_masked = np.ma.masked_values(primary_added_output, 1e20)

# Compute 55-year average for each ensemble member, ignoring fill values
avg_output = np.ma.mean(output_masked, axis=1)      # shape: (nens, nlat, nlon)
avg_n2o_output = np.ma.mean(n2o_output_masked, axis=1)
avg_primary_added_output = np.ma.mean(primary_added_output_masked, axis=1)

# Subtract member 0's average from all other members (exclude member 0 in output)
diff_output = avg_output[1:, :, :] - avg_output[0, :, :]
diff_n2o_output = avg_n2o_output[1:, :, :] - avg_n2o_output[0, :, :]
diff_primary_added_output = avg_primary_added_output[1:, :, :] - avg_primary_added_output[0, :, :]

# Create coordinate arrays
ensemble = np.arange(1, nens)  # Exclude member 0

# Create a NetCDF file
ncfile = Dataset(f"{path_root}/r_sequestration_diff.nc", "w", format="NETCDF4")

# Define dimensions
ncfile.createDimension("ensemble", nens - 1)
ncfile.createDimension("lat", nlat)
ncfile.createDimension("lon", nlon)

# Create coordinate variables
ensemble_var = ncfile.createVariable("ensemble", "i4", ("ensemble",))
lat_var = ncfile.createVariable("lat", "f4", ("lat",))
lon_var = ncfile.createVariable("lon", "f4", ("lon",))

# Assign coordinate values
ensemble_var[:] = ensemble
lat_var[:] = lat
lon_var[:] = lon

# Add attributes
lat_var.units = "degrees_north"
lon_var.units = "degrees_east"

# Create the r_sequestration_diff variable
r_sequestration_diff_var = ncfile.createVariable("r_sequestration_diff", "f4", ("ensemble", "lat", "lon"), fill_value=1e20)
r_sequestration_diff_var.units = "g C/m2/s"
r_sequestration_diff_var.long_name = "Difference in 55-year avg carbon sequestration rate (member - member0)"

# Create the n2o_diff variable
n2o_diff_var = ncfile.createVariable("n2o_diff", "f4", ("ensemble", "lat", "lon"), fill_value=1e20)
n2o_diff_var.units = "g N2O-N/m2/s"
n2o_diff_var.long_name = "Difference in 55-year avg N2O flux (F_N2O_DENIT + F_N2O_NIT, member - member0)"

# Create the primary_added_diff variable
primary_added_diff_var = ncfile.createVariable("primary_added_diff", "f4", ("ensemble", "lat", "lon"), fill_value=1e20)
primary_added_diff_var.units = "g/m2/s"  # Update units as appropriate
primary_added_diff_var.long_name = "Difference in 55-year avg primary_added (summed over mineral, member - member0)"

# Assign data to the variables, filling masked values with 1e20
r_sequestration_diff_var[:] = np.ma.filled(diff_output, 1e20)
n2o_diff_var[:] = np.ma.filled(diff_n2o_output, 1e20)
primary_added_diff_var[:] = np.ma.filled(diff_primary_added_output, 1e20)

# Close the file
ncfile.close()

print("NetCDF file 'r_sequestration_diff.nc' with N2O and primary_added difference created successfully.")

