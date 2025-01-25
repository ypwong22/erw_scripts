import os
import numpy as np
import pandas as pd
import itertools as it
from netCDF4 import Dataset

###########
# Settings
###########
path_root = os.path.join(os.environ['E3SM_ROOT'], 'inputdata', 'lnd', 'clm2','surfdata_map')

# Grain sizes: 10-1000 um
# Application rates: 20-100 t/ha => 2-10 kg/m2
# Day of year: Jan 1 # 244 (Sep-1 in nonleap years)
mpft = 17
planning_start = 2025
planning_end = 2050
day_of_application = 1 # 244
namendspec = 10

# minerals_name = ['Wollastonite_CaSiO3', 'Forsterite_Mg2SiO4', 'Albite_NaAlSi3O8', 
#                 'Anorthite_CaAl2Si2O8', 'Epidote_Ca2FeAl2(SiO4)3(OH)', 'Calcite_CaCO3',
#                 'Labradorite_Ca0.6Na0.4Al1.6Si2.4O8', 'Augite_Ca0.9Mg0.9Na0.1Al0.4Fe0.2Si1.9O6',
#                 'Kfeldspar_KAlSi3O8', 'Enstatite_MgSiO3']
# normal alkali basalt (fast)
pct_basalt = np.array([0, 12, 0, 0, 0, 0, 43, 21, 6, 0])
# normal tholeiitic basalt (slow)
# pct_basalt = np.array([0, 0, 0, 0, 0, 0, 45, 34, 5, 3])

# tuples = list(it.product(np.linspace(10, 1000, 10), np.linspace(2, 10, 10)))
grain_size = [2, 10, 20, 50, 100]
app_rate = [2, 4, 6, 8, 10]
app_freq = [1,3,9,999]
start_year = np.arange(planning_start, planning_end + 1)

###########
# Create a giant array that use 1/0 to indicate rock powder applied/not applied
###########
app_occ_tuples = np.full([len(start_year) * len(app_freq), 2], -1)
app_occurrence = np.full([len(start_year) * len(app_freq), 2100-1850+1], False)
count = 0
for sy, af in it.product(start_year, app_freq):
    selection = np.arange(sy, planning_end + 1, af) - 1850
    app_occurrence[count, selection] = True
    app_occ_tuples[count, :] = np.array([sy, af])
    count = count + 1

# Get the indices of the first occurrence of each unique row
app_occurrence, indices = np.unique(app_occurrence, axis=0, return_index=True)
app_occ_tuples = app_occ_tuples[indices, :]

###########
# Save the grain size, app rate, start year, and app freq to txt needed by OLMT
###########
tuples = [(1,0,-1,-1)]
for gra, app, tup in it.product(grain_size, app_rate, app_occ_tuples):
    tuples.append([gra, app, tup[0], tup[1]])
np.savetxt(os.path.join(path_root, 'erw_ensemble_record.txt'), np.array(tuples))

###########
# Create the "base" files by copying from non-ERW file
###########
skip = True
if not skip:
    file_orig = os.path.join(path_root, 
                             'landuse.timeseries_0.5x0.5_combined_simyr1850-2100_c240508_newlon.nc')
    for count in range(len(tuples) - 1):
        file_dest = os.path.join(path_root, 'erw_ensemble', f'landuse.timeseries_conus_erw_on_combined_simyr1850-2100_c240508_ensemble_{count}.nc')
        os.system(f"cp {file_orig} {file_dest}")

###########
# Loop to create netcdf files
# 4 is 1-190
###########
for count, (gra, app, tup_count) in enumerate(it.product(grain_size, app_rate, 
                                                         range(len(app_occ_tuples)))):
    if (count < ((4-1)*12)) | (count >= (4*12)):
        continue

    file_dest = os.path.join(path_root, 'erw_ensemble', f'landuse.timeseries_conus_erw_on_combined_simyr1850-2100_c240508_ensemble_{count}.nc')

    # Please make sure each time step has valid values, because the streamfile reader
    # is updating every year!

    ds = Dataset(file_dest, 'a')

    ds.createDimension('namendspec', namendspec)
    namendspec = ds.createVariable('namendspec', datatype = 'i2', dimensions = ('namendspec',))
    namendspec.long_name = 'indices of species in soil amendment mixture (e.g. rock powder)'
    namendspec.units = 'index'
    namendspec[:] = range(1, namendspec + 1)

    ds.createDimension('mpft', mpft)
    mpft_data = ds.createVariable('mpft', datatype = 'i2', dimensions = ('mpft',))
    mpft_data.long_name = 'indices of natural PFTs and cfts, if any'
    mpft_data.units = 'index'
    mpft_data[:] = range(1, mpft + 1)

    ds.sync()
    print('dimensions created')

    doy = ds.createVariable('SOIL_AMENDMENTS_DOY', datatype = 'f8', 
                            dimensions = ('time','mpft','lsmlat','lsmlon'), 
                            fill_value = 1e20)
    doy.long_name = 'soil amendment application day of year'
    doy.units = 'day of year'
    doy[:,:,:,:] = day_of_application

    ds.sync()
    print('SOIL_AMENDMENTS_DOY created')

    rate = ds.createVariable('SOIL_AMENDMENTS_RATE', datatype = 'f8', 
                            dimensions = ('time','mpft','lsmlat','lsmlon'), 
                            fill_value = 1e20)
    rate.long_name = 'soil amendment application rate'
    rate.units = 'kg/m2/yr'
    rate[:,:,:,:] = 0.
    rate[app_occurrence[tup_count, :], :, :, :] = app

    ds.sync()
    print('SOIL_AMENDMENTS_RATE created')

    size = ds.createVariable('SOIL_AMENDMENTS_GRAINSIZE', datatype = 'f8', 
                            dimensions = ('time','lsmlat','lsmlon'), 
                            fill_value = 1e20)
    size.long_name = 'grain size of applied soil amendment'
    size.units = 'micrometer'
    size[:,:,:,:] = gra

    ds.sync()
    print('SOIL_AMENDMENTS_GRAINSIZE created')

    pct = ds.createVariable('SOIL_AMENDMENTS_PCT', datatype = 'f8', 
                            dimensions = ('time', 'namendspec', 'lsmlat','lsmlon'), 
                            fill_value = 1e20)
    pct.long_name = 'species fraction of applied soil amendment'
    pct.units = 'percent'
    for i in range(namendspec):
        pct[:,i,:,:] = pct_basalt[i]

    ds.sync()
    print('SOIL_AMENDMENTS_PCT created')

    ds.close()