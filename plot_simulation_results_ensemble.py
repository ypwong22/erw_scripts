""" Draw spatial map of the gridded CDR, total over all the years """
import numpy as np
import pandas as pd
import os
from glob import glob
from netCDF4 import Dataset
from utils_ensemble import *
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import itertools as it
import cartopy.crs as ccrs
from utils_ensemble import  get_total_sequestration, get_pH, ens_files, load_ensemble_tuples


pft = 1
start_year = 2035
gra = 10 # 10 um, has the most number of complete series
rate_list = [2, 4, 6, 8, 10]
freq_list = [3, 9, 999]

ensemble_setups = load_ensemble_tuples()

path_out = os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble', f'pft{pft}')


################################################################################################
# Accumulated total CDR in tCO2 ha-1 util 2080

flist = ens_files(0, pft)
baseline = get_total_sequestration(flist)


with Dataset(flist[0]) as nc:
    lat = nc['lat'][:].data
    lon = nc['lon'][:].data


total_sequestration = {}

for i, (rate, freq) in enumerate(it.product(rate_list, freq_list)):

    ens_id = get_ens_id(ensemble_setups, gra, rate, start_year, freq)
    if ens_id is None or ens_id > 1099:
        continue

    flist = ens_files(ens_id, pft)
    if len(flist) != 56:
        continue

    total_sequestration[(rate,freq)] = get_total_sequestration(flist)


fig, axes = plt.subplots(5, 3, figsize = (20, 15), sharex = True, sharey = True,
                         subplot_kw={'projection': ccrs.PlateCarree()})
for i, (rate, freq) in enumerate(it.product(rate_list, freq_list)):

    ens_id = get_ens_id(ensemble_setups, gra, rate, start_year, freq)
    if ens_id is None or ens_id > 1099:
        continue

    flist = ens_files(ens_id, pft)
    if len(flist) != 56:
        continue

    ax = axes.flat[i]
    ax.set_title(f'{rate} kg/m2, Every {freq} year')

    ax.coastlines()
    ax.set_extent([-126,-66.5,23,50])
    levels = np.linspace(-0.5, 0.5, 21)
    norm = BoundaryNorm(levels, ncolors = 256, extend = 'both')
    cf = ax.pcolormesh(lon, lat, total_sequestration[(rate,freq)] - baseline, 
                       norm = norm, cmap = 'RdBu')
    plt.colorbar(cf, ax = ax, extend = 'both')

fig.suptitle('Total CDR (tCO2 ha-1) 2025-2080')
fig.savefig(os.path.join(path_out, f'map_CDR_{start_year}_{gra}.png'), dpi = 600., 
            bbox_inches = 'tight')
plt.close(fig)


################################################################################################
# pH change
baseline_pH = get_pH(flist, 2080)

ens_pH = {}

for i, (rate, freq) in enumerate(it.product(rate_list, freq_list)):

    ens_id = get_ens_id(ensemble_setups, gra, rate, start_year, freq)
    if ens_id is None or ens_id > 1099:
        continue

    flist = ens_files(ens_id, pft)
    if len(flist) != 56:
        continue

    ens_pH[(rate,freq)] = get_pH(flist, 2080)


fig, axes = plt.subplots(5, 3, figsize = (20, 15), sharex = True, sharey = True,
                         subplot_kw={'projection': ccrs.PlateCarree()})
for i, (rate, freq) in enumerate(it.product(rate_list, freq_list)):

    ens_id = get_ens_id(ensemble_setups, gra, rate, start_year, freq)
    if ens_id is None or ens_id > 1099:
        continue

    flist = ens_files(ens_id, pft)
    if len(flist) != 56:
        continue

    ax = axes.flat[i]
    ax.set_title(f'{rate} kg/m2, Every {freq} year')

    ax.coastlines()
    ax.set_extent([-126,-66.5,23,50])
    levels = np.linspace(-0.5, 0.5, 21)
    norm = BoundaryNorm(levels, ncolors = 256, extend = 'both')
    cf = ax.pcolormesh(lon, lat, ens_pH[(rate,freq)] - baseline_pH, 
                       norm = norm, cmap = 'RdBu')
    plt.colorbar(cf, ax = ax, extend = 'both')

fig.suptitle('$\Delta$ soil pH 2080')
fig.savefig(os.path.join(path_out, f'map_pH_{start_year}_{gra}.png'), dpi = 600., 
            bbox_inches = 'tight')
plt.close(fig)


################################################################################################
# Percentage dissolved primary mineral