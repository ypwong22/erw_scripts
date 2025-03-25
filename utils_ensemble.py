# Utilities for PFT ensembles
import os
from glob import glob
import numpy as np
from netCDF4 import Dataset
import pandas as pd


def ens_files(ens_id, pft, start = None, end = None):

    # pick the latest date if multiple directories exist
    path_ensemble = f'/gpfs/wolf2/cades/cli185/proj-shared/zdr/ERW/output/UQ/pft{pft}/'
    directories = sorted([d for d in glob(path_ensemble + f'/*{ens_id:04d}') if os.path.isdir(d)])
    directory = directories[-1]

    flist = sorted(glob(f'{directory}/run/*.h1.*.nc'))
    if '2081' in flist[-1]:
        flist = flist[:-1]

    if not start is None:
        flist2 = []
        for f in flist:
            year = int(f.split('.')[-2].split('-')[0])
            if (year >= start) & (year <= end):
                flist2.append(f)
        return flist2

    return flist


def load_ensemble_tuples():
    ensemble_setups = np.loadtxt(os.path.join(os.environ['E3SM_ROOT'], 'inputdata', 'lnd', 'clm2',
                                            'surfdata_map', 'erw_ensemble_record.txt'))
    ensemble_setups = pd.DataFrame(ensemble_setups.astype(int), 
                                columns = ['Grain size (um)', 'Rate (kg/m2)', 
                                            'Start year', 'Frequency'])
    return ensemble_setups


def get_ens_id(ensemble_setups, gra = None, rate = None, year = None, freq = None):
    filt = np.full(len(ensemble_setups), True)
    if not gra is None:
        filt = filt & (ensemble_setups['Grain size (um)'] == gra)
    if not rate is None:
        filt = filt & (ensemble_setups['Rate (kg/m2)'] == rate)
    if not year is None:
        filt = filt & (ensemble_setups['Start year'] == year)
    if not year is None:
        filt = filt & (ensemble_setups['Frequency'] == freq)

    if sum(filt) == 0:
        return None

    return ensemble_setups.index[filt][0]


def get_total_sequestration(flist):
    """ Map of total CDR until 2080 """
    total = 0
    for f in flist:
        with Dataset(f) as nc:
            # Compute: gC m-2 s-1 => gC m-2 year-1 => tCO2 ha-1 year-1
            temp = nc['r_sequestration'][0, :, :] * 365 * 24 * 3600 * 44/12 * 0.01
            total = total + temp
    total = np.where(total.mask, np.nan, total.data)
    return total


def get_pH(flist, year):
    """ Top 30cm soil pH """
    dzsoi = np.array([0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360, 0.5539, 0.9133, 
                      1.5058])

    soil_pH = None
    for f in flist:
        if str(year) in f:
            with Dataset(f) as nc:
                soil_pH = np.sum(nc['soil_pH'][0, :5, :, :] * dzsoi[:5].reshape(-1, 1, 1),
                                 axis = 0)

    if soil_pH is None:
        raise Exception('Year not in input file list.')
    return np.where(soil_pH.mask, np.nan, soil_pH.data)


def get_areaTot_sequestration(flist):
    """ Time series of accumulated CDR """
    rates = np.full(len(flist), np.nan)
    for i, f in enumerate(flist):
        with Dataset(f) as nc:
            # Compute: gC m-2 s-1 x km2 => gC year-1 => tCO2 year-1
            temp = nc['r_sequestration'][0, :, :] * nc['area'][:, :] * 1e6 * 365 * 24 * 3600 / 1e6 * 44 / 12
            # Accumulate
            rates[i] = np.ma.sum(temp)
    return rates
