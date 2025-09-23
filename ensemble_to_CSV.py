import os
import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from utils_ensemble import load_ensemble_tuples
import pickle
import itertools as it
from tqdm import tqdm
from sklearn_quantile import RandomForestQuantileRegressor


# Which PFT to use
which_pft = 'pft1'


# Current bad ensemble members - skip these
if which_pft == 'pft1':
    skip = np.arange(2266, 2276)
else:
    skip = np.array([])


def collect_predictand():
    """ Organize the simulated CDR (gC sequestered g basalt -1) into a 
        3D array (ensemble, lat, lon). Log-transform due to long-tail behavior.
    """
    nc = Dataset(os.path.join(os.environ['ZDR'], 'ERW', 'output', 'UQ', which_pft, 
                              'r_sequestration_diff.nc'))
    cdr = (nc['r_sequestration_diff'][:, :, :] - nc['n2o_diff'][:, :, :] * 270) * 86400 # gC/m2/s
    cdr = np.where(cdr.mask, np.nan, cdr.data) * 86400 # g/m2/s
    pmr = nc['primary_added_diff'][:, :, :]
    nc.close()

    # delete skipped runs
    if len(skip) > 0:
        cdr = np.delete(cdr, skip-1, axis = 0)
        pmr = np.delete(pmr, skip-1, axis = 0)

    return cdr, pmr


# add ensemble-specific configurations
ensemble_setups = load_ensemble_tuples().iloc[1:, :] # skip ensemble 0
ensemble_setups.index.name = 'ensemble'
if len(skip) > 0:
    ensemble_setups = ensemble_setups.drop(skip, axis = 0)


cdr, pmr = collect_predictand()
cdr = cdr.reshape(cdr.shape[0], -1)
pmr = pmr.reshape(pmr.shape[0], -1)
valid_pts = ~np.any(np.isnan(cdr), axis = 0)
cdr = cdr[:, valid_pts]
pmr = pmr[:, valid_pts]


cdr = pd.DataFrame(cdr, index = ensemble_setups.index, columns = range(cdr.shape[1]))
pmr = pd.DataFrame(pmr, index = ensemble_setups.index, columns = range(pmr.shape[1]))


# matrix of CDR
combined = pd.concat([ensemble_setups, cdr], axis = 1)
combined.to_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble', which_pft, 'combined_cdr.csv'))


combined = pd.concat([ensemble_setups, pmr], axis = 1)
combined.to_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble', which_pft, 'combined_pmr.csv'))