import os
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from utils_ensemble import load_ensemble_tuples


# Which PFT to use
which_pft = 'pft1' # 'pft15'


# Current bad ensemble members - skip these
if which_pft == 'pft1':
    # weird values, need re-run
    skip = np.array([2241, 2245, 2249, 2253, 2257, 2261, 2264])
else:
    skip = np.array([])


def collect_predictand():
    """ Organize the simulated CDR (gC sequestered g basalt -1) into a 
        3D array (ensemble, lat, lon). Log-transform due to long-tail behavior.
    """
    nc = Dataset(os.path.join(os.environ['ZDR'], 'ERW', 'output', 'UQ', which_pft, 
                              'r_sequestration_diff.nc'))
    cdr = (nc['r_sequestration_diff'][:, :, :] - nc['n2o_diff'][:, :, :] * 270) # gC/m2/s
    cdr = np.where(cdr.mask, np.nan, cdr.data) * 86400 * 365 # gC/m2/s => gC/m2/year; note total is 56 years, 2025-2080
    pmr = nc['primary_added_diff'][:, :, :] * 86400 * 365 # g basalt/m2/s => g basalt/m2/year; note total is 56 years, 2025-2080
    nc.close()

    # delete skipped runs
    if len(skip) > 0:
        cdr[skip-1, :, :] = np.nan

    return cdr, pmr


# add ensemble-specific configurations
ensemble_setups = load_ensemble_tuples().iloc[1:, :] # skip ensemble 0
ensemble_setups.index.name = 'ensemble'


cdr, pmr = collect_predictand()
cdr = cdr.reshape(cdr.shape[0], -1)
pmr = pmr.reshape(pmr.shape[0], -1)
valid_pts = ~np.all(np.isnan(cdr), axis = 0) # remove spatial pts that does not have any data at all
cdr = cdr[:, valid_pts]
pmr = pmr[:, valid_pts]


cdr = pd.DataFrame(cdr, index = ensemble_setups.index, columns = range(cdr.shape[1]))
pmr = pd.DataFrame(pmr, index = ensemble_setups.index, columns = range(pmr.shape[1]))


# matrix of CDR
combined = pd.concat([ensemble_setups, cdr], axis = 1)
combined.to_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble', which_pft, 'combined_cdr.csv'))


combined = pd.concat([ensemble_setups, pmr], axis = 1)
combined.to_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble', which_pft, 'combined_pmr.csv'))