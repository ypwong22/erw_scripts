import os
import xarray as xr
from netCDF4 import Dataset
import numpy as np
import pandas as pd

path_root = os.path.join(os.environ['E3SM_ROOT'], 'inputdata', 'lnd', 'clm2', 'PTCLM')

#data_from_gNATSGO = {
#    'UC_Davis': {'CEC_TOT': 9, 'CEC_EFF': 9, 'CEC_ACID': 1},
#    'HBR_1': {'CEC_TOT': 7, 'CEC_EFF': 7, 'CEC_ACID': 1}
#}

################################################################################################
# Interpolated NCSS - extracted from global data and modified
for site in ['UC_Davis']: 
    path_surffdata = os.path.join(path_root, site, 'surfdata_erw_from_conus.nc')
    path_newsurf = os.path.join(path_root, site, 'surfdata_erw.nc')

    os.system(f'cp {path_surffdata} {path_newsurf}')
    nc = Dataset(path_newsurf, 'r+')

    if site == 'UC_Davis':
        nc.variables['SOIL_PH'][:] = 6.95 * np.ones(10)

    nc.close()

"""
################################################################################################
# Hubbard Brook observation
elm_input = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 
                                     'results', 'HBR_elm_input.csv'), index_col = 0)
logkm = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 
                                     'results', 'HBR_logkm.csv'), index_col = 0)


path_surffdata = os.path.join(path_root, 'HBR', 'surfdata.nc')
os.system(f'cp {path_surffdata} {path_surffdata}_temp')
hr = xr.open_dataset(f'{path_surffdata}_temp')

dims = ['nlevsoi', 'lsmlat', 'lsmlon']
coords = {'nlevsoi': hr['nlevsoi'], 'lsmlat': hr['lsmlat'], 'lsmlon': hr['lsmlon']}


# Organic matter content needs all layers
# layer 9 & 10 - use C
hr['ORGANIC'] = xr.DataArray(
    np.append(
        elm_input['OM'].values, [1.256715, 1.256715]
    ).reshape(-1, 1, 1),
    dims = dims,
    coords = coords,
    attrs = {'long_name': 'organic matter density at soil levels', 
             'units': 'kg/m3 (assumed carbon content 0.58 gC per gOM)'}
)


# Below are for enhanced weathering; only need top layer
hr['SOIL_PH'] = xr.DataArray(
    np.append(
        elm_input['pH'].values, [np.nan, np.nan]
    ).reshape(-1, 1, 1),
    dims = dims,
    coords = coords,
    attrs = {'long_name': 'soil pH', 'units': ''}
)

hr['CEC_TOT'] = xr.DataArray(
    np.append(
        elm_input['CEC_TOT'].values, [np.nan, np.nan]
    ).reshape(-1, 1, 1),
    dims = dims,
    coords = coords,
    attrs = {'long_name': 'acid exchange capacity', 
             'units': 'meq 100g-1 dry soil'}
)

hr['CEC_ACID'] = xr.DataArray(
    np.append(
        elm_input['CEC_ACID'].values, [np.nan, np.nan]
    ).reshape(-1, 1, 1),
    dims = dims,
    coords = coords,
    attrs = {'long_name': 'total cation exchange capacity', 
             'units': 'meq 100g-1 dry soil'}
)

for a,col in enumerate(['CEC_Ca', 'CEC_Mg', 'CEC_Na', 'CEC_K', 'CEC_Al']):
    hr[f'CEC_EFF_{a+1}'] = xr.DataArray(
        np.append(
            elm_input[col].values, [np.nan, np.nan]
        ).reshape(-1, 1, 1),
        dims = dims, 
        coords = coords,
        attrs = {'long_name': 'individual cation exchange capacity',
                 'units': 'meq 100g-1 dry soil'}
    )


for a,col in enumerate(['Ca', 'Mg', 'Na', 'K', 'Al']):
    hr[f'LOG_KM_{a+1}'] = xr.DataArray(
        np.append(
            logkm[col].values, [np.nan, np.nan]
        ).reshape(-1, 1, 1),
        dims = dims, 
        coords = coords,
        attrs = {'long_name': 'individual Gaines-Thomas cation exchange coefficients',
                 'units': ''}
    )


pct_kaolinite = np.array([0.04420373, 0.04420373, 0.04420373, 0.07958049, 0.07958049,
                            0.07958049, 0.07958049, 3.9479778 , 3.9479778 , 3.9479778 ])
pct_calcite = np.array([0.05702039, 0.05702039, 0.05702039, 0.07568372, 0.07568372,
                        0.07568372, 0.07568372, 1.4137915 , 1.4137915 , 1.4137915 ])

hr['PCT_KAOLINITE'] = xr.DataArray(
    pct_kaolinite.reshape(-1, 1, 1), 
    dims = dims, 
    coords = coords,
    attrs = {'long_name': 'percentage naturally occuring kaolinite in soil mineral', 
                'units': 'g 100 g-1 soil'}
)

hr['PCT_CALCITE'] = xr.DataArray(
    pct_calcite.reshape(-1, 1, 1),
    dims = dims, 
    coords = coords,
    attrs = {'long_name': 'percentage naturally occuring CaCO3 in soil mineral', 
                'units': 'g 100 g-1 soil'}
)

# add variable soil thickness
hr['aveDTB'] = xr.DataArray(
    [[1.]], 
    dims = ['lsmlat', 'lsmlon'],
    attrs = {'long_name': 'mean soil depth to bedrock', 
             'units': 'meters below surface', 
             'standard_name': 'aveDTB'}
)

# increase the saturation and inundation fraction
hr['F0'] = xr.DataArray(
    [[0.12]], 
    dims = ['lsmlat', 'lsmlon'],
    attrs = {'long_name': 'maximum gridcell fractional inundated area', 
             'units': 'unitless'}
)

hr['FMAX'] = xr.DataArray(
    [[0.73]], 
    dims = ['lsmlat', 'lsmlon'],
    attrs = {'long_name': 'maximum fractional saturated area', 
             'units': 'unitless'}
)


hr.to_netcdf(path_surffdata.replace('.nc', '_erw.nc'))
hr.close()

os.system(f'rm {path_surffdata}_temp')
"""
