import os
import xarray as xr
import numpy as np
import pandas as pd

path_root = os.path.join(os.environ['PROJDIR'], 'E3SM', 'inputdata', 'lnd', 'clm2', 'PTCLM')

#data_from_gNATSGO = {
#    'UC_Davis': {'CEC_TOT': 9, 'CEC_EFF': 9, 'CEC_ACID': 1},
#    'HBR_1': {'CEC_TOT': 7, 'CEC_EFF': 7, 'CEC_ACID': 1}
#}


################################################################################################
# Interpolated NCSS

data_cec = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results',
                                    'interp_NCSS_results.csv'), 
                       index_col = [0,1])
nlevsoi = 10

for site in ['UC_Davis']: # ['HBR_1', 'HBR_2']:
    path_surffdata = os.path.join(path_root, site, 'surfdata.nc')
    if not os.path.exists(f'{path_surffdata}_orig'):
        os.system(f'cp {path_surffdata} {path_surffdata}_orig')
    hr = xr.open_dataset(f'{path_surffdata}_orig')

    dims = ['nlevsoi', 'lsmlat', 'lsmlon']
    coords = {'nlevsoi': hr['nlevsoi'], 'lsmlat': hr['lsmlat'], 'lsmlon': hr['lsmlon']}

    if site == 'UC_Davis':
        hr['SOIL_PH'] = xr.DataArray(
            np.full([nlevsoi] + [1]*(len(dims)-1), 6.95),
            dims = dims, 
            coords = coords,
            attrs = {'long_name': 'soil pH', 'units': ''}
        )
    else:
        hr['SOIL_PH'] = xr.DataArray(
            np.full([nlevsoi] + [1]*(len(dims)-1), 3.5),
            dims = dims, 
            coords = coords,
            attrs = {'long_name': 'soil pH', 'units': ''}
        )

    for varname, longname, colname in \
        zip(['CEC_TOT', 'CEC_ACID'],
            ['total cation exchange capacity', 'acid exchange capacity'],
            ['cec_nh4_ph_7', 'acidity_bacl2_tea_ph_8_2']):

        hr[varname] = xr.DataArray(
            np.full([nlevsoi] + [1]*(len(dims)-1), np.nan),
            dims = dims, 
            coords = coords,
            attrs = {'long_name': longname, 'units': 'meq 100g-1 dry soil'}
        )
        hr[varname][:3] = data_cec.loc[(site, 'layer_0'), colname]
        hr[varname][3:5] = data_cec.loc[(site, 'layer_1'), colname]
        hr[varname][5] = data_cec.loc[(site, 'layer_2'), colname]
        hr[varname][6] = data_cec.loc[(site, 'layer_3'), colname]
        hr[varname][7] = data_cec.loc[(site, 'layer_4'), colname]
        hr[varname][8:10] = data_cec.loc[(site, 'layer_5'), colname]

    for a, colname in enumerate(
        ['ca_nh4_ph_7', 'mg_nh4_ph_7', 'na_nh4_ph_7', 'k_nh4_ph_7', 'aluminum_kcl_extractable']
    ):
        varname = f'CEC_EFF_{a+1}'
        hr[varname] = xr.DataArray(
            np.full([nlevsoi] + [1]*(len(dims)-1), np.nan),
            dims = dims, 
            coords = coords,
            attrs = {'long_name': 'individual cation exchange capacity',
                     'units': 'meq 100g-1 dry soil'}
        )
        hr[varname][:3] = data_cec.loc[(site, 'layer_0'), colname]
        hr[varname][3:5] = data_cec.loc[(site, 'layer_1'), colname]
        hr[varname][5] = data_cec.loc[(site, 'layer_2'), colname]
        hr[varname][6] = data_cec.loc[(site, 'layer_3'), colname]
        hr[varname][7] = data_cec.loc[(site, 'layer_4'), colname]
        hr[varname][8:10] = data_cec.loc[(site, 'layer_5'), colname]

    # Because H+ exchange capacity is not measured at the same pH as CEC_TOT and
    # includes Al3+, recalculate CEC_TOT as the sum
    hr['CEC_ACID'] = hr['CEC_ACID'] - hr['CEC_EFF_5'].values
    hr['CEC_TOT'] = hr['CEC_ACID'] + hr['CEC_EFF_1'].values + hr['CEC_EFF_2'].values + \
        hr['CEC_EFF_3'].values + hr['CEC_EFF_4'].values + hr['CEC_EFF_5'].values

    hr.to_netcdf(path_surffdata)
    hr.close()

################################################################################################
# Hubbard Brook observation
site = 'HBR_1'


elm_input = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 
                                     'results', 'HBR_elm_input.csv'), index_col = 0)
logkm = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 
                                     'results', 'HBR_logkm.csv'), index_col = 0)


path_surffdata = os.path.join(path_root, site, 'surfdata.nc')
if not os.path.exists(f'{path_surffdata}_orig'):
    os.system(f'cp {path_surffdata} {path_surffdata}_orig')
hr = xr.open_dataset(f'{path_surffdata}_orig')

dims = ['nlevsoi', 'lsmlat', 'lsmlon']
coords = {'nlevsoi': hr['nlevsoi'], 'lsmlat': hr['lsmlat'], 'lsmlon': hr['lsmlon']}


# Organic matter content needs all layers
# layer 7 - use Bs, layer 8 - use (Bh+Bsm)/2, layer 9 & 10 - use C
hr['ORGANIC'] = xr.DataArray(
    np.append(
        elm_input['OM'].values, [8.307600, 8.2178765, 1.502483, 1.502483]
    ).reshape(-1, 1, 1),
    dims = dims,
    coords = coords,
    attrs = {'long_name': 'organic matter density at soil levels', 
             'units': 'kg/m3 (assumed carbon content 0.58 gC per gOM)'}
)


# Below are for enhanced weathering; only need top layer
hr['SOIL_PH'] = xr.DataArray(
    np.append(
        elm_input['pH'].values, [np.nan, np.nan, np.nan, np.nan]
    ).reshape(-1, 1, 1),
    dims = dims,
    coords = coords,
    attrs = {'long_name': 'soil pH', 'units': ''}
)

hr['CEC_TOT'] = xr.DataArray(
    np.append(
        elm_input['CEC_TOT'].values, [np.nan, np.nan, np.nan, np.nan]
    ).reshape(-1, 1, 1),
    dims = dims,
    coords = coords,
    attrs = {'long_name': 'acid exchange capacity', 
             'units': 'meq 100g-1 dry soil'}
)

hr['CEC_ACID'] = xr.DataArray(
    np.append(
        elm_input['CEC_ACID'].values, [np.nan, np.nan, np.nan, np.nan]
    ).reshape(-1, 1, 1),
    dims = dims,
    coords = coords,
    attrs = {'long_name': 'total cation exchange capacity', 
             'units': 'meq 100g-1 dry soil'}
)

for a,col in enumerate(['CEC_Ca', 'CEC_Mg', 'CEC_Na', 'CEC_K', 'CEC_Al']):
    hr[f'CEC_EFF_{a+1}'] = xr.DataArray(
        np.append(
            elm_input[col].values, [np.nan, np.nan, np.nan, np.nan]
        ).reshape(-1, 1, 1),
        dims = dims, 
        coords = coords,
        attrs = {'long_name': 'individual cation exchange capacity',
                 'units': 'meq 100g-1 dry soil'}
    )


for a,col in enumerate(['Ca', 'Mg', 'Na', 'K', 'Al']):
    hr[f'LOG_KM_{a+1}'] = xr.DataArray(
        np.append(
            logkm[col].values, [np.nan, np.nan, np.nan, np.nan]
        ).reshape(-1, 1, 1),
        dims = dims, 
        coords = coords,
        attrs = {'long_name': 'individual Gaines-Thomas cation exchange coefficients',
                 'units': ''}
    )

hr.to_netcdf(path_surffdata)
hr.close()