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

for site in ['UC_Davis']: # ['HBR']
    path_surffdata = os.path.join(path_root, site, 'surfdata.nc')

    os.system(f'cp {path_surffdata} {path_surffdata}_temp')
    hr = xr.open_dataset(f'{path_surffdata}_temp')

    dims = ['nlevsoi', 'lsmlat', 'lsmlon']
    coords = {'nlevsoi': hr['nlevsoi'], 'lsmlat': hr['lsmlat'], 'lsmlon': hr['lsmlon']}

    if site == 'UC_Davis':
        hr['SOIL_PH'] = xr.DataArray(
            np.full([nlevsoi] + [1]*(len(dims)-1), 6.95),
            dims = dims, 
            coords = coords,
            attrs = {'long_name': 'soil pH', 'units': ''}
        )
    elif 'HBR' in site:
        hr['SOIL_PH'] = xr.DataArray(
            np.full([nlevsoi] + [1]*(len(dims)-1), 4.33),
            dims = dims, 
            coords = coords,
            attrs = {'long_name': 'soil pH', 'units': ''}
        )

    hr['CEC_TOT'] = xr.DataArray(
        np.full([nlevsoi] + [1]*(len(dims)-1), np.nan),
        dims = dims, 
        coords = coords,
        attrs = {'long_name': 'total cation exchange capacity', 'units': 'meq 100g-1 dry soil'}
    )
    hr['CEC_TOT'][:3] = data_cec.loc[(site, 'layer_0'), 'CEC_tot']
    hr['CEC_TOT'][3:5] = data_cec.loc[(site, 'layer_1'), 'CEC_tot']
    hr['CEC_TOT'][5] = data_cec.loc[(site, 'layer_2'), 'CEC_tot']
    hr['CEC_TOT'][6] = data_cec.loc[(site, 'layer_3'), 'CEC_tot']
    hr['CEC_TOT'][7] = data_cec.loc[(site, 'layer_4'), 'CEC_tot']
    hr['CEC_TOT'][8:10] = data_cec.loc[(site, 'layer_5'), 'CEC_tot']

    hr['CEC_ACID'] = xr.DataArray(
        np.full([nlevsoi] + [1]*(len(dims)-1), np.nan),
        dims = dims, 
        coords = coords,
        attrs = {'long_name': 'acid cation exchange capacity', 'units': 'meq 100g-1 dry soil'}
    )
    temp = ['beta_Ca','beta_Mg','beta_Na','beta_K','beta_Al']
    hr['CEC_ACID'][:3] = (1 - data_cec.loc[(site, 'layer_0'), temp].sum()) * \
        data_cec.loc[(site, 'layer_0'), 'CEC_tot']
    hr['CEC_ACID'][3:5] = (1 - data_cec.loc[(site, 'layer_1'), temp].sum()) * \
        data_cec.loc[(site, 'layer_1'), 'CEC_tot']
    hr['CEC_ACID'][5] = (1 - data_cec.loc[(site, 'layer_2'), temp].sum()) * \
        data_cec.loc[(site, 'layer_2'), 'CEC_tot']
    hr['CEC_ACID'][6] = (1 - data_cec.loc[(site, 'layer_3'), temp].sum()) * \
        data_cec.loc[(site, 'layer_3'), 'CEC_tot']
    hr['CEC_ACID'][7] = (1 - data_cec.loc[(site, 'layer_4'), temp].sum()) * \
        data_cec.loc[(site, 'layer_4'), 'CEC_tot']
    hr['CEC_ACID'][8:10] = (1 - data_cec.loc[(site, 'layer_5'), temp].sum()) * \
        data_cec.loc[(site, 'layer_5'), 'CEC_tot']

    for a, colname in enumerate(
        ['beta_Ca','beta_Mg','beta_Na','beta_K','beta_Al']
    ):
        varname = f'CEC_EFF_{a+1}'
        hr[varname] = xr.DataArray(
            np.full([nlevsoi] + [1]*(len(dims)-1), np.nan),
            dims = dims, 
            coords = coords,
            attrs = {'long_name': 'individual cation exchange capacity',
                     'units': 'meq 100g-1 dry soil'}
        )
        hr[varname][:3] = data_cec.loc[(site, 'layer_0'), colname] * \
            data_cec.loc[(site, 'layer_0'), 'CEC_tot']
        hr[varname][3:5] = data_cec.loc[(site, 'layer_1'), colname] * \
            data_cec.loc[(site, 'layer_1'), 'CEC_tot']
        hr[varname][5] = data_cec.loc[(site, 'layer_2'), colname] * \
            data_cec.loc[(site, 'layer_2'), 'CEC_tot']
        hr[varname][6] = data_cec.loc[(site, 'layer_3'), colname] * \
            data_cec.loc[(site, 'layer_3'), 'CEC_tot']
        hr[varname][7] = data_cec.loc[(site, 'layer_4'), colname] * \
            data_cec.loc[(site, 'layer_4'), 'CEC_tot']
        hr[varname][8:10] = data_cec.loc[(site, 'layer_5'), colname] * \
            data_cec.loc[(site, 'layer_5'), 'CEC_tot']

    # checked from the global file
    if site == 'UC_Davis':
        pct_kaolinite = np.array([0.8099533, 0.8099533, 0.8099533, 5.2448173, 5.2448173, 
                                  5.2448173, 5.2448173, 2.7055457, 2.7055457, 2.7055457])
        pct_calcite = np.array([2.0673952, 2.0673952, 2.0673952, 1.1968948, 1.1968948, 1.1968948,
                                1.1968948, 0.6235169, 0.6235169, 0.6235169])
    else:
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

    hr.to_netcdf(path_surffdata.replace('.nc', '_erw.nc'))
    hr.close()

    os.system(f'rm {path_surffdata}_temp')

################################################################################################
# Hubbard Brook observation
site = 'HBR'


elm_input = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 
                                     'results', 'HBR_elm_input.csv'), index_col = 0)
logkm = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 
                                     'results', 'HBR_logkm.csv'), index_col = 0)


path_surffdata = os.path.join(path_root, site, 'surfdata.nc')
os.system(f'cp {path_surffdata} {path_surffdata}_temp')
hr = xr.open_dataset(f'{path_surffdata}_temp')

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

hr.to_netcdf(path_surffdata.replace('.nc', '_erw.nc'))
hr.close()

os.system(f'rm {path_surffdata}_temp')