"""
Need to first extract the point data from CONUS data. Then modify on it.
Go to OLMT_ERW (that is: https://github.com/ypwong22/OLMT/tree/yw/elm_erw) to use the makepointdata.py script

#site=UC_Davis
site=HBR
#site=UIEF
inputgen=${E3SM_ROOT}/inputdata
python makepointdata.py --ccsm_input ${inputgen} --keep_duplicates --site ${site} --sitegroup ERW --mysimyr 1850 --use_erw
# rename the produced surfdata.pftdyn.nc to surfdata.pftdyn_erw_ctrl.nc because it is created from erw_ensemble_JRA55 ensemble 0 (control run)
"""
import os
import xarray as xr
from netCDF4 import Dataset
import numpy as np
import pandas as pd

path_root = os.path.join(os.environ['E3SM_ROOT'], 'inputdata', 'lnd', 'clm2', 'PTCLM')


################################################################################################
# UC Davis: interpolated NCSS - extracted from global data and modified
path_surffdata = os.path.join(path_root, 'UC_Davis', 'surfdata.nc')
path_newsurf = os.path.join(path_root, 'UC_Davis', 'surfdata.nc_erw_obs')

os.system(f'cp {path_surffdata} {path_newsurf}')
nc = Dataset(path_newsurf, 'r+')
nc.variables['SOIL_PH'][:] = 6.95 * np.ones(10)
nc.close()


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


##pct_kaolinite = np.array([0.04420373, 0.04420373, 0.04420373, 0.07958049, 0.07958049,
##                           0.07958049, 0.07958049, 3.9479778 , 3.9479778 , 3.9479778 ])
##pct_calcite = np.array([0.05702039, 0.05702039, 0.05702039, 0.07568372, 0.07568372,
##                        0.07568372, 0.07568372, 1.4137915 , 1.4137915 , 1.4137915 ])
##hr['PCT_KAOLINITE'] = xr.DataArray(
##    pct_kaolinite.reshape(-1, 1, 1), 
##    dims = dims, 
##    coords = coords,
##    attrs = {'long_name': 'percentage naturally occuring kaolinite in soil mineral', 
##                'units': 'g 100 g-1 soil'}
##)

##hr['PCT_CALCITE'] = xr.DataArray(
##    pct_calcite.reshape(-1, 1, 1),
##    dims = dims, 
##    coords = coords,
##    attrs = {'long_name': 'percentage naturally occuring CaCO3 in soil mineral', 
##                'units': 'g 100 g-1 soil'}
##)


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

# TOP.nc: FMAX=0.73
# TOP_test2.nc: FMAX=0.33 + SOIL COLOR
# TOP_FMAX_LOW.nc: FMAX=0.33
# TOP_FMAX_UP.nc: FMAX=0.93, preferred
hr['FMAX'] = xr.DataArray(
    [[0.93]], 
    dims = ['lsmlat', 'lsmlon'],
    attrs = {'long_name': 'maximum fractional saturated area', 
             'units': 'unitless'}
)

# 
hr['SINSL_COSAS'] = xr.DataArray(
    [[0.09012882]],
    dims=['lsmlat', 'lsmlon'],
    attrs={
        'long_name': 'sin(slope) * cos(aspect)',
        'units': 'unitless'
    }
)

hr['SINSL_SINAS'] = xr.DataArray(
    [[-0.59361121]],
    dims=['lsmlat', 'lsmlon'],
    attrs={
        'long_name': 'sin(slope) * sin(aspect)',
        'units': 'unitless'
    }
)

hr['SKY_VIEW'] = xr.DataArray(
    [[0.999356244778633]],
    dims=['lsmlat', 'lsmlon'],
    attrs={
        'long_name': 'sky view factor',
        'units': 'unitless'
    }
)

hr['STDEV_ELEV'] = xr.DataArray(
    [[72.226845371136]],
    dims=['lsmlat', 'lsmlon'],
    attrs={
        'long_name': 'standard deviation of elevation',
        'units': 'm'
    }
)

hr['TERRAIN_CONFIG'] = xr.DataArray(
    [[0.00256595338379266]],
    dims=['lsmlat', 'lsmlon'],
    attrs={
        'long_name': 'terrain configuration factor',
        'units': 'unitless'
    }
)

hr['TOPO'] = xr.DataArray(
    [[739.49707]],
    dims=['lsmlat','lsmlon'],
    attrs={
        'long_name': 'mean elevation on land', 
        'units': 'm'
    }
)

hr['STD_ELEV'] = xr.DataArray(
    [[72.226845371136]],
    dims=['lsmlat','lsmlon'],
    attrs={
        'long_name': 'standard deviation of elevation',
        'units': 'm'
    }
)

# TOP_test.nc : decrease soil albedo, not useful
##hr['SOIL_COLOR'] = xr.DataArray(
##    [[19]], 
##    dims = ['lsmlat', 'lsmlon'],
##    attrs = {'long_name': 'soil color', 
##             'units': 'unitless'}
##)

fix = {}
for v in hr.variables:
    if hr[v].dtype == np.int64:        # NetCDF-3 has no int64
        hr[v] = hr[v].astype(np.int32)
    if hr[v].dtype == bool:            # NetCDF-3 has no bool
        hr[v] = hr[v].astype(np.int8)
    # ensure _FillValue type matches variable type
    if '_FillValue' in hr[v].attrs:
        hr[v].encoding["_FillValue"] = hr[v].attrs["_FillValue"].astype(hr[v].dtype)

hr.to_netcdf(path_surffdata.replace('.nc', '_erw_TOP_FMAX_UP.nc'), format='NETCDF3_CLASSIC')
hr.close()

os.system(f'rm {path_surffdata}_temp')

################################################################################################
# UIEF
# 
# Chemistry data from Table 1 in
# Namoi, N., Lin, C.-H., Jang, C., Wasonga, D., Zumpf, C., Arshad, M. U., et al. (2025). Field-scale evaluation of ecosystem service benefits of bioenergy switchgrass. Journal of Environmental Quality, 54(3), 576–589. https://doi.org/10.1002/jeq2.70025
path_surffdata = os.path.join(path_root, 'UIEF', 'surfdata.nc')
os.system(f'cp {path_surffdata} {path_surffdata}_erw_obs')
nc = Dataset(f'{path_surffdata}_erw_obs', 'r+')

nc['ORGANIC'][:3,0,0] = 3.31 * 1.0 # SOM (%) x BD (g cm-3) converted to kg/m3
nc['ORGANIC'][3,0,0] = 3.09 * 1.1 # SOM (%) x BD (g cm-3) converted to kg/m3
nc['ORGANIC'][4,0,0] = 2.89 * 1.1 # SOM (%) x BD (g cm-3) converted to kg/m3
nc['ORGANIC'][5,0,0] = 2.56 * 1.1 # SOM (%) x BD (g cm-3) converted to kg/m3
nc['ORGANIC'][6:,0,0] = 1.87 * 1.2 # SOM (%) x BD (g cm-3) converted to kg/m3

nc['SOIL_PH'][:3,0,0] = 6.11
nc['SOIL_PH'][3,0,0] = 6.13
nc['SOIL_PH'][4,0,0] = 6.34
nc['SOIL_PH'][5,0,0] = 6.63
nc['SOIL_PH'][6:,0,0] = 6.86 # use layer 6 data for all layers below

nc['CEC_TOT'][:3,0,0] = 13.3
nc['CEC_TOT'][3,0,0] = 14.1
nc['CEC_TOT'][4,0,0] = 14.7
nc['CEC_TOT'][5,0,0] = 15.6
nc['CEC_TOT'][6:,0,0] = 14.5 # use layer 6 data for all layers below

# organic matter content is okay compared to observed

nc['CEC_ACID'][:3,0,0] = 0.2 * nc['CEC_TOT'][:3,0,0]
nc['CEC_ACID'][3:5,0,0] = 0.175 * nc['CEC_TOT'][3:5,0,0]
nc['CEC_ACID'][5:,0,0] = 0.175 * nc['CEC_TOT'][5:,0,0] # use deepest data for all the layeres below?

# scale down the CEC_EFF proportionally
# This is not too far away from the Mehlich 3 soil test Ca & Mg & K content
# reported in the paper - the test extracts exchangeable Ca + some carbonate/surface-bound Ca
# Ca = 1720 mg kg-1 => 1720 / 40 * 2 / 10 = 8.6 cmolc kg-1
# Mg = 287 mg kg-1  => 287 / 24 * 2 / 10 = 2.4 cmolc kg-1
# K = 150 mg kg-1 => 150 / 39 / 10 = 0.38 cmolc kg-1
##factor = np.empty(10)
##factor[:5] = (nc['CEC_TOT'][:5,0,0] - nc['CEC_ACID'][:5,0,0]) / \
##             (nc['CEC_EFF_1'][:5,0,0] + nc['CEC_EFF_2'][:5,0,0] + \
##              nc['CEC_EFF_3'][:5,0,0] + nc['CEC_EFF_4'][:5,0,0] + \
##              nc['CEC_EFF_5'][:5,0,0])
##factor[5:7] = nc['CEC_TOT'][5:7,0,0] / \
##    (nc['CEC_ACID'][5:7,0,0] + nc['CEC_EFF_1'][5:7,0,0] + \
##     nc['CEC_EFF_2'][5:7,0,0] + nc['CEC_EFF_3'][5:7,0,0] + \
##     nc['CEC_EFF_4'][5:7,0,0] + nc['CEC_EFF_5'][5:7,0,0])
##for i in range(1,6):
##    for j in range(7):
##        z
##for j in range(5,7):
##    nc['CEC_ACID'][j,0,0] *= factor[j]
factor = (nc['CEC_TOT'][:10,0,0] - nc['CEC_ACID'][:10,0,0]) / \
             (nc['CEC_EFF_1'][:10,0,0] + nc['CEC_EFF_2'][:10,0,0] + \
              nc['CEC_EFF_3'][:10,0,0] + nc['CEC_EFF_4'][:10,0,0] + \
              nc['CEC_EFF_5'][:10,0,0])
for i in range(1,6):
    for j in range(10):
        nc[f'CEC_EFF_{i}'][j,0,0] *= factor[j]

# re-clip CEC_ACID to prevent small numerical errors 
nc['CEC_ACID'][:,0,0] = nc['CEC_TOT'][:,0,0] - \
        ( nc['CEC_EFF_1'][:,0,0] + nc['CEC_EFF_2'][:,0,0] + \
          nc['CEC_EFF_3'][:,0,0] + nc['CEC_EFF_4'][:,0,0] + nc['CEC_EFF_5'][:,0,0] )

# soil texture data from Table 1 of
# Smith, C. M., David, M. B., Mitchell, C. A., Masters, M. D., Anderson-Teixeira, K. J., Bernacchi, C. J., & DeLucia, E. H. (2013). Reduced Nitrogen Losses after Conversion of Row Crop Agriculture to Perennial Biofuel Crops. Journal of Environmental Quality, 42(1), 219–228. https://doi.org/10.2134/jeq2012.0210
# Because Maize & Soybean in Beerling et al. 2024 ERW study, use Corn data in Table 1
nc['PCT_SAND'][:3,0,0] = 18
nc['PCT_SAND'][3:5,0,0] = 16
nc['PCT_SAND'][5,0,0] = 12
nc['PCT_SAND'][6:,0,0] = 17 # use deepest data for all the layeres below?

nc['PCT_CLAY'][:3,0,0] = 22
nc['PCT_CLAY'][3:5,0,0] = 23
nc['PCT_CLAY'][5,0,0] = 30
nc['PCT_CLAY'][6:,0,0] = 32 # use deepest data for all the layeres below?

# default LOG_KM values @20250916
# test @ 20250917
#     logkm_assumed = np.array([-4.7, -4.7, -2.825, -2.238, -5])
# test @ 20250918
logkm_assumed = np.array([-4, -4, -2.825, -1.78, -4])
cols = ['Ca', 'Mg', 'Na', 'K', 'Al']

for a, col in enumerate(cols):
    varname = f"LOG_KM_{a+1}"
    var = nc.createVariable(varname, "f4", ("nlevsoi", "lsmlat", "lsmlon"))
    var[:, :, :] = logkm_assumed[a] * np.ones((10, 1, 1))
    var.long_name = "individual Gaines-Thomas cation exchange coefficients"
    var.units = ""

nc.sync()

nc.close()