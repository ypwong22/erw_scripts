import xarray as xr
import os
import numpy as np
import shutil


nminerals = 2
ncations = 5
nminsec = 1 # hard-code calcite for now


########################################################################################
# Dissolution reaction studies
# --------------------------------------------------------------------------------------
# Hubbard Brook study: Wollatonite (CaSiO3), 116.159 g/mol
# 
# CaSiO3 +2.0000 H+  =  + 1.0000 Ca++ + 1.0000 H2O + 1.0000 SiO2
# --------------------------------------------------------------------------------------
# UC Davis study: olivine (Mg2SiO4), 140.6931 g/mol
# 
# Mg2SiO4 +4.0000 H+  =  + 1.0000 SiO2 + 2.0000 H2O + 2.0000 Mg++
########################################################################################

shutil.copyfile("clm_params.erw_auto.nc", "clm_params.erw_temp.nc")

hr = xr.open_dataset('clm_params.erw_temp.nc')

# --------------------------------------------------------------------------------------
# Rate constants in log10 unit
# --------------------------------------------------------------------------------------
# note: acid, neutral, and basic mechanisms are separate
# fortran dimensions are swapped in reading
hr['log_k_primary'] = xr.DataArray(np.full([3, nminerals], np.nan),
    coords = {'minerals': range(1, 1+nminerals), 'mechanism': range(1,4)}, 
    dims = ['mechanism','minerals'],
    attrs = {'long_name': 'log10 of primary mineral reaction constant at 298.15K', 
             'unit': 'log mol m-2 s-1'}) #  (m-2 is the mineral surface area)
# wollasonite
hr['log_k_primary'][0,0] = -5.37
hr['log_k_primary'][1,0] = -8.88
hr['log_k_primary'][2,0] = -9999
# olivine (forsterite)
hr['log_k_primary'][0,1] = -6.85
hr['log_k_primary'][1,1] = -10.64
hr['log_k_primary'][2,1] = -9999

# --------------------------------------------------------------------------------------
# Activation energy
# --------------------------------------------------------------------------------------
hr['e_primary'] = xr.DataArray(np.full([3, nminerals], np.nan),
    coords = {'minerals': range(1, 1+nminerals), 'mechanism': range(1,4)}, 
    dims = ['mechanism','minerals'],
    attrs = {'long_name': 'primary mineral reaction activation energy constant at 298.15K',
             'unit': 'KJ mol-1'})
# wollastonite
hr['e_primary'][0,0] = 54.7
hr['e_primary'][1,0] = 54.7
hr['e_primary'][2,0] = 0.
# olivine (forsterite)
hr['e_primary'][0,1] = 67.2
hr['e_primary'][1,1] = 79
hr['e_primary'][2,1] = 0.

# --------------------------------------------------------------------------------------
# Reaction order on H+ and OH-
# --------------------------------------------------------------------------------------
hr['n_primary'] = xr.DataArray(np.full([3,nminerals], np.nan),
    coords = {'minerals': range(1, 1+nminerals), 'mechanism': range(1,4)}, 
    dims = ['mechanism','minerals'],
    attrs = {'long_name': 'reaction order of H+ and OH- with respect to acid and basic mechanisms', 
            'unit': ''})
# wollastonite
hr['n_primary'][0,0] = 0.4
hr['n_primary'][1,0] = 0.
hr['n_primary'][2,0] = 0.
# olivine (forsterite)
hr['n_primary'][0,1] = 0.47
hr['n_primary'][1,1] = 0.
hr['n_primary'][2,1] = 0.

# --------------------------------------------------------------------------------------
# Reaction equilibrium constants
# --------------------------------------------------------------------------------------
hr['log_keq_primary'] = xr.DataArray(np.full(nminerals, np.nan),
    coords = {'minerals': range(1, 1+nminerals)}, 
    dims = ['minerals'],
    attrs = {'long_name': 'log10 of equilibrium constants for primary mineral dissolution', 
             'unit': ''})
# wollastonite
hr['log_keq_primary'][0] = 13.7605
# olivine (forsterite)
hr['log_keq_primary'][1] = 27.8626

# --------------------------------------------------------------------------------------
# Reaction stoichiometry following the paradigm of
#  primary mineral + proton + (water) = cations + SiO2 + (water)
# --------------------------------------------------------------------------------------
hr['primary_stoi_proton'] = xr.DataArray(np.full(nminerals, np.nan), 
    coords = {'minerals': range(1, 1+nminerals)}, 
    dims = ['minerals'],
    attrs = {'long_name': 'reaction stoichiometry coefficient in front of H+', 
             'unit': ''})
# wollastonite
hr['primary_stoi_proton'][0] = 2
# olivine (forsterite)
hr['primary_stoi_proton'][1] = 4


hr['primary_stoi_h2o_in'] = xr.DataArray(np.full(nminerals, np.nan), 
    coords = {'minerals': range(1, 1+nminerals)}, 
    dims = ['minerals'],
    attrs = {'long_name': 'reaction stoichiometry coefficient in front of H2O as an reactant',
             'unit': ''})
# wollastonite
hr['primary_stoi_h2o_in'][0] = 0
# olivine (forsterite)
hr['primary_stoi_h2o_in'][1] = 0

# dimensions will be swapped when Fortran tries to read
hr['primary_stoi_cations'] = xr.DataArray(np.full([ncations, nminerals], np.nan), 
    coords = {'minerals': range(1, 1+nminerals), 'cations': range(1, 1+ncations)}, 
    dims = ['cations','minerals'], 
    attrs = {'long_name': 'reaction stoichiometry coefficient in front of cations', 
             'unit': ''})
# wollastonite
hr['primary_stoi_cations'][:, 0] = 0
hr['primary_stoi_cations'][0, 0] = 1
# olivine (forsterite)
hr['primary_stoi_cations'][:, 1] = 0
hr['primary_stoi_cations'][1, 1] = 2


hr['primary_stoi_silica'] = xr.DataArray(np.full(nminerals, np.nan), 
    coords = {'minerals': range(1, 1+nminerals)}, 
    dims = ['minerals'],
    attrs = {'long_name': 'reaction stoichiometry coefficient in front of SiO2', 'unit': ''})
# wollastonite
hr['primary_stoi_silica'][0] = 1
# olivine (forsterite)
hr['primary_stoi_silica'][1] = 1


hr['primary_stoi_h2o_out'] = xr.DataArray(np.full(nminerals, np.nan), 
    coords = {'minerals': range(1, 1+nminerals)}, 
    dims = ['minerals'],
    attrs = {'long_name': 'reaction stoichiometry coefficient in front of H2O as a product', 
             'unit': ''})
# wollastonite
hr['primary_stoi_h2o_out'][0] = 0
# olivine (forsterite)
hr['primary_stoi_h2o_out'][1] = 2

# In the order of wollastonite, forsterite
hr['primary_mass'] = xr.DataArray(np.full(nminerals, np.nan),
    coords = {'minerals': range(1, 1+nminerals)}, 
    dims = ['minerals'],
    attrs = {'long_name': 'molar mass of the primary minerals', 'unit': 'g mol-1'})
hr['primary_mass'][0] = 116.159
hr['primary_mass'][1] = 140.6931

encoding = {}
for data_var in hr.data_vars:
    if '_FillValue' in hr[data_var].encoding.keys():
        continue
    else:
        encoding[data_var] = {'_FillValue': None}

hr.to_netcdf('clm_params.erw_20240203.nc')

hr.close()

os.remove("clm_params.erw_temp.nc")
