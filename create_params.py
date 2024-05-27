import xarray as xr
import os
import numpy as np
import shutil

# In the order of CaSiO3, Mg2SiO4, NaAlSi3O8, CaAl2Si2O8, Ca2FeAl2(SiO4)3(OH)
nminerals = 5
ncations = 5 # Ca Mg Na K Al

########################################################################################
# Dissolution reaction studies
# --------------------------------------------------------------------------------------
# Hubbard Brook study: Wollatonite (CaSiO3), 116.159 g/mol
# 
# CaSiO3 +2.0000 H+  =  + 1.0000 Ca++ + 1.0000 H2O + 1.0000 SiO2
# --------------------------------------------------------------------------------------
# UC Davis study:
#
# Metabasalt, 40 t ha-1, fall 2019, fall 2020
#    33.4% * 1/2 albite, NaAlSi3O8, 262.2230 g/mol
#       	NaAlSi3O8 + 4 H+ = Al+3 + Na+ + 2 H2O + 3 SiO2
#    33.4% * 1/2 anorthite, CaAl2Si2O8, 278.2073 g/mol
#           CaAl2Si2O8 + 8H+ = 2Al3+ + Ca2+ + 4H2O + 2 SiO2
#    14.3% epidote, Ca2FeAl2(SiO4)3(OH), 483.2207 g/mol
#          (suppose the formation of goethite, FeO(OH)
#           Price, J. R., Velbel, M. A., & Patino, L. C. (2005). Allanite and epidote weathering 
#              at the Coweeta Hydrologic Laboratory, western North Carolina, U.S.A. American
#              Mineralogist, 90(1), 101–114. https://doi.org/10.2138/am.2005.1444
#           )
#           Ca2FeAl2(Si2O7)(SiO4)O(OH) + 4 H+ = 2 Ca2+ + FeO(OH) + 3Al3+ + 2 H2O + 3 SiO2
#
# Olivine, 27 t ha-1, fall 2020
#    87.01% Mg2SiO4, 140.6931 g/mol
#       Mg2SiO4 +4.0000 H+  =  + 1.0000 SiO2 + 2.0000 H2O + 2.0000 Mg++
#    6.83% serpentine, not considered
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
# CaSiO3 wollastonite
hr['log_k_primary'][0,0] = -5.37
hr['log_k_primary'][1,0] = -8.88
hr['log_k_primary'][2,0] = -9999
# Mg2SiO4 forsterite
hr['log_k_primary'][0,1] = -6.85
hr['log_k_primary'][1,1] = -10.64
hr['log_k_primary'][2,1] = -9999
# NaAlSi3O8 albite (Palandri 2004 Table 13)
hr['log_k_primary'][0,2] = -10.16
hr['log_k_primary'][1,2] = -12.56
hr['log_k_primary'][2,2] = -15.60
# CaAl2Si2O8 anorthite (Palandri 2004 Table 13)
hr['log_k_primary'][0,3] = -3.50
hr['log_k_primary'][1,3] = -9.12
hr['log_k_primary'][2,3] = -9999
# Ca2FeAl2(SiO4)3(OH) epidote (Palandri 2004 Table 23)
hr['log_k_primary'][0,4] = -10.60
hr['log_k_primary'][1,4] = -11.99
hr['log_k_primary'][2,4] = -17.33

# --------------------------------------------------------------------------------------
# Activation energy
# --------------------------------------------------------------------------------------
hr['e_primary'] = xr.DataArray(np.full([3, nminerals], np.nan),
    coords = {'minerals': range(1, 1+nminerals), 'mechanism': range(1,4)}, 
    dims = ['mechanism','minerals'],
    attrs = {'long_name': 'primary mineral reaction activation energy constant at 298.15K',
             'unit': 'KJ mol-1'})
# CaSiO3 wollastonite
hr['e_primary'][0,0] = 54.7
hr['e_primary'][1,0] = 54.7
hr['e_primary'][2,0] = 0.
# Mg2SiO4 forsterite
hr['e_primary'][0,1] = 67.2
hr['e_primary'][1,1] = 79
hr['e_primary'][2,1] = 0.
# NaAlSi3O8 albite (Palandri 2004 Table 13)
hr['e_primary'][0,2] = 65.
hr['e_primary'][1,2] = 69.8
hr['e_primary'][2,2] = 71.
# CaAl2Si2O8 anorthite (Palandri 2004 Table 13)
hr['e_primary'][0,3] = 16.6
hr['e_primary'][1,3] = 17.8
hr['e_primary'][2,3] = 0
# Ca2FeAl2(SiO4)3(OH) epidote (Palandri 2004 Table 23)
hr['e_primary'][0,4] = 71.1
hr['e_primary'][1,4] = 70.7
hr['e_primary'][2,4] = 79.1


# --------------------------------------------------------------------------------------
# Reaction order on H+ and OH-
# --------------------------------------------------------------------------------------
hr['n_primary'] = xr.DataArray(np.full([3,nminerals], np.nan),
    coords = {'minerals': range(1, 1+nminerals), 'mechanism': range(1,4)}, 
    dims = ['mechanism','minerals'],
    attrs = {'long_name': 'reaction order of H+ and OH- with respect to acid and basic mechanisms', 
            'unit': ''})
# CaSiO3 wollastonite
hr['n_primary'][0,0] = 0.4
hr['n_primary'][1,0] = 0.
hr['n_primary'][2,0] = 0.
# Mg2SiO4 forsterite
hr['n_primary'][0,1] = 0.47
hr['n_primary'][1,1] = 0.
hr['n_primary'][2,1] = 0.
# NaAlSi3O8 albite (Palandri 2004 Table 13)
hr['n_primary'][0,2] = 0.457
hr['n_primary'][1,2] = 0.
hr['n_primary'][2,2] = -0.572
# CaAl2Si2O8 anorthite (Palandri 2004 Table 13)
hr['n_primary'][0,3] = 1.411
hr['n_primary'][1,3] = 0.
hr['n_primary'][2,3] = 0.
# Ca2FeAl2(SiO4)3(OH) epidote (Palandri 2004 Table 23)
hr['n_primary'][0,4] = 0.338
hr['n_primary'][1,4] = 0.
hr['n_primary'][2,4] = -0.556

# --------------------------------------------------------------------------------------
# Reaction equilibrium constants
# --------------------------------------------------------------------------------------
hr['log_keq_primary'] = xr.DataArray(np.full(nminerals, np.nan),
    coords = {'minerals': range(1, 1+nminerals)}, 
    dims = ['minerals'],
    attrs = {'long_name': 'log10 of equilibrium constants for primary mineral dissolution', 
             'unit': ''})
# CaSiO3 wollastonite (llnl.dat)
hr['log_keq_primary'][0] = 13.7605
# Mg2SiO4 forsterite (llnl.dat)
hr['log_keq_primary'][1] = 27.8626
# NaAlSi3O8 albite
hr['log_keq_primary'][2] = 2.7645
# CaAl2Si2O8 anorthite
hr['log_keq_primary'][3] = 26.5780
# Ca2FeAl2(SiO4)3(OH) epidote
hr['log_keq_primary'][4] = 32.9296


# --------------------------------------------------------------------------------------
# Reaction stoichiometry following the paradigm of
#  primary mineral + proton + (water) = cations + SiO2 + (water)
# --------------------------------------------------------------------------------------
hr['primary_stoi_proton'] = xr.DataArray(np.full(nminerals, np.nan), 
    coords = {'minerals': range(1, 1+nminerals)}, 
    dims = ['minerals'],
    attrs = {'long_name': 'reaction stoichiometry coefficient in front of H+', 
             'unit': ''})
# CaSiO3 wollastonite
hr['primary_stoi_proton'][0] = 2
# Mg2SiO4 forsterite
hr['primary_stoi_proton'][1] = 4
# NaAlSi3O8 albite
hr['primary_stoi_proton'][2] = 4
# CaAl2Si2O8 anorthite
hr['primary_stoi_proton'][3] = 8
# Ca2FeAl2(SiO4)3(OH) epidote
hr['primary_stoi_proton'][4] = 4


hr['primary_stoi_h2o_in'] = xr.DataArray(np.full(nminerals, np.nan), 
    coords = {'minerals': range(1, 1+nminerals)}, 
    dims = ['minerals'],
    attrs = {'long_name': 'reaction stoichiometry coefficient in front of H2O as an reactant',
             'unit': ''})
# CaSiO3 wollastonite
hr['primary_stoi_h2o_in'][0] = 0
# Mg2SiO4 forsterite
hr['primary_stoi_h2o_in'][1] = 0
# NaAlSi3O8 albite
hr['primary_stoi_h2o_in'][2] = 0
# CaAl2Si2O8 anorthite
hr['primary_stoi_h2o_in'][3] = 0
# Ca2FeAl2(SiO4)3(OH) epidote
hr['primary_stoi_h2o_in'][4] = 0

# dimensions will be swapped when Fortran tries to read
hr['primary_stoi_cations'] = xr.DataArray(np.full([ncations, nminerals], np.nan), 
    coords = {'minerals': range(1, 1+nminerals), 'cations': range(1, 1+ncations)}, 
    dims = ['cations','minerals'], 
    attrs = {'long_name': 'reaction stoichiometry coefficient in front of cations', 
             'unit': ''})
# CaSiO3 wollastonite
hr['primary_stoi_cations'][:, 0] = 0
hr['primary_stoi_cations'][0, 0] = 1
# Mg2SiO4 forsterite
hr['primary_stoi_cations'][:, 1] = 0
hr['primary_stoi_cations'][1, 1] = 2
# NaAlSi3O8 albite
hr['primary_stoi_cations'][:, 2] = 0
hr['primary_stoi_cations'][2, 2] = 1
hr['primary_stoi_cations'][4, 2] = 1
# CaAl2Si2O8 anorthite
hr['primary_stoi_cations'][:, 3] = 0
hr['primary_stoi_cations'][0, 3] = 1
hr['primary_stoi_cations'][4, 3] = 2
# Ca2FeAl2(SiO4)3(OH) epidote
hr['primary_stoi_cations'][:, 4] = 0
hr['primary_stoi_cations'][0, 4] = 2


hr['primary_stoi_silica'] = xr.DataArray(np.full(nminerals, np.nan), 
    coords = {'minerals': range(1, 1+nminerals)}, 
    dims = ['minerals'],
    attrs = {'long_name': 'reaction stoichiometry coefficient in front of SiO2', 'unit': ''})
# CaSiO3 wollastonite
hr['primary_stoi_silica'][0] = 1
# Mg2SiO4 forsterite
hr['primary_stoi_silica'][1] = 1
# NaAlSi3O8 albite
hr['primary_stoi_silica'][2] = 3
# CaAl2Si2O8 anorthite
hr['primary_stoi_silica'][3] = 2
# Ca2FeAl2(SiO4)3(OH) epidote
hr['primary_stoi_silica'][4] = 3



hr['primary_stoi_h2o_out'] = xr.DataArray(np.full(nminerals, np.nan), 
    coords = {'minerals': range(1, 1+nminerals)}, 
    dims = ['minerals'],
    attrs = {'long_name': 'reaction stoichiometry coefficient in front of H2O as a product', 
             'unit': ''})
# CaSiO3 wollastonite
hr['primary_stoi_h2o_out'][0] = 0
# Mg2SiO4 forsterite
hr['primary_stoi_h2o_out'][1] = 2
# NaAlSi3O8 albite
hr['primary_stoi_h2o_out'][2] = 2
# CaAl2Si2O8 anorthite
hr['primary_stoi_h2o_out'][3] = 4
# Ca2FeAl2(SiO4)3(OH) epidote
hr['primary_stoi_h2o_out'][4] = 2


# remaining solids will be calculated as being subtracted; no stoichiometry


# In the order of CaSiO3, Mg2SiO4, NaAlSi3O8, CaAl2Si2O8, Ca2FeAl2(SiO4)3(OH)
hr['primary_mass'] = xr.DataArray(np.full(nminerals, np.nan),
    coords = {'minerals': range(1, 1+nminerals)}, 
    dims = ['minerals'],
    attrs = {'long_name': 'molar mass of the primary minerals', 'unit': 'g mol-1'})
hr['primary_mass'][0] = 116.159
hr['primary_mass'][1] = 140.6931
hr['primary_mass'][2] = 262.22
hr['primary_mass'][3] = 278.21
hr['primary_mass'][4] = 483.22


#encoding = {}
#for data_var in hr.data_vars:
#    if (hr[data_var].values.shape == ()) or (np.isnan(hr[data_var].values.reshape(-1)).sum() == 0):
#        encoding[data_var] = {'_FillValue': None}
#    else:
#        encoding[data_var] = {'_FillValue': 1e20}

encoding = {}
for data_var in hr.data_vars:
    if '_FillValue' in hr[data_var].encoding.keys():
        continue
    else:
        encoding[data_var] = {'_FillValue': None}

hr.to_netcdf('clm_params.erw_20240404.nc', encoding = encoding, format='NETCDF3_CLASSIC')

hr.close()

os.remove("clm_params.erw_temp.nc")