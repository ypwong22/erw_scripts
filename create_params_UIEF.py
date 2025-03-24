import xarray as xr
import os
import numpy as np
import shutil

########################################################################################
# Dissolution reaction studies
# --------------------------------------------------------------------------------------
# University of Illinois Energy-Farm study: 
#
# Blue Ridge Basalt, 50 t ha-1 yr-1, 2016-2019 November
#    2.6% Calcite
#    11.6% Actinolite (Amphibole): Implement as Tremolite Ca₂Mg₅Si₈O₂₂(OH)₂
#    17.8% Piemontite: Implement as Epidote
#    23.3% Albite
#    34% Chlorite: Implement as Clinochlore since it is so in USGS database
#    Ignore quartz and muscovite
########################################################################################
nminerals = 10
nminsec = 2
ncations = 5
nks = 3
string_length = 40

# Define the data
minerals_name = ['Wollastonite_CaSiO3', 'Forsterite_Mg2SiO4', 'Albite_NaAlSi3O8', 
                 'Anorthite_CaAl2Si2O8', 'Epidote_Ca2FeAl2(SiO4)3(OH)', 'Calcite_CaCO3',
                 'Tremolite_Ca2Mg5Si8O22(OH)2', 'Clinochlore14A_Mg5Al2Si3O10(OH)8',
                 'Kfeldspar_KAlSi3O8', 'Enstatite_MgSiO3']
minsecs_name = ['Calcite_CaCO3', 'Kaolinite_Al2Si2O5(OH)4']
cations_name = ['Ca2+', 'Mg2+', 'Na+', 'K+', 'Al3+']

# Fill arrays with missing values where data is not provided
fill = -9999

# --------------------------------------------------------------------------------------
# Create the dataset
#
# log_keq_primary from llnl.dat and phreeqc.dat
#
# k, e, n from Palandri et al. 2004
#
# Reaction stoichiometry following the paradigm of
#  primary mineral + proton + (water) = cations + SiO2 + (water) + (HCO3-)
#  sign convention for (water): positive if on the right side, negative if on the left
# --------------------------------------------------------------------------------------
ds = xr.Dataset({
    'minerals_name': ('nminerals', np.array(minerals_name, dtype=f'S{string_length}')),
    'minsecs_name': ('nminsecs', np.array(minsecs_name, dtype=f'S{string_length}')),
    'cations_name': ('ncations', np.array(cations_name, dtype=f'S{string_length}')),
    'log_keq_primary': ('nminerals', np.array([13.7605, 27.8626, 2.7645, 26.578, 32.9296, 
                                               1.8487, 61.2367, 67.2391, -0.2753, 11.3269])),
    'log_k_primary': (('nks', 'nminerals'), np.array([
        [ -5.37,  -6.85, -10.16,  -3.5,  -10.6,   -0.3,  -8.40, -11.11, -10.06,  -9.02], # acid
        [ -8.88, -10.64, -12.56, -9.12, -11.99,  -5.81, -10.60, -12.52, -12.41, -12.72], # neutral
        [  fill,   fill,  -15.6,  fill, -17.33,   fill,   fill,   fill,  -21.2,   fill] # base
    ])),
    'e_primary': (('nks', 'nminerals'), np.array([
        [  54.7,   67.2,     65, 16.6,    71.1,   14.4,   18.9,     88,   51.7,    80], # acid
        [  54.7,     79,   69.8, 17.8,    70.7,   23.5,   94.4,     88,     38,    80], # neutral
        [     0,      0,     71,    0,    79.1,      0,   fill,      0,   94.1,     0] # base
    ])),
    'n_primary': (('nks', 'nminerals'), np.array([
        [   0.4,   0.47,  0.457,1.411,   0.338,      1,   0.70,   0.50,    0.5,   0.6], # acid
        [     0,      0,      0,    0,       0,      0,      0,      0,      0,     0], # neutral
        [     0,      0, -0.572,    0,  -0.556,      1,      0,      0,  -0.82,     0] # base
    ])),
    'primary_stoi_proton': ('nminerals', np.array([ # one can see this from valence
        2, 4, 4, 8, 4, 1, 14, 16, 4, 2
    ])),
    'primary_stoi_cations': (('ncations', 'nminerals'), np.array([
        [1,0,0,0,0], [0,2,0,0,0], [0,0,1,0,1], [1,0,0,0,2], [2,0,0,0,0], 
        [1,0,0,0,0], [2,5,0,0,0], [0,5,0,0,2], [0,0,0,1,1],[0,1,0,0,0]
    ]).T), # one can see this from the chemical formula; except epidote
    'primary_stoi_h2o': ('nminerals', np.array([1, 2, 2, 4, 2, 0, 8, 12, 2, 1])),
    'primary_stoi_sio2': ('nminerals', np.array([1, 1, 3, 2, 3, 0, 8, 3, 3, 1])),
    'primary_stoi_hco3': ('nminerals', np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0])),
    'primary_mass': ('nminerals', np.array([116.159, 140.6931,  262.22, 278.21, 483.22, 
                                            100.0872, 812.3665, 555.7973, 278.35, 100.4])),
    'cations_mass': ('ncations', np.array([40.078, 24.305, 22.99, 39.0983, 26.98])),
    'cations_diffusivity': ('ncations', np.array([0.793e-9, 0.705e-9, 1.33e-9, 1.96e-9, 0.559e-9])),
    'bicarbonate_diffusivity': np.array([1.180e-9]),
    'carbonate_diffusivity': np.array([0.955e-9]),
    'cations_valence': ('ncations', np.array([2, 2, 1, 1, 3])),
    'minsecs_mass': ('nminsecs', np.array([100.0872, 258.1604])),
    'log_keq_minsecs': ('nminsecs', np.array([-8.48, -6.8101])),
    'alpha_minsecs': ('nminsecs', np.array([9e-10, 6.4e-14])),
})

# Set variable attributes
ds['log_k_primary'].attrs = {'long_name': 'log10 of primary mineral reaction constant at 298.15K', 
                             'unit': 'log mol m-2 s-1'} # (m-2 is the mineral surface area)
ds['e_primary'].attrs = {'long_name': 'primary mineral reaction activation energy constant at 298.15K', 'unit': 'KJ mol-1'}
ds['n_primary'].attrs = {'long_name': 'reaction order of H+ and OH- with respect to acid and basic mechanisms', 'unit': ''}
ds['log_keq_primary'].attrs = {'long_name': 'log10 of equilibrium constants for primary mineral dissolution', 'unit': ''}
ds['primary_stoi_proton'].attrs = {'long_name': 'reaction stoichiometry coefficient in front of H+', 'unit': ''}
ds['primary_stoi_h2o'].attrs = {'long_name': 'reaction stoichiometry coefficient in front of H2O (positive=right, negative=left)', 'unit': ''}
ds['primary_stoi_cations'].attrs = {'long_name': 'reaction stoichiometry coefficient in front of cations', 'unit': ''}
ds['primary_stoi_sio2'].attrs = {'long_name': 'reaction stoichiometry coefficient in front of SiO2', 'unit': ''}
ds['primary_mass'].attrs = {'long_name': 'molar mass of the primary minerals', 'unit': 'g mol-1'}
ds['cations_mass'].attrs = {'long_name': 'molar mass of the cations', 'unit': 'g mol-1'}
ds['cations_diffusivity'].attrs = {'long_name': 'diffusion coefficients of the cations in water', 'unit': 'm2 s-1'}
ds['bicarbonate_diffusivity'].attrs = {'long_name': 'diffusion coefficients of HCO3- in water', 'unit': 'm2 s-1'}
ds['carbonate_diffusivity'].attrs = {'long_name': 'diffusion coefficients of CO3-- in water', 'unit': 'm2 s-1'}
ds['cations_valence'].attrs = {'long_name': 'valence of the cations', 'unit': ''}
ds['minsecs_mass'].attrs = {'long_name': 'molar mass of the secondary minerals', 'unit': 'g mol-1'}
ds['log_keq_minsecs'].attrs = {'long_name': 'log10 of equilibrium constants for secondary mineral dissolution', 'unit': ''}
ds['log_keq_minsecs'].attrs = {'long_name': 'precipitation rate parameter of the secondary minerals', 'unit': ''}

# Set global attributes
ds.attrs['title'] = 'soil/rock powder weathering constants'
ds.attrs['Created_by'] = 'f9y,edited by ywo'
ds.attrs['Conventions'] = 'CF-1.0'
ds.attrs['Created_on'] = 'Wed Jun 05 13:12:46 EDT 2024'
ds.attrs['NCO'] = '4.6.6'
ds.attrs['history'] = 'Wed Jun 05 12:07:37 2024: created by F.-M. Yuan, ESD/CCSI-ORNL'

#encoding = {}
#for data_var in hr.data_vars:
#    if (hr[data_var].values.shape == ()) or (np.isnan(hr[data_var].values.reshape(-1)).sum() == 0):
#        encoding[data_var] = {'_FillValue': None}
#    else:
#        encoding[data_var] = {'_FillValue': 1e20}

encoding = {}
for data_var in ds.data_vars:
    if '_FillValue' in ds[data_var].encoding.keys():
        continue
    elif ds[data_var].dtype.str.startswith('|S'):
        continue
    else:
        encoding[data_var] = {'_FillValue': fill}

# Save the dataset to a NetCDF file
output_filename = os.path.join(os.environ['E3SM_ROOT'], 'inputdata', 'lnd', 'clm2', 
                               'paramdata', 'clm_erw_UIEF_params_c240718.nc')
ds.to_netcdf(output_filename, encoding = encoding)

print(f'NetCDF file {output_filename} created successfully.')