import xarray as xr
import os
import numpy as np
import shutil

nminerals = 1
ncations = 2
nminsec = 1


shutil.copyfile("clm_params.erw_auto.nc", "clm_params.erw_temp.nc")


hr = xr.open_dataset('clm_params.erw_temp.nc')

# Wollastonite, 116.159 g/mol
# https://pubs.usgs.gov/of/2004/1068/pdf/OFR_2004_1068.pdf, page 43 & eq 11
# logk = -8.88, E = 54.7
# CaSiO3 + 2H2O + 2CO2 -> Ca2+ + 2HCO3- + H4SiO4
#
# Olivine, 140.6931 g mol-1
# 
# Mg2SiO4 + 4*CO2 + 4*H2O -> 2*Mg2+ + 4HCO3- + H4SiO4
hr['k_primary'] = xr.DataArray([10**(-8.88)],
                               dims = ['minerals'], 
                               attrs = {'long_name': 'primary mineral reaction constant at 298.15K', 'unit': 'mol m-2 s-1'}) #  (m-2 is the mineral surface area)
hr['e_primary'] = xr.DataArray([54.7], 
                               coords = {'minerals': range(1, 1+nminerals)}, 
                               dims = ['minerals'],
                               attrs = {'long_name': 'primary mineral reaction activation energy constant at 298.15K', 'unit': 'KJ mol-1'})

hr['primary_stoi_co2'] = xr.DataArray([1.], 
                                      coords = {'minerals': range(1, 1+nminerals)}, 
                                      dims = ['minerals'],
                                      attrs = {'long_name': 'reaction stoichiometry coefficient in front of CO2', 'unit': ''})
hr['primary_stoi_h2o'] = xr.DataArray([1.], 
                                      coords = {'minerals': range(1, 1+nminerals)}, 
                                      dims = ['minerals'],
                                      attrs = {'long_name': 'reaction stoichiometry coefficient in front of water', 'unit': ''})
hr['primary_stoi_cations'] = xr.DataArray([[1.,0.]], 
                                          coords = {'minerals': range(1, 1+nminerals), 'cations': range(1, 1+ncations)}, 
                                          dims = ['minerals', 'cations'], 
                                          attrs = {'long_name': 'reaction stoichiometry coefficient in front of cations', 'unit': ''})
hr['primary_stoi_bicarbonate'] = xr.DataArray([0.], 
                                              coords = {'minerals': range(1, 1+nminerals)}, 
                                              dims = ['minerals'],
                                              attrs = {'long_name': 'reaction stoichiometry coefficient in front of bicarbonate', 'unit': ''})
hr['primary_stoi_silicate'] = xr.DataArray([1.], 
                                           coords = {'minerals': range(1, 1+nminerals)}, 
                                           dims = ['minerals'],
                                           attrs = {'long_name': 'reaction stoichiometry coefficient in front of silicate', 'unit': ''})

hr['primary_weight'] = xr.DataArray([116.159], 
                                    coords = {'minerals': range(1, 1+nminerals)}, 
                                    dims = ['minerals'],
                                    attrs = {'long_name': 'formula mass of the primary minerals', 'unit': 'g mol-1'})
# Ca, Mg
hr['cation_weight'] = xr.DataArray([40.078, 24.305],
                                   coords = {'cations': range(1, 1+ncations)}, 
                                   dims = ['cations'],
                                   attrs = {'long_name': 'mass of the cations', 'unit': 'g mol-1'})
hr['cation_valence'] = xr.DataArray([2., 2.], 
                                    coords = {'cations': range(1, 1+ncations)}, 
                                    dims = ['cations'],
                                    attrs = {'long_name': 'valence of the cations', 'unit': 'g mol-1'})

#
# Calcium carbonate, 100.0869 g/mol
# llnl.dat
# log_k = -7.0017, E = 30.5767
# Ca2+ + 2*HCO3- -> CaCO3 + CO2 + H2O
hr['secondary_weight'] = xr.DataArray([100.0869], 
                                      coords = {'minsec': range(1, 1+nminsec)}, 
                                      dims = ['minsec'],
                                      attrs = {'long_name': 'formula mass of the secondary minerals', 'unit': 'g mol-1'})

hr['k_secondary'] = xr.DataArray([10**(-7.0017)], 
                                 coords = {'minsec': range(1, 1+nminsec)}, 
                                 dims = ['minsec'],
                                 attrs = {'long_name': 'secondary mineral reaction constant at 298.15K', 'unit': 'mol m-2 s-1'}) #  (m-2 is the mineral surface area)
hr['e_secondary'] = xr.DataArray([30.5767], 
                                 coords = {'minsec': range(1, 1+nminsec)}, 
                                 dims = ['minsec'],
                                 attrs = {'long_name': 'secondary mineral reaction activation energy constant at 298.15K', 'unit': 'KJ mol-1'})
hr['secondary_stoi_cations'] = xr.DataArray([[1, 0]], 
                                            coords = {'minsec': range(1, 1+nminsec), 'cations': range(1, 1+ncations)}, 
                                            dims = ['minsec', 'cations'])
hr['secondary_stoi_bicarbonate'] = xr.DataArray([2.], 
                                                coords = {'minsec': range(1, 1+nminsec)}, 
                                                dims = ['minsec'],
                                                attrs = {'long_name': 'reaction stoichiometry coefficient in front of bicarbonate', 'unit': ''})
hr['secondary_stoi_minsec'] = xr.DataArray([1.], 
                                           coords = {'minsec': range(1, 1+nminsec)}, 
                                           dims = ['minsec'],
                                           attrs = {'long_name': 'reaction stoichiometry coefficient in front of secondary mineral', 'unit': ''})
hr['secondary_stoi_co2'] = xr.DataArray([1.], 
                                        coords = {'minsec': range(1, 1+nminsec)}, 
                                        dims = ['minsec'],
                                        attrs = {'long_name': 'reaction stoichiometry coefficient in front of CO2', 'unit': ''})
hr['secondary_stoi_h2o'] = xr.DataArray([1.], 
                                        coords = {'minsec': range(1, 1+nminsec)}, 
                                        dims = ['minsec'],
                                        attrs = {'long_name': 'reaction stoichiometry coefficient in front of H2O', 'unit': ''})

encoding = {}
for data_var in hr.data_vars:
    if '_FillValue' in hr[data_var].encoding.keys():
        continue
    else:
        encoding[data_var] = {'_FillValue': None}

hr.to_netcdf('../data/clm_params.erw_20240203.nc')

hr.close()

os.remove("clm_params.erw_temp.nc")