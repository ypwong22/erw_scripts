""" Functions to read Hubbard Brook data """
import os
import pandas as pd
import numpy as np


def read_lysimeter():
    data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'DATA', 'Weathering', 
                                    'Hubbard Brook', 'knb-lter-hbr.138.11',
                                    'W1Lysim_HB1996-2020.csv'),
                    index_col = [0, 3, 2], parse_dates=True)
    data = data[['pH','Ca2+', 'Mg2+', 'Na+', 'K+', 'Alt']].sort_index() # 'Elevation',
    data[data < 0] = np.nan

    # some sites have measurements at two elevations, average them
    duplicate_sites = data.index[data.index.to_frame().reset_index(drop = True).duplicated()]

    temp = data.loc[duplicate_sites, :]
    temp = temp.groupby(temp.index).mean()
    temp.index = pd.MultiIndex.from_tuples(temp.index, names = ['Site', 'Date'])
    data.loc[temp.index, :] = temp

    data = data.drop_duplicates()

    # average over all the sites in a watershed
    data = data.groupby(data.index.get_level_values(0)).mean()

    # limit to pre-application
    # data = data.loc[data.index.get_level_values(0) < pd.Timestamp('2001-01-30'), :]

    #cation_mass = {'Ca2+': 40.078, 'Mg2+': 24.305, 'Na+': 22.99, 'K+': 39.0983, 'Alt': 26.98}
    #for col in cation_mass.keys():
    #    data[col] = data[col] # * cation_mass[col] / 1e3 # 1e-6 mol/L => g/m^3
    return data


def read_streamChem():
    # unit: umol/L

    data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'DATA', 'Weathering', 
                                    'Hubbard Brook', 'knb-lter-hbr.156.8', 'W1Long_StrmChem_HB1991-2020.csv'),
                       index_col = [0, 2], parse_dates=True)
    data = data[['pH','Ca2', 'Mg2', 'Na', 'K', 'Alt']].sort_index() # 'Elevation',
    data[data < -800] = np.nan

    # some sites have measurements at two elevations, average them
    duplicate_sites = data.index[data.index.to_frame().reset_index(drop = True).duplicated()]

    temp = data.loc[duplicate_sites, :]
    temp = temp.groupby(temp.index).mean()
    temp.index = pd.MultiIndex.from_tuples(temp.index, names = ['Site', 'Date'])
    data.loc[temp.index, :] = temp

    data = data.drop_duplicates()

    # average over all the sites in a watershed
    data = data.groupby(data.index.get_level_values(0)).mean()

    return data


def convert_streamChem(cation_conc):
    """ convert the cation flux from umol/L to g m-3 s-1 using streamflow """
    streamflow = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'DATA', 'Weathering',
                                          'Hubbard Brook', 'knb-lter-hbr.2.14', 'HBEF_DailyStreamflow_1956-2023.csv'),
                             index_col = 0, parse_dates = True)
    streamflow = streamflow.loc[streamflow['WS'] == 1, 'Streamflow']
    # streamflow.plot()
    # plt.ylabel('mm/day')

    # Streamflow rate in mm/day = 10-3 m3/day / Area of watershed (m2)
    # Concentration umol/L = 1e-3 mol/m3
    # 
    # Therefore, total transported amount is: 
    # 
    # (Concentration /1000) mol/m3 * (Streamflow rate / 1000) m3/day / Area of watershed (m2)
    #    = Concentration * Streamflow rate / 1e6 mol / day / m2
    #
    # and convert this further to mol/m2/second and multiply by cation molar mass
    cation_mass = {'Ca2': 40.078, 'Mg2': 24.305, 'Na': 22.99, 'K': 39.0983, 'Alt': 26.98}

    cations = cation_conc.drop('pH', axis = 1)
    export_rates = pd.DataFrame(np.nan, index = cations.index, columns = cations.columns)
    for col in cations.columns:
        export_rates.loc[:, col] = (streamflow.loc[cations.index] * cations[col]) \
                                    / 1e6 / 86400 * cation_mass[col]

    return export_rates