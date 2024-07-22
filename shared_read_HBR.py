""" Functions to read Hubbard Brook data """
import os
import pandas as pd
import numpy as np
from scipy.stats import linregress


def read_lai():
    """ Leaf area index is estimated from litterfall for Watershed 1"""
    data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 
                                    'Hubbard_Brook', 'knb-lter-hbr.293.2',
                                    'HBEF_WS1_LAI_1998-2019.csv'),
                       index_col = [1, 0, 2])
    data[data < -900] = np.nan
    return data


def read_pheno():
    """ Phenology is measured from 0-4 as no leaf-out to summer leaf condition"""
    data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 
                                    'Hubbard_Brook', 'knb-lter-hbr.51.14',
                                    'HBEF_Phenology_longform.csv'),
                       index_col = 1, parse_dates=True)
    data = data.loc[data['SEASON'] == 'SPRING', :]

    # Create a linear regression line and uncertainty intervals
    x = data['DAY']
    y = data['Phenology_Stage']
    res = linregress(x, y)
    return res.slope, res.intercept, res.stderr


def read_snowcourse():
    """ Snowcourse 2 is in Watershed 1 and has continued measurements till present """
    data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 
                                    'Hubbard_Brook', 'knb-lter-hbr.27.20',
                                    'HBEF_snowcourse_1956-2024.csv'),
                       index_col = [2, 1], parse_dates=True)
    data = data.loc['STA2', 'swe']
    return data


def read_lysimeter():
    data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 
                                    'Hubbard_Brook', 'knb-lter-hbr.138.11',
                                    'W1Lysim_HB1996-2020.csv'),
                    index_col = [3, 2, 0], parse_dates=True)
    data = data[['pH','Ca2+', 'Mg2+', 'Na+', 'K+', 'Alt']].sort_index() # 'Elevation',
    data[data < 0] = np.nan

    # umol/L => mol/L
    data.iloc[:, 1:] = data.iloc[:, 1:] * 1e-6

    # some sites have measurements at two elevations, average them
    duplicate_sites = data.index[data.index.to_frame().reset_index(drop = True).duplicated()]

    temp = data.loc[duplicate_sites, :]
    temp = temp.groupby(temp.index).mean()
    temp.index = pd.MultiIndex.from_tuples(temp.index, names = ['Site', 'Date'])
    data.loc[temp.index, :] = temp

    data = data.drop_duplicates().reset_index()

    # The Oe horizon is only for Site 13, the Oa horizon only for Site 1-12
    # The Oae horizon only for Site 9, and is pre-2014 data; the post-2014 data for Site 9
    # only has the Oa horizon
    #
    # The Bh horizon is for Sites 3-7 & 10 & 12; The Bhs horizon is for Sites 8-9
    # The Bhs1 & Bhs3 horizons are for Sites 11-13
    #
    # The Bs horizon is for Sites 1-7 & 10 & 12; The Bs2 horizon is for Sites 8-9
    data['Horizon'] = data['Horizon'].map({
        'Oa': 'Oae', 'Oe': 'Oae', 'Oae': 'Oae', 
        'Bh': 'Bhs', 'Bhs': 'Bhs', 'Bhs1': 'Bhs', 'Bhs3': 'Bhs', 
        'Bs': 'Bhs', 'Bs2': 'Bhs'
    })

    data = data.set_index(['Horizon', 'Site', 'Date'])

    data_mean = data.groupby(['Horizon', 'Date']).mean()
    data_std = data.groupby(['Horizon', 'Date']).std()

    # limit to pre-application
    # data = data.loc[data.index.get_level_values(0) < pd.Timestamp('2001-01-30'), :]

    #cation_mass = {'Ca2+': 40.078, 'Mg2+': 24.305, 'Na+': 22.99, 'K+': 39.0983, 'Alt': 26.98}
    #for col in cation_mass.keys():
    #    data[col] = data[col] # * cation_mass[col] / 1e3 # 1e-6 mol/L => g/m^3
    return data_mean, data_std


def read_streamChem():
    """
    Stream chemistry (umol/L) measured on individual dates at monthly intervals
    """
    data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 
                                    'Hubbard_Brook', 'knb-lter-hbr.156.8', 'W1Long_StrmChem_HB1991-2020.csv'),
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
    data_std = data.groupby(data.index.get_level_values(0)).std()

    return data


def read_runoff():
    """ read streamflow reported in mm/day """
    streamflow = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 
                                          'Hubbard_Brook', 'knb-lter-hbr.2.14', 
                                          'HBEF_DailyStreamflow_1956-2023.csv'),
                             index_col = 0, parse_dates = True)
    streamflow = streamflow.loc[streamflow['WS'] == 1, 'Streamflow']
    return streamflow


def read_cec():
    data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 
                                    'Hubbard_Brook', 'knb-lter-hbr.158.1',
                                    'w1ffexchem.txt'), sep = ",")
    data['Horizon'] = data['Horizon'].map({'min': 'min', 'Oie': 'Oie', 'Oa': 'Oa', 
                                        'min': 'Min', 'cor': 'Min'})
    data = data.set_index(['Horizon', 'Year', 'Plot'])[['ExAcidcmolc_kg', 'ExCacmolc_kg', 'ExMgcmolc_kg', 'ExNacmolc_kg', 'ExKcmolc_kg', 'ExAlcmolc_kg']]
    data[data < 0] = np.nan
    data_mean = data.groupby(['Horizon', 'Year']).mean()
    data_std = data.groupby(['Horizon', 'Year']).std()
    return data_mean, data_std