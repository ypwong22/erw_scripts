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
    data = data.loc[data.index.year >= 2012, :]

    # Create a linear regression line and uncertainty intervals
    data_sp = data.loc[data['SEASON'] == 'SPRING', :]
    x = data_sp['DAY']
    y = data_sp['Phenology_Stage']
    res_sp = linregress(x, y)

    data_fa = data.loc[data['SEASON'] == 'FALL', :]
    x = data_fa['DAY']
    y = data_fa['Phenology_Stage']
    res_fa = linregress(x, y)

    return (res_sp.slope, res_sp.intercept, res_sp.stderr), (res_fa.slope, res_fa.intercept, res_fa.stderr)


def read_snowcourse():
    """ Snowcourse 2 is in Watershed 1 and has continued measurements till present """
    data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 
                                    'Hubbard_Brook', 'knb-lter-hbr.27.20',
                                    'HBEF_snowcourse_1956-2024.csv'),
                        index_col = [2, 1], parse_dates=True)
    data = data.loc[['STA2','STA9','STA17','STA19','STAHQ'], 'swe']
    data[data < 0] = np.nan
    data = data.loc[~data.index.duplicated(keep='first')]
    data = data.unstack().T['STA2'].sort_index()
    data = data.loc[data.index.year >= 2012]
    data = data.reindex(pd.date_range('2012-01-01','2024-04-10'))
    # linearly interpolate missing values
    data = data.interpolate(method = 'linear', limit = 14)
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
        'Oa': 'Oa', 'Oe': 'Oe', 'Oae': 'Oa', 
        'Bh': 'B', 'Bhs': 'B', 'Bhs1': 'B', 'Bhs3': 'B', 'Bs': 'B', 'Bs2': 'B'
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
    data = data[['Elevation','pH','Ca2', 'Mg2', 'Na', 'K', 'Alt']].sort_index()
    data[data < -800] = np.nan

    # merge sites 4a 4b 4c into a single site
    data = data.reset_index()
    data['Site'] = data['Site'].map({'1': '1', '2': '2', '2.5': '2.5', 
                                    '3': '3', '4a': '4', '4b': '4', '4c': '4'})

    # remove occasional entries that are measured at same site, slightly different elevations
    data = data.groupby(['Date','Site']).mean()

    data.drop('Elevation', axis = 1, inplace = True)

    # umol/L => mol/L
    data.columns = ['pH', 'Ca2+', 'Mg2+', 'Na+', 'K+', 'Al3+']
    data.iloc[:, 1:] = data.iloc[:, 1:].values * 1e-6

    # remove two outliers
    data.loc[data['Na+'] > 0.00015, 'Na+'] = np.nan
    data.loc[data['K+'] > 0.00005, 'K+'] = np.nan

    # average over all the sites in a watershed
    data_mean = data.groupby(data.index.get_level_values(0)).mean()
    data_std = data.groupby(data.index.get_level_values(0)).std()

    return data_mean, data_std


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
    data = data.set_index(['Horizon', 'Year', 'Plot'])[['ExAcidcmolc_kg', 'ExCacmolc_kg', 'ExMgcmolc_kg', 'ExNacmolc_kg', 'ExKcmolc_kg', 'ExAlcmolc_kg']].sort_index()
    data[data < 0] = np.nan
    # subtract the Al3+ from Acid exchange to get the H+
    data['ExAcidcmolc_kg'] = data['ExAcidcmolc_kg'] - data['ExAlcmolc_kg'].values
    data_mean = data.groupby(['Horizon', 'Year']).mean()
    data_std = data.groupby(['Horizon', 'Year']).std()
    return data_mean, data_std


def read_tsoi():
    data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 
                                    'Hubbard_Brook', 'knb-lter-hbr.315.2',
                                    'hbr315_snowgrad.soilT.csv'),
                       index_col = 0, parse_dates = True)
    data = data[['E10_soilT_1', 'E10_soilT_2']]

    # the second is the uncertainty between the two sensors
    return data


def read_solar_scan():
    """ Read the solar radiation data at the SCAN site """
    data_all = []
    for year in range(2012,2024):
        data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 
                                        'Hubbard_Brook', 'SCAN',
                                        f'2069_SOLAR_YEAR={year}.csv'),
                           skiprows = 6)
        data.index = pd.DatetimeIndex(data['Date'] + ' ' + data['Time'])
        data = data.loc[:,  'SRADV.H-1 (watt) ']
        data_all.append(data)
    data_all = pd.concat(data_all)
    return data_all

def read_solar():
    solrad = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 
                                      'Hubbard_Brook', 'knb-lter-hbr.237.1', 'HBEF_15minSolarRadiation_WX1.csv'))

    start_date = pd.Timestamp(year = 2014, month = 7, day = 21, hour = 10, minute = 30)
    end_date = pd.Timestamp(year = 2019, month = 4, day = 3, hour = 10, minute = 30)
    date_range = pd.date_range(start_date, end_date, freq='15T')

    solrad.index = date_range

    # average over two stations
    solrad = (solrad['SolRad1'] + solrad['SolRad2']) / 2
    solrad = solrad.resample('1H').mean()
    return solrad

def read_solar_daily():
    """ Solar radiation: Mj m-2 day-1 => W m-2 long-term daily """ 
    solrad_daily = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 
                                            'Hubbard_Brook', 'knb-lter-hbr.60.13', 
                                            'HBEF_DailySolarRadiation_HQ.csv'),
                               index_col = 0, parse_dates=True)['solar_rad']
    solrad_daily = solrad_daily.loc[solrad_daily.index.year >= 2012] * 277.78 / 24
    return solrad_daily

def read_era5(variable):
    filelist = [os.path.join(os.environ['PROJDIR'], "ERW_LDRD", "data", "GEE", "HubbardBrook", 
                             f"HubbardBrook_{variable}_{year}.csv") \
                for year in range(2011, 2024)]
    alldata = []
    for file in filelist:
        alldata = alldata + [pd.read_csv(file, index_col = 0, parse_dates=True)[variable]]
    alldata = pd.concat(alldata)
    if variable == "temperature_2m":
        alldata = alldata - 273.15
    if "radiation" in variable:
        alldata = alldata / 3600
    return alldata
