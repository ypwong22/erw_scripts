import os
import xarray as xr
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed


def read_obs_UIEF():
    """ Read observations """
    path_data = os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 'UIEF')

    # (4) Weathering difference
    weathering_difference = pd.read_excel(os.path.join(path_data, 'pnas.2319436121.sd06.xlsx'),
                                        header = 2, usecols = [0,49,50,52,53],
                                        index_col = 0, 
                                        skiprows = list(range(2)) + list(range(4,34)) + \
                                                   list(range(39,80)))
    weathering_difference.columns = ['Ca weathered', 'stdev.1', 'Mg weathered', 'stdev.2']

    # (5) CEC %, extracted from the plot
    cec_pct = {
        ('base sat', '0-10cm' , 'treatment'): [80.35714285714286, 82.94642857142857, 
                                                89.28571428571428, 87.58928571428572, 
                                                82.94642857142857],
        ('base sat', '0-10cm' , 'control'  ): [78.92857142857143, 80.71428571428572, 85, 
                                                78.92857142857143, 72.67857142857143],
        ('base sat', '10-30cm', 'treatment'): [83.125, 79.55357142857143, 80.89285714285714,
                                                80.44642857142857, 86.51785714285714],
        ('base sat', '10-30cm', 'control'  ): [81.07142857142857, 77.14285714285714, 
                                                76.51785714285714, 75.80357142857143, 
                                                82.76785714285714],
        ('acid'    , '0-10cm' , 'control'  ): [20.416666666666668, 19.09722222222222, 
                                            14.722222222222221, 20.694444444444443, 
                                            26.875000000000004],
        ('acid'    , '0-10cm' , 'treatment'): [18.88888888888889, 16.73611111111111,
                                            10.416666666666664, 12.013888888888888, 
                                            16.73611111111111],
        ('acid'    , '10-30cm', 'control'  ): [18.819444444444443, 22.63888888888889,
                                            23.055555555555554, 23.61111111111111, 
                                            16.80555555555555],
        ('acid'    , '10-30cm', 'treatment'): [16.666666666666668, 20.20833333333333, 
                                            18.88888888888889, 19.23611111111111, 13.26388888888889]
    }
    cec_pct = pd.DataFrame(cec_pct, index = range(2016, 2021))

    # (6) Soil pH, extracted from the plot
    obs_soil_pH = {
        ('0-10cm', 'treatment' ): [6.131292517006803, 6.6020408163265305, 
                                6.661904761904761, 6.92312925170068, 6.925850340136054],
        ('0-10cm', 'control'   ): [6.210204081632653, 6.280952380952381,
                                6.313605442176871, 6.134013605442177, 5.97891156462585],
        ('10-30cm', 'treatment'): [6.109523809523809, 6.406122448979592,
                                6.5639455782312925, 6.485034013605442, 6.512244897959183],
        ('10-30cm', 'control'  ): [6.022448979591837, 6.289115646258503,
                                6.316326530612245, 6.131292517006803, 6.245578231292517]
    }
    obs_soil_pH = pd.DataFrame(obs_soil_pH, index = range(2016, 2021))

    # subset to treatment only, because ensemble simulation is treatment
    obs_soil_pH = obs_soil_pH.reorder_levels([1,0], axis = 1)['treatment']
    cec_pct = cec_pct.reorder_levels([2,0,1], axis = 1)['treatment']

    return obs_soil_pH, cec_pct, weathering_difference


def read_sim_UIEF(ensemble_id):
    """ Function to obtain UIEF data """

    # Global constants
    nlevsoi = 10
    zsoi = np.array([0., 0.01751282, 0.04509179, 0.09056182, 0.16552924, 0.2891296,
                        0.49291214, 0.82889277, 1.3828311, 2.2961211, 3.8018818])
    thickness = np.diff(zsoi)

    depth_ranges = [(0., 0.1), (0.1, 0.3)]
    depth_names = ['0-10cm', '10-30cm']
    tvec = np.arange(2015, 2022)

    # Soil bulk density at the site, obtained from other runs
    bd_col = np.array([1440.936, 1440.936, 1440.936, 1434.132, 1434.132, 1420.524,
                       1437.534, 1437.534, 1437.534, 1437.534, 1437.534, 1437.534,
                       1437.534, 1437.534, 1437.534])

    case_name = '20251002_UIEF_ICB20TRCNPRDCTCBC_erw'
    flist = [os.path.join(os.environ['E3SM_ROOT'], 'output', 'UQ', case_name, f'g{ensemble_id:05d}', 
                          f'{case_name}.elm.h1.{yy}-01-01-00000.nc') \
             for yy in range(2015, 2022)]

    ## debug test
    ##case_name = '20250920_UIEF_ICB20TRCNPRDCTCBC_obsForc'
    ##flist = [os.path.join(os.environ['E3SM_ROOT'], 'output', f'{case_name}erw', 'run', 
    ##                        f'{case_name}erw.elm.h1.{yy}-01-01-00000.nc') \
    ##         for yy in range(2015, 2022)]

    hr = xr.open_mfdataset(flist)

    # Two depths
    soil_pH = pd.DataFrame(np.nan, index = tvec, columns = depth_names)
    for i, (min_depth, max_depth) in enumerate(depth_ranges):
        dname = depth_names[i]
        weight = np.where((zsoi[:-1] < max_depth) & (zsoi[1:] > min_depth),
                          np.minimum(thickness, max_depth - zsoi[:-1], zsoi[1:] - min_depth),
                          0)
        soil_pH[dname] = np.nansum(hr['soil_pH'].resample(time = '1Y').mean().values[:,:nlevsoi,0] * \
                                   weight.reshape(1, -1), axis = 1) / np.sum(weight)

    # CEC: convert from g m-3 soil to meq 100g-1 dry soil
    pct_acid = pd.DataFrame(np.nan, index = tvec, columns = depth_names)
    base_saturation = pd.DataFrame(np.nan, index = tvec, columns = depth_names)
    for i, (min_depth, max_depth) in enumerate(depth_ranges):
        weight = np.where((zsoi[:-1] < max_depth) & (zsoi[1:] > min_depth),
                          np.minimum(thickness, max_depth - zsoi[:-1], zsoi[1:] - min_depth),
                          0)
        dname = depth_names[i]

        soil_cec_sim = {}
        for i, (cation, val_icat, mass_icat) in \
            enumerate(zip(['Ca2+','Mg2+','Na+','K+','Al3+'],
                          [2, 2, 1, 1, 3], [40, 24, 23, 29, 27])):
            soil_cec_sim[cation] = \
                ((hr[f'cec_cation_vr_{1+i}'][:,:nlevsoi,0].resample(time='1Y').mean() * \
                  weight.reshape(1, -1)).sum(axis = 1) * 1000 / 10 * val_icat / mass_icat / \
                 (bd_col[:nlevsoi] * weight).sum()).values

        soil_cec_sim['H+'] = \
            ((hr[f'cec_proton_vr'][:,:nlevsoi,0].resample(time='1Y').mean() * \
              weight.reshape(1, -1)).sum(axis = 1) * 1000 / 10 / \
             (bd_col[:nlevsoi] * weight).sum()).values

        soil_cec_sim = pd.DataFrame(soil_cec_sim, index = tvec)

        pct_acid[dname] = 100 * (soil_cec_sim[['H+','Al3+']].sum(axis = 1) / soil_cec_sim.sum(axis = 1))
        base_saturation[dname] = 100 * (soil_cec_sim.drop(['H+','Al3+'], axis = 1).sum(axis = 1) / soil_cec_sim.sum(axis = 1))

    # Dissolved amount of basalt, converted to amount of lost Ca and Mg
    # g m-2
    # note: need to divide by 89% to get total basalt amount
    total_residing = hr['primary_mineral'][:, :, 0].sum(axis = 1)
    total_residing = total_residing.resample(time = '1Y').last()

    # g m-2 s-1 => cumulative g m-2
    # note: need to divide by 89% to get total basalt amount
    cum_applied = (hr['primary_added'][:, :, 0].sum(axis = 1) * 3600 * 24).cumsum()
    cum_applied = cum_applied.resample(time = '1Y').last()

    # dissolved basalt, assume 58% in the 0-10cm like in the paper
    total_dissolved = (cum_applied - total_residing) * 0.58

    # convert from g m-2 to ug g-1 soil in top 10cm, using bd = 1.2 t m-3
    total_dissolved = total_dissolved * 1e6 / (1.2e6 * 0.1)

    # convert to Ca & Mg, in ug g-1
    mineral_mass = 262.22 * 0.233 + 483.22 * 0.178 + 100.0872 * 0.026 + 812.3665 * 0.116 + 555.7973 * 0.34
    ca_mass = 40 * (2*0.178 + 0.026 + 2*0.116)
    mg_mass = 24 * (5*0.116 + 5*0.34)
    total_dissolved_ca = pd.Series(total_dissolved.values * ca_mass / mineral_mass, index = tvec)
    total_dissolved_mg = pd.Series(total_dissolved * mg_mass / mineral_mass, index = tvec)

    hr.close()

    return soil_pH, pct_acid, base_saturation, total_dissolved_ca, total_dissolved_mg


def read_param_UIEF(ensemble_id):
    """ Function to get the """
    case_name = '20251002_UIEF_ICB20TRCNPRDCTCBC_erw'
    hr = xr.open_dataset(os.path.join(os.environ['E3SM_ROOT'], 'output', 'UQ', case_name,
                                      f'g{ensemble_id:05d}', f'surfdata_{ensemble_id:05d}.nc'))
    params = {}
    for i in range(1,6):
        var = f'LOG_KM_{i}'
        params[('param',var)] = float(hr[var][0,0,0])
    hr.close()

    return pd.Series(params)


def calc_kge(ensemble_id):
    """ Use Kling-Gupta Efficiency to simultaneously obtain correlation, mean bias, and variability """

    def _kge(obs, sim):
        obs, sim = np.array(obs), np.array(sim)

        # correlation
        r = np.corrcoef(obs, sim)[0, 1]

        # mean bias ratio
        beta = np.mean(sim) / np.mean(obs)

        # variability ratio (coefficient of variation)
        gamma = (np.std(sim) / np.mean(sim)) / (np.std(obs) / np.mean(obs))

        # Kling-Gupta Efficiency
        return 1 - np.sqrt((r - 1)**2 + (beta - 1)**2 + (gamma - 1)**2)

    obs_soil_pH, cec_pct, weathering_difference = read_obs_UIEF()
    soil_pH, pct_acid, base_saturation, total_dissolved_ca, total_dissolved_mg = read_sim_UIEF(ensemble_id)

    rmse = pd.DataFrame(np.nan, index = ['soil_pH', 'pct_acid', 'base_saturation', 'total_dissolved_ca', 'total_dissolved_mg'], 
                        columns = ['0-10cm', '10-30cm'])
    for depth in rmse.columns:
        rmse.loc['soil_pH', depth] = _kge(obs_soil_pH[depth], soil_pH.loc[obs_soil_pH.index,depth])
        rmse.loc['pct_acid', depth] = _kge(cec_pct[('acid',depth)], pct_acid.loc[cec_pct.index,depth])
        rmse.loc['base_saturation', depth] = _kge(cec_pct[('base sat',depth)], base_saturation.loc[cec_pct.index,depth])
        rmse.loc['total_dissolved_ca', depth] = _kge(weathering_difference['Ca weathered'], total_dissolved_ca.loc[weathering_difference.index]) # No depth in this one
        rmse.loc['total_dissolved_mg', depth] = _kge(weathering_difference['Mg weathered'], total_dissolved_mg.loc[weathering_difference.index])
    rmse = rmse.unstack()

    params = read_param_UIEF(ensemble_id)

    return pd.concat([params, rmse])


def parallel_process(ens_list, max_workers=8):
    results = {}

    # Run in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_folder = {executor.submit(calc_kge, f): f for f in ens_list}

        for future in as_completed(future_to_folder):
            folder_id = future_to_folder[future]
            try:
                series_result = future.result()
                results[folder_id] = series_result
            except Exception as e:
                print(f"Error processing folder {folder_id}: {e}")

    # Combine results into a DataFrame
    df = pd.DataFrame.from_dict(results, orient="index")
    return df


df = parallel_process(range(1,1001), max_workers=25)
df.to_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble_UIEF_kge.csv'))
