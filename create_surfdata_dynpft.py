import os
import xarray as xr
import numpy as np

path_root = os.path.join(os.environ['PROJDIR'], 'E3SM', 'inputdata', 'lnd', 'clm2', 'PTCLM')

mpft = 17
namendspec = 10

for site in ['UC_Davis', 'HBR']:

    path_surffdata = os.path.join(path_root, site, 'surfdata.pftdyn.nc')

    os.system(f'cp {path_surffdata} {path_surffdata}_temp')
    hr = xr.open_dataset(f'{path_surffdata}_temp')

    # minerals_name = ['Wollastonite_CaSiO3', 'Forsterite_Mg2SiO4', 'Albite_NaAlSi3O8', 
    #                 'Anorthite_CaAl2Si2O8', 'Epidote_Ca2FeAl2(SiO4)3(OH)', 'Calcite_CaCO3',
    #                 'Labradorite_Ca0.6Na0.4Al1.6Si2.4O8', 
    #                 'Augite_Ca0.9Mg0.9Na0.1Al0.4Fe0.2Si1.9O6',
    #                 'Kfeldspar_KAlSi3O8', 'Enstatite_MgSiO3']
    if site == 'UC_Davis':
        doy = 274. # Oct-1
        rate = 4. # 40 t ha-1 = 4 kg / m2
        size = 105. # um
        pct = np.array([0, 0, 33.4, 33.4, 14.3, 0, 0, 0, 0, 0])
    elif site == 'HBR':
        doy = 292. # Oct-19
        rate = 0.466 # 55 tons / 11.8 ha = 0.466 kg / m2
        size = 9.6 # um
        pct = np.array([100., 0, 0, 0, 0, 0, 0, 0, 0, 0])
    year_of_application = np.where(hr['time'].values == 2000)[0][0]

    temp = np.full([len(hr['time']), mpft, 1,1], np.nan)
    temp[year_of_application, :, :, :] = doy
    hr['SOIL_AMENDMENTS_DOY'] = xr.DataArray(temp,
        dims = ['time','mpft','lsmlat','lsmlon'],
        attrs = {'long_name': 'soil amendment application day of year',
                 'units': 'day of year'})

    temp = np.full([len(hr['time']), mpft, 1,1], np.nan)
    temp[year_of_application, :, :, :] = rate
    hr['SOIL_AMENDMENTS_RATE'] = xr.DataArray(temp,
        dims = ['time','mpft','lsmlat','lsmlon'],
        attrs = {'long_name': 'soil amendment application rate',
                 'units': 'kg/m2/yr'})

    temp = np.full([len(hr['time']), 1,1], np.nan)
    temp[year_of_application, :, :] = size
    hr['SOIL_AMENDMENTS_GRAINSIZE'] = xr.DataArray(temp, 
        dims = ['time','lsmlat','lsmlon'],
        attrs = {'long_name': 'grain size of applied soil amendment',
                 'units': 'micro meters'})

    temp = np.full([len(hr['time']), namendspec, 1,1], np.nan)
    temp[year_of_application, :, 0, 0] = pct
    hr['SOIL_AMENDMENTS_PCT'] = xr.DataArray(temp, 
        dims = ['time','namendspec','lsmlat','lsmlon'],
        attrs = {'long_name': 'species fraction of applied soil amendment',
                 'units': 'percent'})

    hr['namendspec'] = xr.DataArray(
          np.arange(1, namendspec+1), dims = ['namendspec'],
          attrs = {'long_name': 'indices of species in soil amendment mixture (e.g. rock powder)',
                   'units': 'index'}
    )

    hr['mpft'] = xr.DataArray(
          np.arange(1, mpft+1), dims = ['mpft'],
          attrs = {'long_name': 'indices of natural PFTs and cfts, if any',
                   'units': 'index'}
    )

    hr.to_netcdf(path_surffdata.replace('.nc', '_erw.nc'))
    hr.close()

    os.system(f'rm {path_surffdata}_temp')