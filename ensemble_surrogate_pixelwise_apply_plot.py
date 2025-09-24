import pandas as pd
import os
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from matplotlib.colors import BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def get_path_root(which_pft) -> str:
    return os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble', which_pft)


def read_data(path_root, which_pft):
    """ read the CDR rate in gC/m2/year """
    if which_pft == 'pft1':
        cdr = pd.read_csv(os.path.join(path_root, 'combined_cdr.csv'), index_col = 0)
        is_predicted = np.full(cdr.shape[0], False)

        temp = (pd.read_csv(os.path.join(path_root, f'predictions_random_forest.csv'), 
                            index_col = 0) + \
                pd.read_csv(os.path.join(path_root, f'predictions_gradient_boosting.csv'),
                            index_col = 0)) * 0.5

        skip = np.array([2241, 2245, 2249, 2253, 2257, 2261, 2264])
        cdr.loc[skip, cdr.columns[4:]] = temp.loc[skip, :].values

        is_predicted[skip-1] = True
    else:
        cdr = pd.read_csv(os.path.join(path_root, 'predictions.csv'), index_col = 0)
        is_predicted = cdr['is_predicted'] # mask, True = surrogate model predicted
        cdr = cdr.drop(['is_predicted'], axis = 1)

    """ read the primary mineral application rate in g basalt/m2/year
        get from the benchmark PFT=1 since it is the same for all models
        verified identical for the simulated members for PFT=15 """
    pmr = pd.read_csv(os.path.join(path_root.replace(which_pft, 'pft1'), 'combined_pmr.csv'), index_col = 0)
    pmr = pmr.iloc[:, 4].values # this is identical to all columns

    """ calculate CDR rate in gC g basalt-1 """
    cdr_rate = cdr.copy()
    cdr_rate.iloc[:, 4:] = cdr.iloc[:, 4:] / pmr.reshape(-1, 1)

    ### remove negative values - they are not meaningful
    ##cdr_rate = cdr_rate.clip(lower = 0.)

    return cdr_rate, is_predicted


def plot_CONUS(cdr_rate, path_root):
    """ plot the CONUS-average CDR rate for the ensemble members """
    fig, ax = plt.subplots(figsize = (20, 10))
    cdr_rate.iloc[:, 4:].mean(axis = 1).plot(ax = ax)
    ax.set_title('CONUS average CDR rate for each ensemble member (gC g basalt-1)')
    ax.set_ylabel('')
    fig.savefig(os.path.join(path_root, 'ensemble_CDR_rate.png'), dpi = 600., bbox_inches = 'tight')
    plt.close(fig)


def df_to_map(df, which_pft):
    """ convert the pixel IDs to spatial points """

    if 'SSA (m2 g-1)' in df.columns or 'Grain size (um)' in df.columns:
        df = df.iloc[:, 4:] # remove the label columns to keep just data

    # Obtain spatial mask
    nc = Dataset(os.path.join(os.environ['ZDR'], 'ERW', 'output', 'UQ', which_pft,
                              'r_sequestration_diff.nc'))
    mask_data = nc['r_sequestration_diff'][:,:,:].mask  # where no data exist
    ensemble = nc['ensemble'][:]
    lat = nc['lat'][:]
    lon = nc['lon'][:]
    nc.close()

    n_ens, n_lat, n_lon = mask_data.shape
    # invalid spatial pts must be where no ensemble has data
    mask2d = np.all(mask_data, axis = 0)

    # Decide whether df columns correspond to True or False in the mask
    n_cols = df.shape[1]
    valid_pts = np.flatnonzero(~mask2d) # verified to be the same as obtained in "ensemble_to_CSV.py"
    if n_cols != valid_pts.size:
        raise ValueError(f"df has {n_cols} columns, but mask has {valid_pts.size} valid spatial pixels; cannot align.")

    # Reindex df rows to match NetCDF ensemble order
    try:
        df_aligned = df.reindex(ensemble)
    except Exception:
        # If index types don't match, coerce to the same dtype
        df_aligned = df.copy()
        df_aligned.index = df_aligned.index.astype(ensemble.dtype, copy=False)
        df_aligned = df_aligned.reindex(ensemble)

    new_data = np.full((n_ens, n_lat, n_lon), np.nan, dtype=float)

    for i in range(n_ens):
        row_vals = df_aligned.iloc[i].to_numpy().astype(float, copy=False)
        if row_vals.size != valid_pts.size:
            raise ValueError(f"Row {i} length {row_vals.size} does not match target positions {valid_pts.size}.")
        flat = new_data[i].ravel()
        flat[valid_pts] = row_vals
        new_data[i] = flat.reshape(n_lat, n_lon)

    return new_data, ensemble, lat, lon


def save_new_data_to_netcdf(new_data, ensemble, lat, lon, out_path, var_name='cdr_rate'):
    """Save a [ensemble, lat, lon] array to NetCDF with coordinates."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    ens_dtype = 'i4' if np.issubdtype(ensemble.dtype, np.integer) else 'f8'

    with Dataset(os.path.join(out_path, 'cdr_rate.nc'), 'w', format='NETCDF4') as ds:
        # Dimensions
        ds.createDimension('ensemble', len(ensemble))
        ds.createDimension('lat', len(lat))
        ds.createDimension('lon', len(lon))

        # Coordinates
        vens = ds.createVariable('ensemble', ens_dtype, ('ensemble',))
        vlat = ds.createVariable('lat', 'f8', ('lat',))
        vlon = ds.createVariable('lon', 'f8', ('lon',))
        vens[:] = ensemble
        vlat[:] = lat
        vlon[:] = lon
        vlat.standard_name = 'latitude'
        vlat.units = 'degrees_north'
        vlon.standard_name = 'longitude'
        vlon.units = 'degrees_east'

        # Data variable
        v = ds.createVariable(var_name, 'f4', ('ensemble', 'lat', 'lon'),
                              zlib=True, complevel=4, fill_value=np.nan)
        v[:] = new_data.astype('f4', copy=False)

        # Minimal attrs
        ds.title = f'{var_name} on [ensemble, lat, lon]'
        ds.history = 'Created by save_new_data_to_netcdf'


def map_new_data(out_path):
    """ Create maps of the max & min CDR rates over conus """
    with Dataset(os.path.join(out_path, 'cdr_rate.nc'), 'r', format='NETCDF4') as ds:
        cdr_rate = ds['cdr_rate'][:, :, :]
        lon = ds['lon'][:]
        lat = ds['lat'][:]

        levels = np.linspace(-0.16, 0.16, 41)
        norm = BoundaryNorm(boundaries=levels, ncolors=256, clip=True)

        fig, axes = plt.subplots(2, 1, figsize=(20, 17), subplot_kw={'projection': ccrs.PlateCarree()})

        for ax in axes:
            ax.add_feature(
                cfeature.STATES.with_scale('110m'),
                edgecolor='black',
                linewidth=0.6
            )

        ax = axes.flat[0]

        spatial_max = cdr_rate.max(axis=0)
        cf = ax.pcolormesh(lon, lat, spatial_max, cmap='Spectral', norm=norm)
        ax.set_title('Maximum CDR across ensemble members (gC g basalt$^{{-1}}$)')
        plt.colorbar(cf, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.7)

        # Find where the largest sequstration rate occurs
        flat  = np.argmax(spatial_max)
        row, col = divmod(flat, spatial_max.shape[1])
        lat_max = lat[row]
        lon_max = lon[col]

        ax.scatter(lon_max, lat_max, color='r', s=50, marker='*', fc = 'none', transform=ccrs.PlateCarree())
        ax.text(0.05, 0.95, f'Maximum at {lat_max:.2f}°N, {lon_max:.2f}°W',
                color='black', fontsize=12, transform=ax.transAxes)

        ax = axes.flat[1]
        spatial_min = cdr_rate.min(axis=0)
        cf = ax.pcolormesh(lon, lat, spatial_min, cmap='Spectral', norm=norm)
        ax.set_title('Minimum CDR across ensemble members (gC g basalt$^{{-1}}$)')
        plt.colorbar(cf, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.7)

        plt.savefig(os.path.join(out_path, 'map_CDR_rate.png'), dpi = 600., bbox_inches = 'tight')
        plt.close(fig)


if __name__ == '__main__':
    which_pft = 'pft1' # 'pft15'
    path_root = get_path_root(which_pft)
    cdr_rate, is_predicted = read_data(path_root, which_pft)
    plot_CONUS(cdr_rate, path_root)
    new_data, ensemble, lat, lon = df_to_map(cdr_rate, which_pft)
    save_new_data_to_netcdf(new_data, ensemble, lat, lon, path_root, var_name='cdr_rate')
    map_new_data(path_root)