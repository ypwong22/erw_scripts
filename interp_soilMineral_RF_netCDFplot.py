""" Interpolate to 0.5 degree and plot the topmost, middle, and bottom layer of interpolated for all the minerals. For visualization in the paper. """
import rasterio as rio
import os
import numpy as np
from scipy.interpolate import griddata
import itertools as it
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from string import ascii_lowercase


"""Perform the interpolation """
re_interp = True
if re_interp:
    for varname in ['Calcite','Kaolinit','Feldspar','Gibbsite']:

        h = rio.open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', f'RF_{varname}_filled.tif'))

        # tif file is in EPSG:5070
        left, bottom, right, top = h.bounds
        rows, cols = np.indices(h.shape)
        src_x, src_y = rio.transform.xy(h.transform, rows, cols)
        lons = np.array(src_x)
        lats = np.array(src_y)

        # conduct reprojection
        lats_target = np.arange(23.25, 54.26, 0.5)
        lons_target = np.arange(234.75, 293.26, 0.5) - 360.
        lons_mesh, lats_mesh = np.meshgrid(lons_target, lats_target)
        src_coords = np.array([lons.flatten(), lats.flatten()]).T

        # read data
        ds = np.full([10, len(lats_target), len(lons_target)], np.nan)
        for band in range(1, 11):
            data = h.read(band, masked = True)
            src_data_flat = data.flatten()

            data_reproj = griddata(src_coords, src_data_flat, (lons_mesh, lats_mesh),
                                method='linear')
            # use the nearest neighbor for the NE patch of missing data
            search_rows, search_cols = np.where(~np.isnan(data_reproj))
            for row, col in it.product(range(45, 49), range(110, 114)):
                # Northeastern corner
                if (row == 48) and (col == 110 or col == 111):
                    continue
                if (row == 47) and col == 110:
                    continue
                if np.isnan(data_reproj[row,col]):
                    ind = np.argmin(np.power(search_rows-row,2) + np.power(search_cols-col,2))
                    find_row, find_col = search_rows[ind], search_cols[ind]
                    data_reproj[row,col] = data_reproj[find_row, find_col]
            row = 4 # Florida 
            for col in [89, 90]:
                data_reproj[row, col] = data_reproj[row+1, col]
            col = 88
            for row in [4,5]:
                data_reproj[row, col] = data_reproj[6, col]

            ds[band-1, :, :] = data_reproj

        h.close()

        # Limit the values to >= 0
        

        # save to netcdf
        ds = xr.DataArray(ds, coords = {'lat': lats_target, 'lon': lons_target}, 
                            dims = ['nlevsoi','lat','lon'])
        filename = os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 
                                f'RF_{varname}_filled_interp.nc')
        if os.path.exists(filename):
            os.remove(filename)
        ds.to_dataset(name = f'PCT_{varname}').to_netcdf(filename)


""" Plot the interpolated result """
elm_nodes = np.array([0.0071, 0.0279, 0.0623, 0.1189, 0.2122, 0.3661, 0.6198, 1.0380, 1.7276, 2.8646])
selected_layers = np.array([0, 5, 9])
#elm_interface = np.array([0, 1.75, 4.51, 9.06, 16.55, 28.91, 49.29, 82.89, 138.28, 229.61, 380.19, 628.45])

fig, axes = plt.subplots(4, 3, figsize = (20, 18), subplot_kw = {'projection': ccrs.PlateCarree()})
fig.subplots_adjust(hspace = 0.15, wspace = 0.02)
for i, varname in enumerate(['Calcite','Kaolinit','Feldspar','Gibbsite']):
    for j, layer in enumerate(selected_layers):
        ax = axes[i,j]

        hr = xr.open_dataset(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results',
                                          f'RF_{varname}_filled_interp.nc'))
        levels = np.linspace(0, np.ceil(float(hr[f'PCT_{varname}'].max())), 21)
        cf = ax.contourf(hr['lon'], hr['lat'], hr[f'PCT_{varname}'][layer, :, :], cmap = 'Reds',
                         levels = levels)
        if i == 0:
            ax.set_title(f'{elm_nodes[layer]} m')
        if j == 0:
            ax.text(-0.05, 0.5, f'{varname} (g 100g-1 soil)', rotation = 90, transform = ax.transAxes, va = 'center')
        ax.text(0.05, 0.95, ascii_lowercase[i + j*4], transform = ax.transAxes, fontdict={'weight': 'bold'})
    cax = fig.add_axes([0.1, 0.7 - i*0.2, 0.8, 0.01])
    plt.colorbar(cf, cax = cax, orientation = 'horizontal')
plt.savefig(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', f'RF_mineral_all.png'), 
            dpi = 600., bbox_inches = 'tight')
plt.close(fig)
