import numpy as np
import rasterio as rio
from pyproj import Transformer, CRS
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.windows import from_bounds
from scipy.stats import gaussian_kde
import os


def vert_interp(target_nodes, input_nodes, input_data, target_single_level, 
    target_interfaces, input_interfaces):
    """
    Linearly interpolate soil moisture/soil temperature from input_nodes to target_nodes. 
    If the target depths are single-level, returns weighted average based on the distance 
        between the input nodes and target node. 
    If the target depths are defined by bounds (target_interfaces != None), returns weighted
        average based on the overlapping lengths between the target_interfaces and 
        input_interface. 

    Parameters:
    -----------
    target_nodes : Union[List[float], np.ndarray]
        List or numpy array of target node depths in meters.
    
    input_nodes : Union[List[float], np.ndarray]
        List or numpy array of input node depths in meters.

    input_data : np.ndarray
        2D numpy array of input data with shape (time, len(input_nodes)).

    target_single_level : bool
        Indicates whether the target nodes are single level. If true, target_interface and input_interfaces are un-used. 

    Returns:
    --------
    np.ndarray
        Processed data as a 2D numpy array with shape (time, len(target_nodes)).
    """
    # unifying data types
    target_nodes = np.array(target_nodes)
    input_nodes = np.array(input_nodes)

    # sanity checks
    if not input_data.shape[1] == len(input_nodes):
        raise Exception('Mismatch between specified inputs depths and available data')
    if not target_single_level:
        if (target_interfaces is None or input_interfaces is None):
            raise Exception('Must specify depth bounds if target is not single level')
        if not ((len(target_interfaces) - len(target_nodes)) == 1):
            raise Exception('Number of soil layers mismatched between target interface and node depths')
        if not ((len(input_interfaces) - len(input_nodes)) == 1):
            raise Exception('Number of soil layers mismatched between input interface and node depths')
        # unifying data types
        target_interfaces = np.array(target_interfaces)
        input_interfaces = np.array(input_interfaces)

    # actual calculations
    output_data = np.full([input_data.shape[0], len(target_nodes)], np.nan)
    if target_single_level:
        for i, d in enumerate(target_nodes):
            if d < input_nodes[0]:
                output_data[:, i] = input_data[:, 0]
            elif d > input_nodes[-1]:
                output_data[:, i] = input_data[:, -1]
            else:
                d_matched = np.where(np.isclose(input_nodes, d))[0]
                if len(d_matched) > 1:
                    raise Exception('Input nodes have duplicate values')
                elif len(d_matched) == 1:
                    # just apply the nearest node
                    output_data[:, i] = input_data[:, d_matched[0]]
                else:
                    # interpolate between two nearby nodes
                    d_up = np.where(input_nodes < d)[0][-1]
                    d_down = np.where(input_nodes > d)[0][0]
                    f1 = (input_nodes[d_down] - d) / (input_nodes[d_down] - input_nodes[d_up]) 
                    f2 = (d - input_nodes[d_up]) / (input_nodes[d_down] - input_nodes[d_up])
                    output_data[:, i] = input_data[:, d_up] * f1 + input_data[:, d_down] * f2
    else:
        for i, d1 in enumerate(target_interfaces[:-1]):
            d2 = target_interfaces[i+1]
            if d2 <= input_interfaces[1]:
                output_data[:, i] = input_data[:, 0]
            elif d1 >= input_interfaces[-2]:
                output_data[:, i] = input_data[:, -1]
            else:
                output_data[:, i] = 0.
                sum_weight = 0.
                for j, dd1 in enumerate(input_interfaces[:-1]):
                    dd2 = input_interfaces[j+1]
                    if (dd2 <= d1) or (dd1 >= d2):
                        continue
                    else:
                        if (dd1 >= d1):
                            if (dd2 <= d2):
                                sum_weight += (dd2 - dd1)
                                output_data[:, i] += input_data[:, j] * (dd2 - dd1)
                            else:
                                sum_weight += (d2 - dd1)
                                output_data[:, i] += input_data[:, j] * (d2 - dd1)
                        else:
                            if (dd2 <= d2):
                                sum_weight += (dd2 - d1)
                                output_data[:, i] += input_data[:, j] * (dd2 - d1)
                            else:
                                sum_weight = 1.
                                output_data[:, i] = input_data[:, j]
                                break
                output_data[:, i] /= sum_weight
                # print(i, d1, d2, sum_weight)

    return output_data


def read_raster(src_file, lat_list, lon_list, band):
    """ Extract from geotiff file values at given lat lon locations """
    h = rio.open(src_file)
    array = h.read(band, masked = True)
    array = np.where(array.mask, np.nan, array)
    if h.crs != CRS.from_epsg(4326):
        transformer = Transformer.from_crs(CRS.from_epsg(4326), h.crs, always_xy=True)
    # List to store the extracted values
    extracted_values = []
    for lat, lon in zip(lat_list, lon_list):
        if h.crs != CRS.from_epsg(4326):
            # Transform lat/lon to x/y coordinates in the GeoTIFF's CRS
            x, y = transformer.transform(lon, lat)
        else:
            x, y = lon, lat
        # Get the pixel coordinates
        row, col = h.index(x, y)
        # Check if the point is within the raster bounds
        if 0 <= row < h.height and 0 <= col < h.width:
            # Extract values for all bands
            values = array[row, col]
            extracted_values.append(values)
        else:
            extracted_values.append(np.nan)  # Out of bounds
    h.close()
    return np.array(extracted_values)


def reproject_raster(src_path, dst_path, dst_crs='EPSG:4326', resolution=0.5, 
                     resampling=Resampling.bilinear):
    """ Resample a geotiff file to 0.5 degrees lat lon """
    with rio.open(src_path) as src:
        # Calculate the transform for the new CRS and resolution
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds, resolution = resolution
        )
        # Update the metadata for the output file
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })
        # Create the reprojected raster file
        with rio.open(dst_path, 'w', **kwargs) as dst:
            # Reproject and save each band
            for i in range(1, src.count + 1):
                reproject(
                    source=rio.band(src, i),
                    destination=rio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=resampling
                )
    print(f"Reprojected raster saved to {dst_path}")


def local_density(x, y):
    """ For scatter plots, return the vector that can be used
        to color the dots by their local density
    """
    # Calculate the point density
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    return x, y, z


def get_reproj_bounds():
    """ Determine the overlapping boundary of the reprojected file series
        for predictors use in "interp_NCSS_RF.ipynb" and "interp_soilMineral_RF.ipynb"
    """
    # Determine the overlapping boundary of the reprojected file series
    h1 = rio.open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 'SSURGO', 
                            'reproj_caco3_kg_sq_m.tif'))
    h2 = rio.open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 
                            'reproj_USGS_3DEP_800m_EPSG5070.tif'))
    h3 = rio.open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 
                            'reproj_DAYMET_V4_Climatology_1980_2010.tif'))
    h4 = rio.open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 'NLCD',
                            'reproj_NLCD_2001.tif'))

    bounds1 = h1.bounds
    bounds2 = h2.bounds
    bounds3 = h3.bounds
    bounds4 = h4.bounds

    # Calculate the overlapping bounds
    overlap_bounds = (
        max(bounds1.left, bounds2.left, bounds3.left, bounds4.left),
        max(bounds1.bottom, bounds2.bottom, bounds3.bottom, bounds4.bottom),
        min(bounds1.right, bounds2.right, bounds3.right, bounds4.right),
        min(bounds1.top, bounds2.top, bounds3.top, bounds4.top)
    )

    h1.close()
    h2.close()
    h3.close()
    h4.close()

    return overlap_bounds


def read_limited_raster(file_path, band, bounds):
    h = rio.open(file_path)
    window = from_bounds(*bounds, h.transform)
    data = h.read(band, window = window, masked=True)
    data = np.where(data.mask, np.nan, data.data)
    h.close()
    return data