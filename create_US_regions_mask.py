import numpy as np
import geopandas as gpd
from shapely.geometry import Point
import xarray as xr
import datetime as dt

# ------------------------------------------------------------------
# 0.  Lat/Lon arrays (exactly what you supplied)
# ------------------------------------------------------------------
lat_centres = np.arange(25.25, 48.25 + 0.01, 0.5)        # 47 values
lon_centres = np.arange(235.75, 292.25 + 0.01, 0.5)      # 114 values
lon_centres = (lon_centres + 180) % 360 - 180            # wrap to –180…180

lon2d, lat2d = np.meshgrid(lon_centres, lat_centres)     # shape (47, 114)

# ------------------------------------------------------------------
# 1.  Build a GeoDataFrame of every grid-point
# ------------------------------------------------------------------
points = gpd.GeoDataFrame(
    geometry=[Point(xy) for xy in zip(lon2d.ravel(), lat2d.ravel())],
    crs="EPSG:4326"
)

# ------------------------------------------------------------------
# 2.  Download US state polygons (48-state CONUS only)
# ------------------------------------------------------------------
url = "https://www2.census.gov/geo/tiger/GENZ2022/shp/cb_2022_us_state_500k.zip"
states = (
    gpd.read_file(url)
      .to_crs("EPSG:4326")
      .loc[~lambda df: df["STUSPS"].isin(["AK", "HI", "PR", "GU", "MP", "VI", "AS"])]
)

# ------------------------------------------------------------------
# 3.  Attach a “region” field to each state
#     – tweak the lists if you prefer another convention –
# ------------------------------------------------------------------
region_map = {
    "Northeast":  ["ME","NH","VT","MA","RI","CT","NY","NJ","PA"],
    "Midwest":    ["OH","MI","IN","IL","WI","MN","IA","MO","ND","SD","NE","KS"],
    "Southwest":  ["AZ","NM","TX","OK"],
    "West":       ["CA","OR","WA","NV","ID","UT","CO","MT","WY"],
    "Southeast":  ["DE","MD","DC","VA","WV","KY","TN",
                   "NC","SC","GA","FL","AL","MS","AR","LA"],
}
states["region"] = None
for reg, st_list in region_map.items():
    states.loc[states["STUSPS"].isin(st_list), "region"] = reg

# ------------------------------------------------------------------
# 4.  Spatial join: point → state → region
# ------------------------------------------------------------------
points = gpd.sjoin(
    points,                         # left = every grid cell
    states[["region", "geometry"]], # right = state polygons w/ region
    how="left",
    predicate="within"
)

# ------------------------------------------------------------------
# 5.  Build boolean masks (47 × 114) for each region
# ------------------------------------------------------------------
grid_shape = lat2d.shape           # (47, 114)
masks = {
    reg: (points["region"] == reg).values.reshape(grid_shape)
    for reg in region_map
}

# ------------------------------------------------------------------
# 6.  Quick sanity-check: number of land cells per region
# ------------------------------------------------------------------
for reg, m in masks.items():
    print(f"{reg:10s}  {m.sum():4d} cells")

# ------------------------------------------------------------------
# 6.  Package the five masks in an xarray Dataset and write NetCDF
# ------------------------------------------------------------------

# masks  ➜  3-D boolean array (region, lat, lon)
regions   = list(region_map.keys())                 # ['Northeast', …, 'West']
stack_arr = np.stack([masks[r] for r in regions])   # shape (5, 47, 114)

ds = xr.Dataset(
    data_vars={
        reg: (("lat", "lon"), masks[reg].astype("int8"))
        for reg in regions
    },
    coords=dict(
        lat=lat_centres,
        lon=lon_centres,
    ),
    attrs=dict(
        title       ="U.S. region masks on 0.5° grid",
        description ="Boolean (0/1) masks derived from 2022 U.S. Census TIGER/Line state polygons. "
                     "Regions follow NCEI/NOAA convention (West, Southwest, Midwest, Southeast, Northeast).",
        grid        ="Centers at 25.25°N–48.25°N, 235.75°E–292.25°E (-124.25°W to -67.75°W)",
        created_on  =dt.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC"),
        software    ="Python, GeoPandas, Shapely, xarray, NumPy",
    )
)

# Write the file
out_path = "us_region_masks_0.5deg.nc"
ds.to_netcdf(out_path, encoding={reg: {"_FillValue": 0} for reg in regions})

print(f"✓ NetCDF written to {out_path}")
