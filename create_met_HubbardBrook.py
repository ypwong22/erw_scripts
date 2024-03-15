import numpy as np
from netCDF4 import Dataset
import os, sys
import write_elm_met
import pandas as pd


# ------- user input -------------
site = "HBR_1"
start_year = 2003
end_year = 2023
time_offset = 0  # Standard time offset from UTC - but the CSV file is already aligned
npd = 24  # number of time steps per day (48 = half hourly)
if site == 'HBR_1':
    mylon = 360-71.728825  # site longitude (0 to 360)
    mylat = 43.9556695  # site latitude
elif site == 'HBR_2':
    mylon = 360-71.478825
    mylat = 43.9556695
measurement_height = 5  # tower height (m)
filename = os.path.join('./temp', f'downscaling_result_test.csv')
calc_flds = False  # use T and RH to comput FLDS (use if data missing or sparse)
leapdays = True  # input data has leap days (to be removed for ELM)
outdir = os.path.join('./temp', f'1x1pt_{site}')  # Desired directory for ELM met inputs

# outvars   - met variables used as ELM inputs
# invars    - corresponding variables to be read from input file
# conv_add  - offset for converting units (e.g. C to K)
# conv_mult - multiplier for converting units (e.g. hPa to Pa, PAR to FSDS)
# valid_min - minimum acceptable value for this variable (set as NaN outside range)
# valid_max - maximum acceptable value for this variable (set as NaN outside range)

outvars = ["TBOT", "RH", "WIND", "PSRF", "FSDS", "FLDS", "PRECTmms"]
invars = ["Tair", "RH (%)", "Wind", "Pressure (Pa)", "RSDS", "RLDS (W/m2)", "Precip (mm/hour)"]

conv_add = [273.15, 0, 0, 0, 0, 0, 0, 0]
conv_mult = [1, 1, 1, 1, 1, 1, 1 / 3600] # mm/h -> kg/m-2/s
valid_min = [180.00, 0, 0, 8e4, 0, 0, 0]
valid_max = [350.00, 100, 80, 1.5e5, 2500, 800, 0.2]

# ELM Variable names and units
# TBOT:     Air temperature at measurement (tower) height (K)
# RH:       Relative humidity at measurment height (%)
# WIND:     Wind speeed at measurement height (m/s)
# PSRF:     air pressure at surface  (Pa)
# FSDS:     Incoming Shortwave radiation  (W/m2)
# FLDS:     Incoming Longwave radiation   (W/m2)
# PRECTmms: Precipitation       (kg/m2/s)


os.system("mkdir -p " + outdir)


# Load the data
data = pd.read_csv(filename, index_col=0, header=0, parse_dates=True)
data.columns = [c.rstrip().lstrip() for c in data.columns.get_level_values(0)]

# drop 2.29
if leapdays:
    data = data.loc[~((data.index.month == 2) & (data.index.day == 29)), :]

data = data.loc[(data.index.year >= start_year) & (data.index.year <= end_year), :]


metdata = {}
for v in range(0, len(invars)):
    temp = data[invars[v]].astype(float).values * conv_mult[v] + conv_add[v]
    temp[(temp < valid_min[v]) | (temp > valid_max[v])] = np.nan
    metdata[outvars[v]] = list(temp)

out_fname = outdir + "/all_hourly.nc"
write_elm_met.bypass_format(
    out_fname,
    metdata,
    mylat,
    mylon,
    start_year,
    end_year,
    edge=0.1,
    time_offset=time_offset,
    calc_lw=calc_flds,
    zbot=measurement_height,
)