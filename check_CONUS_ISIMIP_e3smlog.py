# Check the lat & lon of the broken point
# Check which year it ended at

import os
import numpy as np
from glob import glob
import pandas as pd
import subprocess
import re

base_path = "/gpfs/wolf2/cades/cli185/proj-shared/zdr/ERW/output/UQ/pft1"

summary = []

#for ens in range(0, 2276):
for ens in range(63, 100):
    ens_str = f"{ens:04d}"
    pattern = os.path.join(base_path, f"*_conus_ICB20TRCNPRDCTCBC_erw_pft1_ens{ens_str}")
    dirs = glob(pattern + "/")

    if not dirs:
        continue
    folder = sorted(dirs)[-1]

    # Gap between restart file and last year
    sim_list = sorted(glob(os.path.join(folder, "run", "*.elm.h1.*.nc")))

    if len(sim_list) == 57:
        # the simulation completed without error
        continue

    if len(sim_list) > 0:
        last_sim = int(sim_list[-1].split(".")[-2].split("-")[0])
    else:
        last_sim = 2025

    rest_list = sorted(glob(os.path.join(folder, "run", "*.elm.r.*.nc")))
    if len(rest_list) > 0:
        last_rest = int(rest_list[-1].split(".")[-2].split("-")[0])
    else:
        last_rest = 2025

    nyear = last_sim - last_rest

    filename = os.path.join(folder, "run", "e3sm_log.txt")
    res =  subprocess.run(["tail", "-n 15", f"{filename}"], capture_output=True, text=True)
    lastlines = res.stdout

    for line in lastlines.split("\n"):
        if "Latdeg,Londeg" in line:
            res = re.split("\s+", line)
            lat,lon = float(res[2]),float(res[3])
            break

    summary.append([ens, lat, lon, nyear])

summary = pd.DataFrame(summary, columns = ['ens_id', 'lat', 'lon', 'nyear'])
summary.set_index('ens_id', inplace = True)

print(summary)