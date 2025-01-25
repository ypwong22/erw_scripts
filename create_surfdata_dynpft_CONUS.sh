#!/bin/csh
#SBATCH --time=24:0:00
#SBATCH -J CONUSdynpftREPLACE
#SBATCH --nodes=1
#SBATCH -A CLI185
#SBATCH -p batch_ccsi
#SBATCH --ntasks-per-node 1

cd ${HOME}/Git/erw_scripts
conda activate myCondaEnv
python temp_create_surfdata_dynpft_CONUS_REPLACE.py