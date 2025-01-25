#!/bin/csh
#SBATCH --time=24:0:00
#SBATCH -J CONUSdynpft4
#SBATCH --nodes=1
#SBATCH -A CLI185
#SBATCH -p batch_ccsi
#SBATCH --ntasks-per-node 1

cd ${HOME}/Git/erw_scripts
conda activate myCondaEnv
python temp_create_surfdata_dynpft_CONUS_4.py