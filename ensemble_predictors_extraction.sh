#!/bin/bash -e

#SBATCH -t 24:00:00
#SBATCH -J process_ERW_predictors
#SBATCH --nodes=1
#SBATCH -A CLI185
#SBATCH -p batch_ccsi

# Change to submission directory
cd ~/Git/erw_scripts/

# Activate your conda environment
source activate myCondaEnv

# Run the script
python ensemble_predictors_extraction.py
