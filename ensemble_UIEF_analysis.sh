#!/bin/bash -e

#SBATCH -t 12:00:00
#SBATCH -J ensemble_RF
#SBATCH --nodes=1
#SBATCH -A CLI185
#SBATCH -p batch_ccsi

# Change to submission directory
cd ~/Git/erw_scripts/

# Activate your conda environment
source activate olmt

# Run the script
python ensemble_UIEF_analysis.py
