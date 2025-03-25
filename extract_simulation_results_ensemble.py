""" Calculate and plot the annual total CDR for the entire PFT over CONUS """

import os
import numpy as np
from tqdm import tqdm
import itertools as it
import multiprocessing as mp
import matplotlib.pyplot as plt
from utils_ensemble import  get_areaTot_sequestration, ens_files, load_ensemble_tuples


path_out = os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble', 'pft1')


def process_ens_id(ens_id):
    """
    Process one ensemble id:
      - Get the file list with ens_files(ens_id)
      - For each file, compute the rate from the netCDF Dataset
    Returns a tuple (ens_id, rates) where rates is an array of computed values or None.
    """
    flist = ens_files(ens_id, 'pft1')
    if flist is None:
        return ens_id, None

    rates = get_areaTot_sequestration(flist)

    return ens_id, rates


recreate_txt = False
if recreate_txt:
    # Initialize the array to hold all rates
    all_rates = np.full([1100, 56], np.nan)

    # Create a pool with the number of available CPUs
    with mp.Pool(mp.cpu_count()) as pool:
        # Use imap for lazy evaluation and wrap it with tqdm for a progress bar
        results = list(tqdm(pool.imap(process_ens_id, range(1100)), total=1100))

    # Fill the all_rates array with the results
    for ens_id, rates in results:
        if rates is not None:
            all_rates[ens_id, :len(rates)] = rates

    np.savetxt(os.path.join(path_out, 'r_sequestration.txt'), all_rates)


##########################################################################################
ensemble_setups = load_ensemble_tuples()

frequency_list = ensemble_setups['Frequency'].unique()[1:]
rate_list = ensemble_setups['Rate (kg/m2)'].unique()[1:]
gra_list = ensemble_setups['Grain size (um)'].unique()[1:]
gcolor_list = ['#DAF7A6', '#FFC300', '#FF5733', '#C70039', '#900C3F', '#581845']
year_index = range(2025, 2081)

all_rates = np.loadtxt(os.path.join(path_out, 'r_sequestration.txt'))

# Control is plotted in all
for year in list(range(2025, 2051)):

    fig, axes = plt.subplots(5, 4, figsize = (20, 15), sharex = True, sharey = True)

    for i, (rate, freq) in enumerate(it.product(rate_list, frequency_list)):
        ax = axes.flat[i]

        for j, gra in enumerate(gra_list):

            ens_id = get_ens_id(ensemble_setups, gra, rate, year, freq)

            if ens_id is None or ens_id > 1099:
                continue

            ax.plot(year_index, all_rates[ens_id, :] - all_rates[0, :], '-o', 
                    color = gcolor_list[j])

        # ax.plot(year_index, all_rates[0, :], '-ok', label = 'Ctrl')

        if i < len(frequency_list):
            ax.set_title(f'Every {freq} year')
        if np.mod(i, len(frequency_list)) == 0:
            ax.set_ylabel(f'{rate} kg/m2')

    for j, gra in enumerate(gra_list):
        ax.plot(np.nan, np.nan, '-o', color = gcolor_list[j], label = f'{gra} um')
    ax.legend(loc = (-3, -0.5), ncol = 5)

    fig.suptitle('Cumulative CDR (tCO2) (relative to Control run)')

    fig.savefig(os.path.join(path_out, f'all_rates_{year}.png'), dpi = 600., bbox_inches = 'tight')
    plt.close()