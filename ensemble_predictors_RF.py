"""Test whether Random Forest can predict CDR from meteorological forcing and
   ensemble 0 hydrology.
   Split train-test in two ways: (1) in space, (2) in ensemble
"""
import os
import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from utils_ensemble import load_ensemble_tuples
import pickle
import itertools as it
from tqdm import tqdm
from sklearn_quantile import RandomForestQuantileRegressor

# Current bad ensemble members - skip these
skip = np.arange(1637, 2276)

# Reduce the spatial predictors
subset_predictors = True

# Retrain the random forest model
retrain = False

# The model object
rf = RandomForestRegressor(n_estimators=200, random_state=42, n_jobs=-1)
#rf = RandomForestQuantileRegressor(n_estimators=200, random_state=42,
#                                   n_jobs=-1, q = 0.5)

def collect_predictors(subset):
    """ Organize all the predictors into a 4D array (ensemble, lat, lon, variable)
        in original units.
    """
    # add ensemble-specific configurations
    ensemble_setups = load_ensemble_tuples().iloc[1:, :] # skip ensemble 0
    nens = ensemble_setups.shape[0]

    # additionally calculate the number of applications during the course of [start year, 2050]
    ensemble_setups['Ntimes'] = ensemble_setups.apply(
        lambda df: len(np.arange(df['Start year'], 2050 + 1, df['Frequency'])),
        axis = 1
    )

    nc = Dataset(os.path.join(os.environ['ZDR'], 'ERW', 'output', 'UQ', 'pft1', 
                              'predictors.nc'))
    data_vars = [var for var in nc.variables.keys() if not var in nc.dimensions]

    collect = []
    pred_names = []
    for var in data_vars:
        # skip the soil texture and base cation CEC due to unimportance
        if subset:
            if 'CEC_EFF' in var or 'PCT_SAND' in var or 'PCT_CLAY' in var:
                continue # skip the 

        if nc[var].ndim == 2:
            # the same for all ensemble members, broadcast to all ensemble members
            temp = nc[var][:,:]
            temp = np.where(temp.mask, np.nan, temp.data)
            temp = np.broadcast_to(temp[np.newaxis, ...], 
                                [nens, temp.shape[0], temp.shape[1]])
            collect.append(temp)
        else:
            # ensemble 0, broadcast to all ensemble members
            temp = nc[var][0, :, :]
            temp = np.where(temp.mask, np.nan, temp.data)
            temp = np.broadcast_to(temp[np.newaxis, ...],
                                   [nens, temp.shape[0], temp.shape[1]])
            collect.append(temp)
        pred_names.append(var)

    for var in ensemble_setups.columns:
        temp = np.broadcast_to(ensemble_setups[var].values.reshape(-1,1,1), 
                               temp.shape)
        collect.append(temp)

    # concatenate all the variables to 4D arrays
    collect = np.stack(collect, axis = 3)

    nc.close()

    pred_names = pred_names + list(ensemble_setups.columns)

    # delete skipped runs
    if len(skip) > 0:
        collect = np.delete(collect, skip-1, axis = 0)

    return collect, pred_names


def collect_predictand():
    """ Organize the simulated CDR (gC sequestered g basalt -1) into a 
        3D array (ensemble, lat, lon). Log-transform due to long-tail behavior.
    """
    nc = Dataset(os.path.join(os.environ['ZDR'], 'ERW', 'output', 'UQ', 'pft1', 
                              'r_sequestration_diff.nc'))
    cdr = (nc['r_sequestration_diff'][:, :, :] - nc['n2o_diff'][:, :, :] * 270) / nc['primary_added_diff'][:, :, :]
    cdr = np.where(cdr.mask, np.nan, cdr.data)
    nc.close()

    # delete skipped runs
    if len(skip) > 0:
        cdr = np.delete(cdr, skip-1, axis = 0)

    return cdr


def create_train_mask(spacing):
    """ Create boolean mask for the training set
        for regularly spaced thinning of lat-lon grids
    """
    lat = np.arange(23.25, 54.26, 0.5)
    lon = np.arange(234.75, 293.26, 0.5)

    bool_lat = np.full(len(lat), False)
    bool_lat[::spacing] = True
    bool_lon = np.full(len(lon), False)
    bool_lon[::spacing] = True

    return bool_lat, bool_lon, lat, lon


def evaluate(y_true, y_pred):
    """ Calculate standard metrics """
    correlation, _ = pearsonr(y_true, y_pred)
    bias = np.mean(y_pred - y_true)
    r2 = np.sqrt(r2_score(y_true, y_pred))
    return correlation, bias, r2


def size_test_spatial(func, x, y, spacing, x_names):
    """ Train on a subset of the spatial data points, but all the ensemble members.
        Test on the remaining spatial points. """
    nens = x.shape[0]

    bool_lat, bool_lon, lat, lon = create_train_mask(spacing)

    train_x = x[:,bool_lat,:,:][:,:,bool_lon,:].reshape(-1, x.shape[-1])
    train_y = y[:,bool_lat,:][:,:,bool_lon].reshape(-1)
    # skip NaN
    train_isnan = np.isnan(train_y) | np.any(np.isnan(train_x), axis = 1)
    train_x = train_x[~train_isnan, :]
    train_y = train_y[~train_isnan]

    test_x = x[:,~bool_lat,:,:][:,:,~bool_lon,:].reshape(-1, x.shape[-1])
    test_y = y[:,~bool_lat,:][:,:,~bool_lon].reshape(-1)
    # skip NaN
    test_isnan = np.isnan(test_y) | np.any(np.isnan(test_x), axis = 1)
    test_x = test_x[~test_isnan, :]
    test_y = test_y[~test_isnan]

    func.fit(train_x, train_y)

    y_pred = func.predict(test_x)

    corr, bias, r2 = evaluate(test_y, y_pred)

    # Create map of the predicted values
    y_pred_full = np.full(y.shape, np.nan)
    temp = np.full(len(test_isnan), np.nan)
    temp[~test_isnan] = y_pred
    temp2 = np.full([nens, sum(~bool_lat), len(bool_lon)], np.nan)
    temp2[:, :, ~bool_lon] = temp.reshape([nens, sum(~bool_lat), sum(~bool_lon)])
    y_pred_full[:, ~bool_lat, :] = temp2

    # Create map of the test values
    y_test_full = np.full(y.shape, np.nan)
    temp2 = np.full([nens, sum(~bool_lat), len(bool_lon)], np.nan)
    temp2[:, :, ~bool_lon] = y[:,~bool_lat,:][:,:,~bool_lon]
    y_test_full[:, ~bool_lat, :] = temp2

    results_setup = {
        "Sample %": len(train_y)/(len(test_y)+len(train_y))*100,
        "N Samples": len(train_y)
    }
    results_eval = {
        "Pearson r": corr,
        "Mean Bias": bias,
        "R$^2$": r2
    }
    results_importances = pd.Series(rf.feature_importances_, x_names)
    results_pred = {
        "Predicted": xr.DataArray(y_pred_full, coords = {'ensemble': np.arange(1,nens+1),
                                  'lat': lat, 'lon': lon}, dims = ['ensemble','lat','lon']),
        "True": xr.DataArray(y_test_full, coords = {'ensemble': np.arange(1,nens+1),
                                  'lat': lat, 'lon': lon}, dims = ['ensemble','lat','lon'])
    }

    return results_setup, results_eval, results_importances, results_pred



def size_test_ens(func, x, y, frac_train, x_names):
    """ Train on a subset of the ensemble members but all the spatial points. 
        Test on the rest. """
    nens = x.shape[0]
    _, _, lat, lon = create_train_mask(1)

    np.random.seed(43)
    bool_ens = np.random.rand(nens) <= (frac_train/100)

    train_x = x[bool_ens,:,:,:].reshape(-1, x.shape[-1])
    train_y = y[bool_ens,:,:].reshape(-1)
    # skip NaN
    train_isnan = np.isnan(train_y) | np.any(np.isnan(train_x), axis = 1)
    train_x = train_x[~train_isnan, :]
    train_y = train_y[~train_isnan]

    test_x = x[~bool_ens,:,:,:].reshape(-1, x.shape[-1])
    test_y = y[~bool_ens,:,:].reshape(-1)
    # skip NaN
    test_isnan = np.isnan(test_y) | np.any(np.isnan(test_x), axis = 1)
    test_x = test_x[~test_isnan, :]
    test_y = test_y[~test_isnan]

    func.fit(train_x, train_y)

    y_pred = func.predict(test_x)

    corr, bias, r2 = evaluate(test_y, y_pred)

    # Create map of the predicted values
    y_pred_full = np.full(y.shape, np.nan)
    temp = np.full(len(test_isnan), np.nan)
    temp[~test_isnan] = y_pred
    y_pred_full[~bool_ens,:,:] = temp.reshape([sum(~bool_ens), x.shape[1], x.shape[2]])

    # Create map of the test values
    y_test_full = np.full(y.shape, np.nan)
    y_test_full[~bool_ens,:,:] = y[~bool_ens,:,:]

    results_setup = {
        "Sample %": frac_train,
        "N Samples": sum(bool_ens),
    }
    results_eval = {
        "Pearson r": corr,
        "Mean Bias": bias,
        "R$^2$": r2
    }
    results_importances = pd.Series(rf.feature_importances_, x_names)
    results_pred = {
        "Predicted": xr.DataArray(y_pred_full, coords = {'ensemble': np.arange(1,nens+1),
                                  'lat': lat, 'lon': lon}, dims = ['ensemble','lat','lon']),
        "True": xr.DataArray(y_test_full, coords = {'ensemble': np.arange(1,nens+1),
                                  'lat': lat, 'lon': lon}, dims = ['ensemble','lat','lon'])
    }

    return results_setup, results_eval, results_importances, results_pred


def size_test_SVD(func, y, frac_train):
    """ 
    Reduce the spatial pattern by applying SVD on the [ensemble, spatial pts] matrix, 
        create a truncated matrix [ensemble, k << pts] by keeping only the leading k EOFs.
    Train on a subset of the ensemble members and predict the truncated matrix. 
    Test on the rest ensemble members.
    """
    nens = y.shape[0]

    # add ensemble-specific configurations
    ensemble_setups = load_ensemble_tuples().iloc[1:, :] # skip ensemble 0
    # additionally calculate the number of applications during the course of [start year, 2050]
    ensemble_setups['Ntimes'] = ensemble_setups.apply(
        lambda df: len(np.arange(df['Start year'], 2050 + 1, df['Frequency'])),
        axis = 1
    )
    if len(skip) > 0:
        ensemble_setups = ensemble_setups.drop(skip, axis = 0)

    _, _, lat, lon = create_train_mask(1)

    np.random.seed(43)
    bool_ens = np.random.rand(nens) <= (frac_train/100)

    # reshape in space
    y_2D = y.reshape(y.shape[0], -1)

    # spatial valid subset
    valid_pts = np.all(~np.isnan(y_2D), axis = 0)

    # subset to the spatial valid points
    y_2D = y_2D[:, valid_pts]

    # Need to transform the predictand because it is too skewed
    # this power transformation seems to get the distribution more symmetric
    y_2D = np.power(y_2D+1, 0.1)

    train_x = ensemble_setups.loc[bool_ens, :].values
    train_y = y_2D[bool_ens,:]
    test_x = ensemble_setups.loc[~bool_ens, :].values
    test_y = y_2D[~bool_ens,:]

    # do SVD
    # // Y = U S V'
    # y = np.matmul(np.matmul(U, np.diag(S)), Vh)
    U, S, Vh = np.linalg.svd(train_y, full_matrices = False)
    # reduced k-rank matrix
    k = 1
    train_ytil = np.matmul(Vh[:k,:], train_y.T).T
    train_ytil = train_ytil[:,0] # remove 

    func.fit(train_x, train_ytil)

    pred_ytil = func.predict(test_x)

    # transform back
    y_pred = np.matmul(pred_ytil.reshape(-1,1), Vh[:k,:])

    # Transform back the predicted values and test values
    test_y = np.power(test_y,10) - 1
    y_pred = np.power(y_pred,10) - 1

    # still need valid-pts only for evaluation
    test_y_ = test_y.reshape(-1)
    y_pred_ = y_pred.reshape(-1)
    filt = ~(np.isnan(test_y_) | np.isnan(y_pred_))
    corr, bias, r2 = evaluate(test_y_[filt], y_pred_[filt])

    # Create map of the predicted values
    y_pred_full = np.full(y.shape, np.nan)
    temp = np.full([y_pred.shape[0], len(valid_pts)], np.nan)
    temp[:, valid_pts] = y_pred
    y_pred_full[~bool_ens,:,:] = temp.reshape([temp.shape[0], y.shape[1], y.shape[2]])

    # Create map of the test values
    y_test_full = np.full(y.shape, np.nan)
    y_test_full[~bool_ens,:,:] = y[~bool_ens,:,:]

    results_setup = {
        "Sample %": frac_train,
        "N Samples": sum(bool_ens),
    }
    results_eval = {
        "Pearson r": corr,
        "Mean Bias": bias,
        "R$^2$": r2
    }

    x_names = ensemble_setups.columns
    results_importances = pd.Series(rf.feature_importances_, x_names)
    results_pred = {
        "Predicted": xr.DataArray(y_pred_full, coords = {'ensemble': np.arange(1,nens+1),
                                  'lat': lat, 'lon': lon}, dims = ['ensemble','lat','lon']),
        "True": xr.DataArray(y_test_full, coords = {'ensemble': np.arange(1,nens+1),
                                  'lat': lat, 'lon': lon}, dims = ['ensemble','lat','lon'])
    }

    return results_setup, results_eval, results_importances, results_pred


###########################################################################################
# Load data
###########################################################################################
predictors, pred_names = collect_predictors(subset_predictors)
cdr = collect_predictand()

if retrain:
    ###########################################################################################
    # Ablation test: spatial
    ###########################################################################################
    results_setup = []
    results_eval = []
    results_importance = []
    results_pred = []
    for spacing in [4,3,2]:
        setup, eval, importance, pred = size_test_spatial(rf, predictors, cdr, spacing, pred_names)
        results_setup.append(setup)
        results_eval.append(eval)
        results_importance.append(importance)
        results_pred.append(pred)
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                           'spatial_setup.pkl'), 'wb') as f:
        pickle.dump(results_setup, f)
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                           'spatial_eval.pkl'), 'wb') as f:
        pickle.dump(results_eval, f)
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                           'spatial_importance.pkl'), 'wb') as f:
        pickle.dump(results_importance, f, protocol=pickle.HIGHEST_PROTOCOL)
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                           'spatial_pred.pkl'), 'wb') as f:
        pickle.dump(results_pred, f)

    ###########################################################################################
    # Ablation test: ensemble
    ###########################################################################################
    results_setup = []
    results_eval = []
    results_importance = []
    results_pred = []
    for frac_train in [8, 16, 25]:
        setup, eval, importance, pred = size_test_ens(rf, predictors, cdr, frac_train, pred_names)
        results_setup.append(setup)
        results_eval.append(eval)
        results_importance.append(importance)
        results_pred.append(pred)
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                        'ens_setup.pkl'), 'wb') as f:
        pickle.dump(results_setup, f)
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                        'ens_eval.pkl'), 'wb') as f:
        pickle.dump(results_eval, f)
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                        'ens_importance.pkl'), 'wb') as f:
        pickle.dump(results_importance, f, protocol=pickle.HIGHEST_PROTOCOL)
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                        'ens_pred.pkl'), 'wb') as f:
        pickle.dump(results_pred, f)

    ###########################################################################################
    # Ablation test: SVD
    ###########################################################################################
    results_setup = []
    results_eval = []
    results_importance = []
    results_pred = []
    for frac_train in [8, 16, 25]:
        setup, eval, importance, pred = size_test_SVD(rf, cdr, frac_train)
        results_setup.append(setup)
        results_eval.append(eval)
        results_importance.append(importance)
        results_pred.append(pred)
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                           'svd_setup.pkl'), 'wb') as f:
        pickle.dump(results_setup, f)
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                           'svd_eval.pkl'), 'wb') as f:
        pickle.dump(results_eval, f)
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                           'svd_importance.pkl'), 'wb') as f:
        pickle.dump(results_importance, f, protocol=pickle.HIGHEST_PROTOCOL)
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                           'svd_pred.pkl'), 'wb') as f:
        pickle.dump(results_pred, f)


###########################################################################################
# Performance diagnostics
###########################################################################################
for style in ['spatial', 'ens', 'svd']:
    # Read the fraction of samples
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                           f'{style}_setup.pkl'), 'rb') as f:
        results_setup = pickle.load(f)
    xticks = [f'{r["N Samples"]}({r["Sample %"]:.1f}%)' for r in results_setup]

    #--------------------------------------------------------------------------------------
    # Plot the summary performance
    #--------------------------------------------------------------------------------------
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                           f'{style}_eval.pkl'), 'rb') as f:
        results_eval = pickle.load(f)

    # (1) Assemble results into a DataFrame
    df_eval = pd.DataFrame(results_eval)
    x = np.arange(df_eval.shape[0])
    # (2) Build the figure with three y-axes
    fig, ax1 = plt.subplots(figsize=(8, 6))

    ln1 = ax1.plot(x, df_eval["Pearson r"], color = '#1f77b4', marker="o", label="Pearson r")
    ax2 = ax1.twinx()
    ln2 = ax2.plot(x, df_eval["Mean Bias"], color = '#ff7f0e', marker="s", label="Mean Bias", linestyle="--")
    ax3 = ax1.twinx()
    ax3.spines["right"].set_position(("axes", 1.15))   # push spine outward
    ln3 = ax3.plot(x, df_eval["R$^2$"], color = '#2ca02c', marker="^", label="R$^2$", linestyle=":")

    ax1.set_xticks(x)
    ax1.set_xticklabels(xticks)
    ax1.set_xlabel("Training-set size")

    ax1.spines["left"].set_color('#1f77b4')
    ax1.tick_params(axis="y", colors='#1f77b4')
    ax1.set_ylabel("Pearson r", color='#1f77b4')

    ax2.spines["right"].set_color('#ff7f0e')
    ax2.tick_params(axis="y", colors='#ff7f0e')
    ax2.set_ylabel("Mean Bias")

    ax3.spines["left"].set_color('#2ca02c')
    ax3.tick_params(axis="y", colors='#2ca02c')
    ax3.set_ylabel("R$^2$")

    fig.tight_layout()
    fig.savefig(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results',
                             'ensemble', f'RF_{style}_summary.png'), dpi = 600., bbox_inches = 'tight')
    plt.close(fig)

    #--------------------------------------------------------------------------------------
    # Plot the importance of predictors
    #--------------------------------------------------------------------------------------
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                           f'{style}_importance.pkl'), 'rb') as f:
        results_importance = pickle.load(f)
    df_importance = pd.concat(results_importance, axis = 1, keys = xticks)
    fig, axes = plt.subplots(3, 1, sharex = True, sharey = True, figsize = (15, 10))
    for i, ax in enumerate(axes.flat):
        ax.bar(np.arange(df_importance.shape[0]), df_importance.iloc[:, i])
        ax.set_title(df_importance.columns[i])
        ax.set_xticks(np.arange(df_importance.shape[0]))
        ax.set_xticklabels(df_importance.index, rotation = 90)
        ax.set_xlim([-0.5, df_importance.shape[0]-0.5])
    fig.tight_layout()
    fig.savefig(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results',
                             'ensemble', f'RF_{style}_importance.png'), dpi = 600., bbox_inches = 'tight')
    plt.close(fig)

    #--------------------------------------------------------------------------------------
    # Plot the spatial performance
    #--------------------------------------------------------------------------------------
    with open(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble',
                           f'{style}_pred.pkl'), 'rb') as f:
        results_pred = pickle.load(f)
    fig, axes = plt.subplots(len(results_pred), 3, figsize = (24, 12))
    for i, pred in enumerate(results_pred):
        lat = pred['Predicted'].lat.values
        lon = pred['Predicted'].lon.values

        corr = np.full([len(lat), len(lon)], np.nan)
        bias = np.full([len(lat), len(lon)], np.nan)
        r2 = np.full([len(lat), len(lon)], np.nan)
        for j,k in tqdm(it.product(range(len(lat)), range(len(lon)))):
            if np.isnan(pred['Predicted'].values[1,j,k]):
                continue

            x = pred['Predicted'].values[:,j,k]
            y = pred['True'].values[:,j,k]
            filt = ~(np.isnan(x) | np.isnan(y))
            if sum(filt) > 0:
                x = x[filt]
                y = y[filt]
                corr[j,k], _ = pearsonr(x,y)
                bias[j,k] = np.mean(x) - np.mean(y)
                r2[j,k] = 1 - np.sum(np.power(x-y,2)) / np.sum(np.power(y-np.mean(y),2))

        ax = axes[i, 0]

        cf = ax.pcolormesh(lon, lat, corr, cmap = 'Spectral', 
                           norm = BoundaryNorm(np.linspace(0, 1, 21), ncolors=256, extend='both'))
        ax.set_title(f'{xticks[i]} Pearson r')
        plt.colorbar(cf)

        ax = axes[i, 1]
        cf = ax.pcolormesh(lon, lat, bias, cmap = 'Spectral',
                           norm = BoundaryNorm(np.linspace(-0.005, 0.0051, 21), ncolors=256, extend='both'))
        ax.set_title(f'{xticks[i]} Mean bias')
        plt.colorbar(cf)

        ax = axes[i, 2]
        cf = ax.pcolormesh(lon, lat, r2, cmap = 'Spectral',
                           norm = BoundaryNorm(np.linspace(-1, 1, 21), ncolors=256, extend='both'))
        ax.set_title(f'{xticks[i]} R$^2$')
        plt.colorbar(cf)

    fig.tight_layout()
    fig.savefig(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble', 
                             f'RF_{style}_maps.png'), dpi = 600., bbox_inches = 'tight')
    plt.close(fig)

    #--------------------------------------------------------------------------------------
    # Plot the ensemble-wise performance
    # Metrics: Pearsonr, mean absolute bias (%), R2 (%)
    #--------------------------------------------------------------------------------------
    fig, axes = plt.subplots(len(results_pred), 1, figsize = (24, 12), sharex = True, sharey = False)
    for i, pred in enumerate(results_pred):
        lat = pred['Predicted'].lat.values
        lon = pred['Predicted'].lon.values

        corr = np.full(pred['Predicted'].shape[0], np.nan)
        bias = np.full(pred['Predicted'].shape[0], np.nan)
        r2 = np.full(pred['Predicted'].shape[0], np.nan)
        for j in tqdm(range(len(corr))):
            x = pred['Predicted'].values[j,:,:].reshape(-1)
            y = pred['True'].values[j,:,:].reshape(-1)
            filt = ~(np.isnan(x) | np.isnan(y))
            if sum(filt) > 0:
                x = x[filt]
                y = y[filt]
                corr[j], _ = pearsonr(x,y)
                bias[j] = np.mean(x) - np.mean(y)
                r2[j] = r2_score(y,x)

        x = np.arange(pred['Predicted'].shape[0])

        ax1 = axes[i]
        ln1 = ax1.plot(x, corr, color = '#1f77b4', marker="o", label="Pearson r")
        ax2 = ax1.twinx()
        ln2 = ax2.plot(x, bias, color = '#ff7f0e', marker="s", label="Mean Bias", linestyle="--")
        ax3 = ax1.twinx()
        ax3.spines["right"].set_position(("axes", 1.15))   # push spine outward
        ln3 = ax3.plot(x, r2, color = '#2ca02c', marker="^", label="R$^2$", linestyle=":")

        ax1.set_xticks(x[::100])
        ax1.set_xlabel("Ensemble member")

        ax1.spines["left"].set_color('#1f77b4')
        ax1.tick_params(axis="y", colors='#1f77b4')
        ax1.set_ylabel("Pearson r", color='#1f77b4')

        ax2.spines["right"].set_color('#ff7f0e')
        ax2.tick_params(axis="y", colors='#ff7f0e')
        ax2.set_ylabel("Mean Bias")

        ax3.spines["left"].set_color('#2ca02c')
        ax3.tick_params(axis="y", colors='#2ca02c')
        ax3.set_ylabel("R$^2$")

        ax1.set_title(xticks[i])

    fig.tight_layout()
    fig.savefig(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble', 
                             f'RF_{style}_ensemble.png'), dpi = 600., bbox_inches = 'tight')
    plt.close(fig)
