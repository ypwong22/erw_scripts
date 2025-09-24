import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from itertools import product

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score

from xgboost import XGBRegressor
from sklearn.ensemble import RandomForestRegressor # , GradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.neighbors import KNeighborsRegressor

from joblib import Parallel, delayed


def get_path_root(which_pft) -> str:
    return os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble', which_pft)


def load_and_prepare_data(path_root: str) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    df = pd.read_csv(os.path.join(path_root, 'combined_cdr.csv'), index_col=0)  # replace with your actual path
    df['Grain size (um)'] = 69.18 * df['Grain size (um)']**(-1.24) # Convert the grain size to SSA
    df['Frequency'] = df.apply(
        lambda df: len(np.arange(df['Start year'], 2050 + 1, df['Frequency'])),
        axis=1
    ) # Convert the frequency to the number of times applied
    df = df.rename(columns={"Grain size (um)": "SSA (m2 g-1)", 'Frequency': 'Ntimes'})

    # Columns 0-3: input variables
    # Columns 4-end: output variables
    return df


def _latin_hypercube(n_samples: int, n_dims: int, random_state: int | None = None) -> np.ndarray:
    """
    Classic LHS in [0,1) for n_dims dimensions.
    Returns an (n_samples, n_dims) array.
    """
    rng = np.random.default_rng(random_state)
    U = np.empty((n_samples, n_dims), dtype=float)
    # For each dimension: random permutation of strata {0..n_samples-1} plus uniform jitter inside each stratum
    for j in range(n_dims):
        strata = rng.permutation(n_samples)
        U[:, j] = (strata + rng.random(n_samples)) / n_samples
    return U


def lhs_from_levels(
    levels: dict[str, np.ndarray | list],
    n: int,
    random_state: int | None = None,
    enforce_unique: bool = False,
) -> pd.DataFrame:
    """
    Latin Hypercube sample over *discrete* treatment levels.

    Parameters
    ----------
    levels : dict[name -> list/array]
        Discrete levels for each treatment (categorical or numeric).
    n : int
        Number of rows to sample.
    random_state : int | None
        Seed for reproducibility.
    enforce_unique : bool
        If True, deduplicate and top-up from the full factorial (without replacement).

    Returns
    -------
    pd.DataFrame
    """
    if n <= 0:
        return pd.DataFrame({k: [] for k in levels})

    keys = list(levels.keys())
    arrays = [np.asarray(levels[k], dtype=object) for k in keys]
    n_dims = len(keys)

    U = _latin_hypercube(n_samples=n, n_dims=n_dims, random_state=random_state)

    # Map each LHS draw in [0,1) to a discrete level index for that variable
    cols = {}
    for j, k in enumerate(keys):
        k_levels = arrays[j]
        idx = np.minimum((U[:, j] * len(k_levels)).astype(int), len(k_levels) - 1)
        cols[k] = k_levels[idx]

    df = pd.DataFrame(cols)

    if not enforce_unique:
        return df

    # Enforce uniqueness by topping up from the full factorial if needed
    df_unique = df.drop_duplicates(ignore_index=True)
    need = n - len(df_unique)
    if need <= 0:
        return df_unique.iloc[:n].reset_index(drop=True)

    # Build the full factorial (careful if gigantic; your case is fine: 5*5*26*26 = 16,900)
    full = pd.DataFrame(list(product(*arrays)), columns=keys)

    # Remove combos already present
    merged = full.merge(df_unique.assign(_present=1), how="left", on=keys)
    remaining = merged[merged["_present"].isna()].drop(columns="_present")

    take = min(need, len(remaining))
    if take > 0:
        rng = np.random.default_rng(random_state)
        extra = remaining.sample(n=take, replace=False, random_state=int(rng.integers(0, 2**31 - 1)))
        df_unique = pd.concat([df_unique, extra], ignore_index=True)

    return df_unique.reset_index(drop=True).iloc[:min(n, len(df_unique))]


def lhs_train_test_split(df: pd.DataFrame, path_root: str, train_frac: float = 0.3):
    treatments = {key: df[key].drop_duplicates().values for key in df.columns[:4]}
    n_allcomb = 16900

    df_na = df.loc[df[df.columns[4:]].isna().all(axis=1), :]
    df_valid = df.dropna(axis = 0, how = 'all', subset = df.columns[4:]) # drop unused rows

    # Because the search space is not the full space, use LHS to sample
    # the full space and remove those that are not in the search space. 
    # Iteratively find the sample fraction needed for the full space.
    lhs_frac = 1
    real_frac = 1
    while real_frac > train_frac:
        lhs_n = int(np.ceil(n_allcomb * lhs_frac))
        df_lhs = lhs_from_levels(treatments, n=lhs_n, random_state=42, enforce_unique=True)
        # Boolean mask for whether the ensemble IDs are in the sampled space
        mask = df_valid.iloc[:, :4].apply(tuple, axis=1).isin(df_lhs.apply(tuple, axis=1))
        real_frac = np.sum(mask) / df_valid.shape[0]
        lhs_frac *= 0.9
        print(f'Latin hypercube sampled fraction = {real_frac:.4f}')

    idx_train = df_valid.index[mask]
    X_train_df = df.loc[idx_train, :].iloc[:, :4]
    Y_train_df = df.loc[idx_train, :].iloc[:, 4:]
    # Save the training ensemble IDs to CSV
    X_train_df.to_csv(os.path.join(path_root, 'training_samples.csv'))

    idx_test = df_valid.index[~mask].append(df_na.index)
    X_test_df = df.loc[idx_test, :].iloc[:, :4]
    Y_test_df = df.loc[idx_test, :].iloc[:, 4:]

    return idx_train, idx_test, X_train_df, Y_train_df, X_test_df, Y_test_df


def standardize_data(X_train_df: pd.DataFrame, X_test_df: pd.DataFrame, Y_train_df: pd.DataFrame, Y_test_df: pd.DataFrame):
    scaler_X = StandardScaler().fit(X_train_df)
    X_train = scaler_X.transform(X_train_df)
    X_test = scaler_X.transform(X_test_df)

    scaler_Y = StandardScaler().fit(Y_train_df)
    Y_train = scaler_Y.transform(Y_train_df)
    Y_test = scaler_Y.transform(Y_test_df)

    return X_train, X_test, Y_train, Y_test, scaler_X, scaler_Y


def evaluate_single(icol, model_cls, params, name, X_train, Y_train, X_test, Y_test, scaler_Y):
    """Fit a fresh model on Y[:,icol], compute un-scaled MSE & R²."""
    model = model_cls(**params)
    model.fit(X_train, Y_train[:, icol])
    y_pred = model.predict(X_test)
    
    mean_i, scale_i = scaler_Y.mean_[icol], scaler_Y.scale_[icol] # invert scaling for this column
    y_pred_uns = y_pred * scale_i + mean_i
    y_true_uns = Y_test[:, icol] * scale_i + mean_i
    ###y_pred_uns = np.expm1(y_pred).clip(min = 0.)
    ##y_true_uns = np.expm1(Y_test[:,icol]).clip(min = 0.)

    # metrics
    filt = ~np.isnan(y_true_uns)
    mse = mean_squared_error(y_true_uns[filt], y_pred_uns[filt])
    r2  = r2_score(y_true_uns[filt], y_pred_uns[filt])

    if hasattr(model, "feature_importances_"):
        return {"Model": name, "Column": icol, "MSE": mse, "R2": r2, "YPred": y_pred_uns, "Importance": model.feature_importances_}
    else:
        return {"Model": name, "Column": icol, "MSE": mse, "R2": r2, "YPred": y_pred_uns, "Importance": None}


def build_and_save_metrics(results, path_root):
    """ collect performance into DataFrame and save """
    metrics_records = [
        {k: r[k] for k in ("Model", "Column", "MSE", "R2")}
        for r in results
    ]
    results_df = pd.DataFrame(metrics_records)
    results_df.to_csv(os.path.join(path_root, "model_performance.csv"), index=False)
    print(f"Saved performance for {len(results_df)} fits to model_performance.csv")
    return results_df


def save_random_forest_importances(results, df, path_root):
    imp_mats = np.full([df.shape[1]-4, 4], np.nan)
    for r in results:
        m = r["Model"]
        j = r["Column"]    # column index in Y
        if m == "Random Forest":
            imp_mats[j, :] = r["Importance"]
    imp_mats = pd.DataFrame(imp_mats, index = df.columns[4:], columns = df.columns[:4])
    imp_mats.to_csv(os.path.join(path_root, "model_importances_random_forest.csv"), index=False)
    print(f"Saved importance values for {imp_mats.shape[0]} outputs to model_importances_random_forest.csv")


def _slug(s: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in s).strip("_").lower()


def save_predictions(results, df, idx_test, path_root):
    models = sorted({r["Model"] for r in results})
    pred_mats = {m: pd.DataFrame(np.nan, index=df.index, columns = df.columns[4:]) for m in models}

    for r in results:
        m = r["Model"]
        j = r["Column"]               # column index in Y
        yhat = np.asarray(r["YPred"]).ravel()
        if yhat.shape[0] != len(idx_test):
            raise ValueError(f"Unexpected prediction length for {m}, col {j}: {yhat.shape}")
        pred_mats[m].loc[idx_test, df.columns[4+j]] = yhat

    for m, df_pred in pred_mats.items():
        pred_path = os.path.join(path_root, f"predictions_{_slug(m)}.csv")
        df_pred.to_csv(pred_path)
        print(f"Saved predictions for {m} to {pred_path}")


def plot_performance(results_df, path_root):
    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    for i, name in enumerate(results_df['Model'].unique()):
        ax = axes.flat[i]
        sub = results_df[results_df['Model'] == name]
        ax.plot(sub['Column'], sub['MSE'], '-o', color = 'tab:blue')
        ax.set_xlabel("Output Column Index")
        ax.set_ylabel("Mean Squared Error")
        ax.set_title(name)
        ax.tick_params(axis='y', labelcolor='tab:blue')
        ax2 = ax.twinx()
        ax2.plot(sub['Column'], sub['R2'], '-s', color = 'tab:orange')
        ax2.set_ylabel("R$^2$ Score")
        ax2.tick_params(axis='y', labelcolor='tab:orange')
    fig.tight_layout()
    plt.savefig(os.path.join(path_root, 'surrogate_pixelwise.png'), dpi = 600., bbox_inches = 'tight')


if __name__ == "__main__":
    path_root = get_path_root('pft1')

    # ========== 1. Load Data ==========
    df = load_and_prepare_data(path_root)

    # ========== 5. Do LHS on the treatment levels for Train/Test Split ==========
    idx_train, idx_test, X_train_df, Y_train_df, X_test_df, Y_test_df = lhs_train_test_split(df, path_root, train_frac=0.3)
    X_train, X_test, Y_train, Y_test, _, scaler_Y = standardize_data(X_train_df, X_test_df, Y_train_df, Y_test_df)

    # ========== 9. Evaluate each model per output column in parallel ==========
    # define base regressors (no multi‐output wrapper here)
    model_defs = [
        ("Random Forest",      RandomForestRegressor,   {"n_estimators":200, "random_state":42}),
        ("Gradient Boosting",  XGBRegressor,            {"n_estimators":150, "random_state":42}),
        ("Ridge Regression",   Ridge,                   {"alpha":1.0}),
        ("KNN Regressor",      KNeighborsRegressor,     {"n_neighbors":5}),
    ]
    # build one job per (model, column)
    jobs = (
        delayed(evaluate_single)(icol, cls, params, name, X_train, Y_train, X_test, Y_test, scaler_Y)
        for name, cls, params in model_defs
        for icol in range(Y_train.shape[1])
    )
    # run in parallel (–1 uses all cores; verbose=10 shows progress)
    results = Parallel(n_jobs=-1, verbose=10)(jobs)

    results_df = build_and_save_metrics(results, path_root)
    save_random_forest_importances(results, df, path_root)
    save_predictions(results, df, idx_test, path_root)

    # ========== 10. Visualize performance across columns ==========
    plot_performance(results_df, path_root)
