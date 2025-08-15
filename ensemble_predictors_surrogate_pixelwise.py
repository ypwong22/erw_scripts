import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from itertools import product

from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score

from xgboost import XGBRegressor
from sklearn.ensemble import RandomForestRegressor # , GradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.neighbors import KNeighborsRegressor

from joblib import Parallel, delayed

path_root = os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble', 'pft1')

# ========== 1. Load Data ==========
df = pd.read_csv(os.path.join(path_root, 'combined_cdr.csv'), index_col = 0)  # replace with your actual path
df['Grain size (um)'] = 69.18 * df['Grain size (um)']**(-1.24) # Convert the grain size to SSA
df['Frequency'] = df.apply(
    lambda df: len(np.arange(df['Start year'], 2050 + 1, df['Frequency'])),
    axis = 1
) # Convert the frequency to the number of times applied
df = df.rename(columns={"Grain size (um)": "SSA (m2 g-1)", 'Frequency': 'Ntimes'})
X = df.iloc[:, :4].values  # Columns B-E: input variables
Y = df.iloc[:, 4:].values   # Columns F-end: output variables

# ========== 5. Do LHS on the treatment levels for Train/Test Split ==========
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

# Generate and then reduce
treatments = {key: df[key].drop_duplicates().values for key in df.columns[:4]}
n_allcomb = 16900

lhs_frac = 1
real_frac = 1
train_frac = 0.3
while real_frac > train_frac:
    lhs_n = int(np.ceil(n_allcomb * lhs_frac))
    df_lhs = lhs_from_levels(treatments, n=lhs_n, random_state=42, enforce_unique=True)
    mask = df.iloc[:, :4].apply(tuple, axis=1).isin(df_lhs.apply(tuple, axis=1))
    real_frac = np.sum(mask) / df.shape[0]
    lhs_frac *= 0.9
    print(f'Latin hypercube sampled fraction = {real_frac:.4f}')

idx_train = df.index[mask]
X_train = df.loc[idx_train, :].iloc[:, :4]
Y_train = df.loc[idx_train, :].iloc[:, 4:]

idx_test = df.index[~mask]
X_test = df.loc[idx_test, :].iloc[:, :4]
Y_test = df.loc[idx_test, :].iloc[:, 4:]

# ========== 4. Standardize Inputs & Outputs ==========
scaler_X = StandardScaler().fit(X_train)
X_train = scaler_X.transform(X_train)
X_test = scaler_X.transform(X_test)

scaler_Y = StandardScaler().fit(Y)
Y_train = scaler_Y.transform(Y_train)
Y_test = scaler_Y.transform(Y_test)

# ========== 9. Evaluate each model per output column in parallel ==========

def evaluate_single(icol, model_cls, params, name):
    """Fit a fresh model on Y[:,icol], compute un-scaled MSE & R²."""
    # instantiate and train
    model = model_cls(**params)
    model.fit(X_train, Y_train[:, icol])
    # predict
    y_pred = model.predict(X_test)
    # invert scaling for this column
    mean_i, scale_i = scaler_Y.mean_[icol], scaler_Y.scale_[icol]
    y_pred_uns = y_pred * scale_i + mean_i
    y_true_uns = Y_test[:, icol] * scale_i + mean_i
    ###y_pred_uns = np.expm1(y_pred).clip(min = 0.)
    ##y_true_uns = np.expm1(Y_test[:,icol]).clip(min = 0.)
    # metrics
    mse = mean_squared_error(y_true_uns, y_pred_uns)
    r2  = r2_score(y_true_uns, y_pred_uns)

    if hasattr(model, "feature_importances_"):
        return {"Model": name, "Column": icol, "MSE": mse, "R2": r2, "YPred": y_pred_uns, "Importance": model.feature_importances_}
    else:
        return {"Model": name, "Column": icol, "MSE": mse, "R2": r2, "YPred": y_pred_uns, "Importance": None}

# define base regressors (no multi‐output wrapper here)
model_defs = [
    ("Random Forest",      RandomForestRegressor,   {"n_estimators":200, "random_state":42}),
    ("Gradient Boosting",  XGBRegressor,            {"n_estimators":150, "random_state":42}),
    ("Ridge Regression",   Ridge,                   {"alpha":1.0}),
    ("KNN Regressor",      KNeighborsRegressor,     {"n_neighbors":5}),
]

# build one job per (model, column)
jobs = (
    delayed(evaluate_single)(icol, cls, params, name)
    for name, cls, params in model_defs
    for icol in range(Y_train.shape[1])
)

# run in parallel (–1 uses all cores; verbose=10 shows progress)
results = Parallel(n_jobs=-1, verbose=10)(jobs)

# collect performance into DataFrame and save
metrics_records = [
    {k: r[k] for k in ("Model", "Column", "MSE", "R2")}
    for r in results
]
results_df = pd.DataFrame(metrics_records)
results_df.to_csv(os.path.join(path_root, "model_performance.csv"), index=False)
print(f"Saved performance for {len(results_df)} fits to model_performance.csv")

# for Random Forest, collect importance values into DataFrame and save
imp_mats = np.full([df.shape[1]-4, 4], np.nan)
for r in results:
    m = r["Model"]
    j = r["Column"]    # column index in Y
    if m == "Random Forest":
        imp_mats[j, :] = r["Importance"]
imp_mats = pd.DataFrame(imp_mats, index = df.columns[4:], columns = df.columns[:4])
imp_mats.to_csv(os.path.join(path_root, "model_importances_random_forest.csv"), index=False)
print(f"Saved importance values for {len(results_df)} fits to model_importances_random_forest.csv")


# collect preedicted values into DataFrame and save
models = sorted({r["Model"] for r in results})
pred_mats = {m: pd.DataFrame(np.nan, index=df.index, columns = df.columns[4:]) for m in models}

for r in results:
    m = r["Model"]
    j = r["Column"]               # column index in Y
    yhat = np.asarray(r["YPred"]).ravel()
    if yhat.shape[0] != len(idx_test):
        raise ValueError(f"Unexpected prediction length for {m}, col {j}: {yhat.shape}")
    pred_mats[m].loc[idx_test, df.columns[4+j]] = yhat

# build DataFrames with same columns & index as original Y (test subset)
predictions_by_model = {
    m: pd.DataFrame(pred_mats[m], index=df.index, columns=df.columns[4:])
    for m in models
}

# save each model's predictions
# --- helper to make safe filenames from model names
def _slug(s: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in s).strip("_").lower()
for m, df_pred in predictions_by_model.items():
    pred_path = os.path.join(path_root, f"predictions_{_slug(m)}.csv")
    df_pred.to_csv(pred_path)
    print(f"Saved predictions for {m} to {pred_path}")


# ========== 10. Visualize performance across columns ==========

fig, axes = plt.subplots(2, 2, figsize=(16, 16))

for i, name in enumerate(results_df['Model'].unique()):
    ax = axes.flat[i]

    # plot MSE on the left y-axis
    sub = results_df[results_df['Model'] == name]
    ax.plot(sub['Column'], sub['MSE'], '-o', color = 'tab:blue')
    ax.set_xlabel("Output Column Index")
    ax.set_ylabel("Mean Squared Error")
    ax.set_title(name)
    ax.tick_params(axis='y', labelcolor='tab:blue')

    # create a second y-axis sharing the same x
    ax2 = ax.twinx()
    ax2.plot(sub['Column'], sub['R2'], '-s', color = 'tab:orange')
    ax2.set_ylabel("R$^2$ Score")
    ax2.tick_params(axis='y', labelcolor='tab:orange')

fig.tight_layout()
plt.savefig(os.path.join(path_root, 'surrogate_pixelwise.png'), dpi = 600., bbox_inches = 'tight')