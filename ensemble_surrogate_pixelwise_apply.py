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
    return df


def split_and_standardize_data(df: pd.DataFrame):
    # Columns 0-3: input variables
    # Columns 4-end: output variables
    idx_train = df.index[~df.iloc[:, 4].isna()]
    idx_test = df.index[df.iloc[:, 4].isna()]

    X_train_df = df.loc[idx_train, :].iloc[:, :4]
    X_test_df = df.loc[idx_test, :].iloc[:, :4]

    scaler_X = StandardScaler().fit(X_train_df)
    X_train = scaler_X.transform(X_train_df)
    X_test = scaler_X.transform(X_test_df)

    Y_train_df = df.loc[idx_train, :].iloc[:, 4:]

    scaler_Y = StandardScaler().fit(Y_train_df)
    Y_train = scaler_Y.transform(Y_train_df)

    return X_train, Y_train, idx_train, idx_test, X_test, scaler_Y



def evaluate_single(icol, X_train, Y_train, X_test, scaler_Y):
    """Fit a fresh model on Y[:,icol], predict the missing values."""
    params = {"n_estimators":200, "random_state":42}
    model = RandomForestRegressor(**params)
    model.fit(X_train, Y_train[:, icol])
    y_pred = model.predict(X_test)

    mean_i, scale_i = scaler_Y.mean_[icol], scaler_Y.scale_[icol] # invert scaling for this column
    y_pred_uns = y_pred * scale_i + mean_i

    return {"Column": icol, "YPred": y_pred_uns}


def collect_prediction(results, idx_test, df_columns):
    pred_mats = pd.DataFrame(np.nan, index=idx_test, columns = df_columns[4:])

    for r in results:
        j = r["Column"]               # column index in Y
        yhat = np.asarray(r["YPred"]).ravel()
        if yhat.shape[0] != len(idx_test):
            raise ValueError(f"Unexpected prediction length for col {j}: {yhat.shape}")
        pred_mats.loc[:, df_columns[4+j]] = yhat

    return pred_mats


def save_predictions(df, df_preds, path_root):
    df_copy = df.copy()
    df_copy['is_predicted'] = df_copy.iloc[:,4].isna()
    df_copy.loc[df_copy['is_predicted'], df.columns[4:]] = df_preds

    pred_path = os.path.join(path_root, f"predictions.csv")
    df_copy.to_csv(pred_path)
    print(f"Saved predictions to {pred_path}")


if __name__ == "__main__":
    path_root = get_path_root('pft15')

    # ========== 1. Load Data ==========
    df = load_and_prepare_data(path_root)

    # ========== 5. Identify the training data (not NaNs) and to-be-predicted data (NaN) ==========
    X_train, Y_train, idx_train, idx_test, X_test, scaler_Y = split_and_standardize_data(df)

    # ========== 9. Predict the un-simulated data ==========
    jobs = (
        delayed(evaluate_single)(icol, X_train, Y_train, X_test, scaler_Y)
        for icol in range(Y_train.shape[1])
    )
    # run in parallel (–1 uses all cores; verbose=10 shows progress)
    results = Parallel(n_jobs=-1, verbose=10)(jobs)

    pred_mats = collect_prediction(results, idx_test, df.columns)
    save_predictions(df, pred_mats, path_root)
 