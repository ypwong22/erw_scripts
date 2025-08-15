import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score

from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.neighbors import KNeighborsRegressor
from sklearn.multioutput import MultiOutputRegressor

from tensorflow import keras
from keras.models import Sequential
from keras.layers import Dense, Input
from keras.optimizers import Adam

#from tensorflow.keras.models import Sequential
#from tensorflow.keras.layers import Dense, Input
#from tensorflow.keras.optimizers import Adam

# ========== 1. Load Data ==========
df = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble', 'pft1', 'combined.csv'))  # replace with your actual path
X = df.iloc[:, 1:5]  # Columns B-E: input variables
Y = df.iloc[:, 5:]   # Columns F-end: output variables

# ========== 2. Filter Samples by Output Range ==========
mask = (Y <= 0.2) & (Y >= -0.2)
rows_to_keep = mask.all(axis=1)
rows_to_remove = ~rows_to_keep

X_clean = X[rows_to_keep].reset_index(drop=True)
Y_clean = Y[rows_to_keep].reset_index(drop=True)
X_removed = X[rows_to_remove].reset_index(drop=True)

# ========== 3. Visualize Removed Inputs ==========
sns.pairplot(X_removed)
plt.suptitle("Input Variable Distributions for Removed Samples", y=1.02)
plt.tight_layout()
plt.show()

# ========== 4. Standardize Inputs & Outputs ==========
scaler_X = StandardScaler().fit(X_clean)
X_scaled = scaler_X.transform(X_clean)

scaler_Y = StandardScaler().fit(Y_clean)
Y_scaled = scaler_Y.transform(Y_clean)

# ========== 5. Split Train/Test ==========
X_train, X_test, Y_train, Y_test = train_test_split(
    X_scaled, Y_scaled, test_size=0.2, random_state=42
)

# ========== 6. PCA on Outputs ==========
pca_full = PCA().fit(Y_train)
cum_var = np.cumsum(pca_full.explained_variance_ratio_)

# Plot cumulative explained variance
plt.figure(figsize=(6, 4))
plt.plot(cum_var)
plt.axhline(0.95, linestyle='--', color='gray')
plt.xlabel('Number of Components')
plt.ylabel('Cumulative Explained Variance')
plt.title('PCA Explained Variance')
plt.grid(True)
plt.tight_layout()
plt.show()

# Choose k where 95% variance is retained
k_opt = np.searchsorted(cum_var, 0.99) + 1
print(f"Selected k = {k_opt} components for 95% variance")

pca = PCA(n_components=k_opt).fit(Y_train)
Y_train_pc = pca.transform(Y_train)
Y_test_pc = pca.transform(Y_test)

# ========== 7. Define Evaluation Function ==========
def evaluate_model(model, name):
    model.fit(X_train, Y_train_pc)
    Y_pc_pred = model.predict(X_test)
    Y_pred = scaler_Y.inverse_transform(pca.inverse_transform(Y_pc_pred))
    Y_true = scaler_Y.inverse_transform(Y_test)
    mse = mean_squared_error(Y_true, Y_pred)
    r2 = r2_score(Y_true, Y_pred, multioutput='variance_weighted')
    return {"Model": name, "MSE": mse, "R2": r2}

# ========== 8. ML Models ==========
models = [
    (MultiOutputRegressor(RandomForestRegressor(n_estimators=100, random_state=42)), "Random Forest"),
    (MultiOutputRegressor(GradientBoostingRegressor(n_estimators=100, random_state=42)), "Gradient Boosting"),
    (MultiOutputRegressor(Ridge(alpha=1.0)), "Ridge Regression"),
    (MultiOutputRegressor(KNeighborsRegressor(n_neighbors=5)), "KNN Regressor"),
]

results = []
for model, name in models:
    results.append(evaluate_model(model, name))

# ========== 9. Neural Network ==========
size = 32 # 64
nn = Sequential([
    Input(shape=(X_train.shape[1],)),
    Dense(size, activation='relu'),
    Dense(size, activation='relu'),
    Dense(k_opt)
])
nn.compile(optimizer=Adam(1e-3), loss='mse')
nn.fit(X_train, Y_train_pc, epochs=30, batch_size=32, validation_split=0.1, verbose=0)

# Predict and evaluate NN
Y_pc_pred_nn = nn.predict(X_test)
Y_pred_nn = scaler_Y.inverse_transform(pca.inverse_transform(Y_pc_pred_nn))
Y_true_nn = scaler_Y.inverse_transform(Y_test)
mse_nn = mean_squared_error(Y_true_nn, Y_pred_nn)
r2_nn = r2_score(Y_true_nn, Y_pred_nn, multioutput='variance_weighted')
results.append({"Model": "Neural Network", "MSE": mse_nn, "R2": r2_nn})

# ========== 10. Compare Results ==========
results_df = pd.DataFrame(results)
print("\nModel Comparison:\n", results_df.sort_values("R2", ascending=False))

# Optional: Plot results
sns.barplot(data=results_df.melt(id_vars="Model", var_name="Metric", value_name="Value"),
            x="Model", y="Value", hue="Metric")
plt.title("Model Comparison (PCA + ML on 4 Inputs → 1542 Outputs)")
plt.xticks(rotation=30)
plt.tight_layout()
plt.show()
plt.savefig(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble', 'pft1', 'surrogate.png'), dpi = 600., bbox_inches = 'tight')