""" Visualize the outcome of ensemble_UIEF_analysis.py """
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble_UIEF_rmse.csv'),
                   index_col = 0, header = [0,1])

# column-wise rank of the RMSE relative to other values in the same column
data_rank = data.drop(['param'], axis = 1, level = 0).rank(axis = 0, ascending=True)

# median rank over all the values
median_data_rank = data_rank.median(axis = 1)

# parameter values
params = data['param']

# thresholds
quantiles = [0.1, 0.3, 0.5]
labels = ["top10%", "top30%", "top50%"]
colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

# assign group based on rank quantiles
groups = pd.cut(median_data_rank, 
                bins=[-np.inf] + [median_data_rank.quantile(q) for q in quantiles] + [np.inf],
                labels=labels + ["all"])

##############################################################################
# Plot the parameter values for the top x% performing ensemble members
##############################################################################
fig, ax = plt.subplots(figsize=(10,6))

positions = []
pos = 1

for cov in params.columns:
    for i, label in enumerate(labels):
        bp = ax.boxplot(params.loc[groups == label, cov].dropna(), positions=[pos+i*0.3], widths=0.25, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor(colors[i])

    positions.append(pos+0.3)
    pos += len(labels) + 1

legend_patches = [mpatches.Patch(color=c, label=l) for c, l in zip(colors, ["top10%", "top30%", "top50%"])]
ax.legend(handles=legend_patches)
ax.set_xticks(positions)
ax.set_xticklabels(params.columns)
ax.set_title('CEC parameter value range in \ntop x% performing ensemble members')
plt.tight_layout()
plt.savefig(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble_UIEF_rmse.png'), dpi = 600., bbox_inches = 'tight')


##############################################################################
# Plot the RMSE values of each variable in the top 10% performing ensemble members
##############################################################################
rmse = data.drop(['param'], axis = 1, level = 0)
rmse = rmse.loc[groups == 'top10%']
rmse.columns = rmse.columns.reorder_levels([1,0]).remove_unused_levels()
rmse = rmse.sort_index(axis = 1)

colors = ["#1f77b4", "#2ca02c"]

fig, axes = plt.subplots(2, 3, figsize = (10, 6))
for i, var in enumerate(rmse.columns.levels[0]):
    subset = rmse[var]
    ax = axes.flat[i]

    for j, depth in enumerate(subset.columns):
        bp = ax.boxplot(subset[depth], positions=[j], widths=0.25, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor(colors[j])

    ax.set_title(var)
    ax.set_xticks([])
    if ax == 0:
        legend_patches = [mpatches.Patch(color=c, label=l) for c, l in zip(colors, subset.columns)]
        ax.legend(handles=legend_patches)
axes.flat[-1].axis('off')
fig.text(-0.07, 0.5, 'RMSE values in top x% performing ensemble members', rotation = 90, va = 'center')
plt.tight_layout()
plt.savefig(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble_UIEF_rmse2.png'), dpi = 600., bbox_inches = 'tight')