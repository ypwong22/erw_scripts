""" Visualize the outcome of ensemble_UIEF_analysis.py """
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches


def get_data_rank(depth, labels):
    data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble_UIEF_kge.csv'),
                    index_col = 0, header = [0,1])
    metric = data.drop(['param'], axis = 1, level = 0).loc[:, depth]

    # column-wise rank of the Kling-Gupta efficiency relative to other values in the same column
    # larger values = better
    data_rank = metric.drop(['total_dissolved_ca','total_dissolved_mg'], axis = 1).rank(axis = 0, ascending=False)

    # median rank over all the values
    median_data_rank = data_rank.median(axis = 1)

    # parameter values
    params = data['param']

    # assign group based on rank quantiles
    groups = pd.cut(median_data_rank, 
                    bins=[-np.inf] + [median_data_rank.quantile(q) for q in quantiles] + [np.inf],
                    labels=labels + ["all"])

    # also identify index of the best ensemble members
    best_ens = median_data_rank.index[np.argmin(median_data_rank)]
    print(f'Depth {depth} best ensemble is {best_ens}')

    return metric, params, groups, best_ens


# thresholds
quantiles = [0.05, 0.1, 0.3]
colors = ["#3182bd", "#9ecae1", "#bcbddc"]
labels = ["top5%", "top10%", "top30%"]


##############################################################################
# Plot the parameter values for the top x% performing ensemble members
##############################################################################
fig, axes = plt.subplots(1, 2, figsize=(10,6))

for i, depth in enumerate(['0-10cm','10-30cm']):

    _, params, groups, best_ens = get_data_rank(depth, labels)

    ax = axes.flat[i]

    for pos, cov in enumerate(params.columns):
        for i,label in enumerate(labels):
            bp = ax.boxplot(params.loc[groups == label, cov].dropna(), positions=[pos+i*0.3], widths=0.25, patch_artist=True)
            for patch in bp['boxes']:
                patch.set_facecolor(colors[i])

    h = ax.scatter(np.arange(0.3, params.shape[1]), params.loc[best_ens, :], color = 'r', label = 'Best', zorder = 10)

    legend_patches = [mpatches.Patch(color=c, label=l) for c, l in zip(colors, labels)]
    ax.legend(handles=[h] + legend_patches)
    ax.set_xticks(np.arange(0.3, params.shape[1]))
    ax.set_xticklabels(params.columns)
    ax.set_title('CEC parameter value range in \ntop x% performing ensemble members')
plt.tight_layout()
plt.savefig(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble_UIEF_kge.png'), dpi = 600., bbox_inches = 'tight')


##############################################################################
# Plot the RMSE values of each variable in the top 10% performing ensemble members
##############################################################################
fig, axes = plt.subplots(1, 2, figsize=(10,6))

for i, depth in enumerate(['0-10cm','10-30cm']):
    ax = axes.flat[i]

    metric, params, groups, best_ens = get_data_rank(depth, labels)
    metric = metric.loc[groups == 'top5%']
    metric = metric.sort_index(axis = 1)

    bp = ax.boxplot(metric, positions = range(metric.shape[1]), widths = 0.25, patch_artist=True)
    h = ax.scatter(np.arange(metric.shape[1]), metric.loc[best_ens, :], color = 'r', label = 'Best', zorder = 10)

    ax.set_title(depth)
    ax.set_xticks(np.arange(metric.shape[1]))
    ax.set_xticklabels(metric.columns, rotation = 90)

fig.text(-0.07, 0.5, 'Kling-Gupta Efficiency in top 5% performing ensemble members', rotation = 90, va = 'center')
plt.tight_layout()
plt.savefig(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'ensemble_UIEF_kge2.png'), dpi = 600., bbox_inches = 'tight')