import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from contacTreesIE.plotting.plot_ct_vs_bt import declutter_violin_plot, add_grid_lines


SUBPLOT_TITLE_ARGS = dict(fontsize=14, fontweight='semibold')

PARAMETER_NAME_MAPPING = {
    'clockRate': 'clock rate',
    'rootHeight': 'root height',
    'treeLength': 'tree length $L$',
    'pairedTreeLength': 'paired tree length $L^{(2)}$',
    'conversionRate': 'contact rate $\gamma$',
    'expectedConversions': 'expected cont. edges $\Gamma$',
    'convCount': 'cont. edge count $|\mathcal{C}|$',
    'meanConvHeight': 'mean cont. edge height',
    'pMove': 'pBorrow $\\beta$',
    'moveCount': 'borrow count $\Vert Z \Vert_0$',
    'movesPerConv': 'borrowings per cont. edge',
}


def plot_errors(bt_stats, ct_stats, groundtruth,
                stat_id='', stat_label='', stat_range=(None,None),
                diff_to_gt=False, ax=None):
    ax = ax or plt.gca()
    stat_mean_id = f'{stat_id}.mean'

    if diff_to_gt:
        bt_stats['error'] = bt_stats[stat_mean_id].to_numpy() - groundtruth[stat_id].to_numpy()
        ct_stats['error'] = ct_stats[stat_mean_id].to_numpy() - groundtruth[stat_id].to_numpy()

    bt_stats['mode'] = 'noCT'
    ct_stats['mode'] = 'CT'
    stats = pd.concat([bt_stats, ct_stats])

    sns.violinplot(
        data=stats,
        x='mode',
        y='error' if diff_to_gt else stat_mean_id,
        inner='quartile',
        cut=1,
        ax=ax,
    )

    print(f'BT {stat_label} mean:\t {bt_stats["error" if diff_to_gt else stat_mean_id].mean()}')
    print(f'CT {stat_label} mean:\t {ct_stats["error" if diff_to_gt else stat_mean_id].mean()}')

    ax.tick_params(axis='x', which='major', labelsize=13)
    ax.set_xlim(-0.7, 1.5)
    ax.set_ylim(*stat_range)
    ax.set_title(stat_label, fontsize=14, pad=20, fontweight='semibold')
    ax.set_xlabel('')
    ax.set_ylabel('')

    # ax.axhline(0.0, color='lightgrey', zorder=0)
    ax.axhline(0.0, color='darkgrey', zorder=0)

    add_grid_lines(ax, grid_values=ax.get_yticks(), clip_on=True, zorder=0)
    declutter_violin_plot(ax)


def plot_coverage(stats, groundtruth,
                  stat_id='', stat_label='', stat_range=(None,None),
                  diff_to_gt=False, ticks=None, lim=None,
                  show_title=False, ax=None):

    ax = ax or plt.gca()

    mid = stats[stat_id + '.mean'].astype(float).to_numpy()
    lo = stats[stat_id + '.95%HPDlo'].astype(float).to_numpy()
    up = stats[stat_id + '.95%HPDup'].astype(float).to_numpy()
    gt = groundtruth[stat_id.partition('.')[0]].to_numpy()

    covered_colors = np.array(['red', 'grey'])
    covered_alpha = np.array([1.0, 0.6])
    covered = np.logical_and(lo <= gt, gt <= up).astype(int)
    color = covered_colors[covered]
    alpha = covered_alpha[covered]

    for i in range(len(mid)):
        ax.plot([gt[i], gt[i]], [lo[i], up[i]], c=color[i], lw=2, alpha=alpha[i], zorder=alpha[i])

    ax.scatter(gt, mid, color='k', s=20, marker='_', zorder=2)

    xy_min = min([*up, *gt])
    xy_max = max([*up, *gt])
    ax.plot([0.0, xy_max], [0.0, xy_max], c='lightgrey', zorder=0)

    if lim:
        ax.set_xlim(*lim)
        ax.set_ylim(*lim)
    else:
        ax.set_xlim(xy_min, xy_max)
        ax.set_ylim(xy_min, xy_max)

    if show_title:
        ax.set_title(PARAMETER_NAME_MAPPING[stat_label], fontweight='semibold')

    ax.set_xlabel('simulated')
    ax.set_ylabel('estimated')
    # ax.tick_params(axis='x', which='major', labelsize=12)
    # ax.set_ylabel(stat_label)

    if ticks:
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)

    ax.text(1.0, 0.0,
            s=f'{sum(covered)}/{len(covered)}',
            # s=f'coverage={sum(covered)}/{len(covered)}',
            ha='right', va='bottom',
            fontsize=11,
            transform = ax.transAxes)


def main_plots():
    bt_stats = pd.read_csv('results/simulation_study/basictrees.tsv', sep='\t')
    ct_stats = pd.read_csv('results/simulation_study/contactrees.tsv', sep='\t')
    groundtruth = pd.read_csv('results/simulation_study/SAMPLING.log', sep='\t')

    # Drop the 0th sample
    groundtruth = groundtruth.drop(0)
    n_runs = groundtruth.shape[0]

    bt_stats['run'] = np.arange(1, 1 + n_runs)
    ct_stats['run'] = np.arange(1, 1 + n_runs)
    groundtruth['run'] = np.arange(1, 1 + n_runs)

    bt_stats.set_index('run')
    ct_stats.set_index('run')
    groundtruth.set_index('run')

    fig, axes = plt.subplot_mosaic(layout=[['height_lbl', 'bt_height', 'ct_height', '', 'err_topo', 'err_height', 'err_rate'],
                                           ['rate_lbl',   'bt_rate',   'ct_rate',   '', 'err_topo', 'err_height', 'err_rate']],
                                   figsize=(17, 6),
                                   gridspec_kw={'width_ratios': [0.03, 1, 1, 0.03, 1, 1, 1]})

    # Column + row titles
    for lbl in ['height_lbl', 'rate_lbl', '']: axes[lbl].axis('off')
    axes['height_lbl'].text(0.0, 0.5, 'root height', rotation=90, va='center', ha='center', **SUBPLOT_TITLE_ARGS)
    axes['rate_lbl'].text(0.0, 0.5, 'clock rate', rotation=90, va='center', ha='center', **SUBPLOT_TITLE_ARGS)

    axes['bt_height'].set_title('without contacTrees', pad=20, **SUBPLOT_TITLE_ARGS)
    axes['ct_height'].set_title('with contacTrees', pad=20, **SUBPLOT_TITLE_ARGS)

    settings = [
        dict(
            stat_id='treeDistance.treeDistance',
            stat_label='tree distance (RNNI)',
            stat_range=(-0.01, 7.75),
            diff_to_gt=False,
            ax=axes['err_topo'],
        ),
        dict(
            stat_id='rootHeight',
            stat_label='root height',
            stat_range=(-0.5, 0.5),
            diff_to_gt=True,
            ax=axes['err_height'],
        ),
        dict(
            stat_id='clockRate',
            stat_label='clock rate',
            stat_range=(-0.154, 0.158),
            diff_to_gt=True,
            ax=axes['err_rate'],
        ),
    ]

    for cfg in settings:
        plot_errors(bt_stats, ct_stats, groundtruth, **cfg)

    plot_coverage(bt_stats, groundtruth,
                  stat_id='rootHeight', stat_label='noCT: root height',
                  ticks=[0.5, 1, 1.5, 2],
                  lim=(0.5, 2),
                  ax=axes['bt_height'])
    plot_coverage(ct_stats, groundtruth,
                  stat_id='rootHeight', stat_label='CT: root height',
                  ticks=[0.5, 1, 1.5, 2],
                  lim=(0.5, 2),
                  ax=axes['ct_height'])

    plot_coverage(bt_stats, groundtruth,
                  stat_id='clockRate', stat_label='noCT: clock rate',
                  ticks=[0.3, 0.4, 0.5, 0.6, 0.7],
                  lim=(0.3, 0.7),
                  ax=axes['bt_rate'])
    plot_coverage(ct_stats, groundtruth,
                  stat_id='clockRate', stat_label='CT: clock rate',
                  ticks=[0.3, 0.4, 0.5, 0.6, 0.7],
                  lim=(0.3, 0.7),
                  ax=axes['ct_rate'])

    plt.tight_layout(pad=0, h_pad=2, w_pad=0.6)
    plt.show()


def main_error_plots():
    bt_stats = pd.read_csv('results/simulation_study/basictrees.tsv', sep='\t')
    ct_stats = pd.read_csv('results/simulation_study/contactrees.tsv', sep='\t')
    groundtruth = pd.read_csv('results/simulation_study/SAMPLING.log', sep='\t')

    # Drop the 0th sample
    groundtruth = groundtruth.drop(0)
    n_runs = groundtruth.shape[0]

    bt_stats['run'] = np.arange(1, 1 + n_runs)
    ct_stats['run'] = np.arange(1, 1 + n_runs)
    groundtruth['run'] = np.arange(1, 1 + n_runs)
    bt_stats.set_index('run')
    ct_stats.set_index('run')
    groundtruth.set_index('run')

    fig, axes = plt.subplot_mosaic(
        layout=[['err_topo', 'err_height', 'err_rate']],
        figsize=(7, 4))

    settings = [
        dict(
            stat_id='treeDistance.treeDistance',
            stat_label='tree distance (RNNI)',
            stat_range=(-0.02, 7.5),
            diff_to_gt=False,
            ax=axes['err_topo'],
        ),
        dict(
            stat_id='rootHeight',
            stat_label='root height',
            stat_range=(-0.5, 0.5),
            diff_to_gt=True,
            ax=axes['err_height'],
        ),
        dict(
            stat_id='clockRate',
            stat_label='clock rate',
            stat_range=(-0.152, 0.156),
            diff_to_gt=True,
            ax=axes['err_rate'],
        ),
    ]

    for cfg in settings:
        plot_errors(bt_stats, ct_stats, groundtruth, **cfg)

    plt.tight_layout(pad=0, w_pad=1)
    plt.show()


def all_coverage_plots(ct=True):
    groundtruth = pd.read_csv('results/simulation_study/SAMPLING.log', sep='\t')
    if ct:
        stats = pd.read_csv('results/simulation_study/contactrees.tsv', sep='\t')
    else:
        stats = pd.read_csv('results/simulation_study/basictrees.tsv', sep='\t')

    # Drop the 0th sample
    groundtruth = groundtruth.drop(0)
    n_runs = groundtruth.shape[0]

    stats['run'] = np.arange(1, 1 + n_runs)
    groundtruth['run'] = np.arange(1, 1 + n_runs)
    stats.set_index('run')
    groundtruth.set_index('run')

    fig, axes = plt.subplots(nrows=3 if ct else 1, ncols=4,
                             figsize=(12, 8 if ct else 2.6))
    axes = axes.flatten()

    params = [
        'clockRate',
        'rootHeight',
        'treeLength',
        'pairedTreeLength',
    ]
    if ct:
        params += [
            'conversionRate',
            'expectedConversions',
            'convCount',
            'meanConvHeight',
            'pMove',
            'moveCount',
            'movesPerConv',
        ]

    for i, param in enumerate(params):
        plot_coverage(stats, groundtruth,
                      stat_id=param,
                      stat_label=param,
                      ax=axes[i],
                      show_title=True)

    for ax in axes[i+1:]:
        ax.axis('off')

    plt.tight_layout(pad=0.1, h_pad=2, w_pad=1)
    plt.show()


if __name__ == '__main__':
    # main_plots()
    # main_error_plots()
    all_coverage_plots(ct=False)
