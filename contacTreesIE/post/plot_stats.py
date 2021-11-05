import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns

def plot_simulation_results():
    fig, ax = plt.subplots(figsize=(10, 6))

    bt_stats = pd.read_csv('/home/nico/workspace/contactrees_simulation_study/contactrees_vs_basictrees/remote/2021_09_26/basictrees/results.tsv', sep='\t')
    ct_stats = pd.read_csv('/home/nico/workspace/contactrees_simulation_study/contactrees_vs_basictrees/remote/2021_09_26/contactrees/results.tsv', sep='\t')
    groundtruth = pd.read_csv('/home/nico/workspace/contactrees_simulation_study/contactrees_vs_basictrees/remote/2021_09_26/sampling/SAMPLING.log', sep='\t')
    groundtruth = groundtruth.drop(0)

    print(bt_stats.shape)
    print(ct_stats.shape)
    print(groundtruth.shape)

    stat_id = 'treeDistance.treeDistance'
    stat_label = 'tree distance (RNNI)'
    stat_range = (0, 6)
    diff_to_gt = False

    # # stat_id = 'clockRate'
    # stat_id = 'rootHeight'
    # stat_label = 'clock rate'
    # stat_range = (0, 4)
    # diff_to_gt = True

    bt_mean = bt_stats[stat_id + '.mean'].astype(float).to_numpy()
    bt_lo = bt_stats[stat_id + '.95%HPDlo'].astype(float).to_numpy()
    bt_up = bt_stats[stat_id + '.95%HPDup'].astype(float).to_numpy()
    ct_mean = ct_stats[stat_id + '.mean'].astype(float).to_numpy()
    ct_lo = ct_stats[stat_id + '.95%HPDlo'].astype(float).to_numpy()
    ct_up = ct_stats[stat_id + '.95%HPDup'].astype(float).to_numpy()
    # gt = groundtruth[stat_id.partition('.')[0]].to_numpy()

    # covered_colors = np.array(['red', 'green'])
    # bt_covered = covered_colors[np.logical_and(bt_lo <= gt, gt <= bt_up).astype(int)]
    # ct_covered = covered_colors[np.logical_and(ct_lo <= gt, gt <= ct_up).astype(int)]
    #
    # for i in range(len(ct_mean)):
    #     # ax.plot([gt[i], gt[i]], [ct_lo[i], ct_up[i]], c=ct_covered[i])
    #     ax.plot([gt[i], gt[i]], [bt_lo[i], bt_up[i]], c=bt_covered[i])
    # ax.plot([0.0, 1], [0.0, 1], c='lightgrey', zorder=0)

    # if diff_to_gt:
    #     bt_stat -= gt
    #     ct_stat -= gt

    sns.kdeplot(ct_mean, shade=True, label='With contacTrees')
    sns.kdeplot(bt_mean, shade=True, label='Without contacTrees')

    plt.legend(prop={'size': 14})

    # plt.xlim(0.0, 1)
    # plt.ylim(0.0, 1)
    plt.xlim(*stat_range)
    plt.xlabel(stat_label, fontdict={'size': 14})
    plt.xticks(range(stat_range[0], stat_range[1] + 1))
    ax.tick_params(axis='x', which='major', labelsize=12)
    plt.ylabel('')
    plt.yticks([])

    # plt.axvline(0.0)


if __name__ == '__main__':
    # plot_simulation_results()


    bt = pd.read_csv('results/workstation/BT_fixTopo/covarion_clockStdv01/BT_fixTopo_covarion_clockStdv01_1.log', sep='\t')
    bt_noloan = pd.read_csv('results/workstation/BT_fixTopo/covarion_noLoans_clockStdv01/BT_fixTopo_covarion_noLoans_clockStdv01_1.log', sep='\t')
    ct = pd.read_csv('results/workstation/CT_fixTopo/covarion/CT_fixTopo_covarion_1.log', sep='\t')

    ct['mode'] = 'with contacTrees'
    bt['mode'] = 'without contacTrees'
    # bt_noloan['mode'] = 'without contacTrees\nloan words removed'

    df = pd.concat([ct, bt, bt_noloan])

    # stat = 'rootHeight'
    stat = 'clockRate.medium'
    # stat = 'clockStdev'
    # stat = 'freqParameter.2'

    # sns.set(rc={'figure.figsize': (8, 6)})

    rcParams['figure.figsize'] = 9, 7
    sns.set_style('whitegrid')

    sns.violinplot(data=df, x='mode', y=stat, inner="quartile")

    # palette = sns.cubehelix_palette(4, rot=-0.2, light=.6, hue=0.25)
    # grid = sns.FacetGrid(df, row='mode', hue='mode', aspect=5, height=2.5,
    #                      palette=palette, sharex=True, sharey=True)
    #
    # grid.map(sns.kdeplot, stat, clip_on=True, shade=True, alpha=1, lw=.5, bw='scott')
    # grid.map(sns.kdeplot, stat, clip_on=True, color='w', lw=1, bw='scott')
    #
    # grid.map(plt.axhline, y=0, lw=1, clip_on=False)
    #
    # grid.fig.subplots_adjust(hspace=0.02, left=0.075, bottom=0.07, top=0.93)
    #
    #
    # grid.set(yticks=[])
    # grid.set(xticks=[0, 2500, 5000, 7500, 10000])
    # grid.set(xlim=(0, 12000))
    # # grid.set_xlabels('')
    # grid.set_titles('')


    plt.xlabel('')
    # plt.ylabel('tree height', fontdict={'size': 14})
    plt.ylabel('clock rate', fontdict={'size': 14})
    # plt.ylabel('frequency of presence', fontdict={'size': 14})

    # plt.ylim(0, 13000)
    # plt.ylim(0, 0.0001)

    ax = plt.gca()
    ax.tick_params(axis='x', labelsize=14)

    plt.tight_layout()
    plt.show()

