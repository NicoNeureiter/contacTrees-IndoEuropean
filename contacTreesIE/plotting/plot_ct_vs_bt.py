#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import string
from pathlib import Path

import numpy as np
import pandas as pd
from nexus import NexusReader

import matplotlib.pyplot as plt
import seaborn as sns

from contacTreesIE.newick_util import get_age

from contacTreesIE.plotting import plot_tree_topology, assign_node_coordinates
from contacTreesIE.plotting.plot_contact_edges import read_tree, add_grid_lines, plot_contactedges, place_contactedges_in_tree, read_contactedges

COLOR_PALETTE = [
    (0.65, 0.25, 0.85, 1.0),
    (0.0, 0.65, 0.8, 1.0),
    (1.0, 0.5, 0.05, 1.0),
]
CEDGE_COLOR = COLOR_PALETTE[2]


def plot_language_tree(nexus, ax=None, max_y=None):
    ax = ax or plt.gca()

    # Parse tree from nexus reader
    tree = read_tree(nexus)

    # Plot the tree topology
    assign_node_coordinates(tree,
                            children_sort_key=lambda c: min(c.get_leaf_names()))
    plot_tree_topology(
        node=tree,
        annotate_leafs=True,
        leaf_label_args=dict(
            rotation=60,
            horizontalalignment='right',
            x_offset=0.2,
        ),
        ax=ax
    )

    # Plot contact edges
    contactedges = read_contactedges(nexus)
    if contactedges:
        place_contactedges_in_tree(tree, contactedges)
        # plot_contactedges(contactedges, c='orange', arrow=False, ax=ax)
        plot_contactedges(contactedges, c='orange', arrow=True, lw=20, ax=ax)

    # Adjust plot layout
    add_grid_lines(ax, max_y=(max_y or get_age(tree)))
    ax.set_ylim(None, max_y)
    ax.axis('off')


def add_subplot_labels(axes, x_relative=-0.05, y_relative=1.0):
    for ax, label in zip(axes.flatten(), string.ascii_lowercase):
        bbox = ax.get_tightbbox(fig.canvas.get_renderer())
        x = (1-x_relative)*bbox.x0 + x_relative*bbox.x1
        y = (1 - y_relative) * bbox.y0 + y_relative * bbox.y1
        fig.text(
            x, y, label,
            fontsize=18,
            fontweight="bold",
            va="top",
            ha="left",
            transform=None
        )


def declutter_violin_plot(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_yaxis().set_ticks([])
    ax.tick_params(axis=u'x', which=u'both', length=0)
    ax.set_xlabel('')
    ax.set_ylabel('')


def format_scientific(x):
    if x == 0:
        return '0'

    b = int(np.floor(np.log10(x)))
    a = int(x / (10**b))
    return f'${a} · 10^{{{b}}}$'


if __name__ == '__main__':
    BASE_DIR = Path('results/fix_clock_stdev')

    subplot_label_args = dict(fontsize=20, fontweight="semibold", va="top", ha="left")
    subplot_title_args = dict(fontsize=20, pad=30)

    # Load the summary tree files with NexusReader
    ct_summary_nexus = NexusReader.from_file(BASE_DIR / 'CT_fixTopo' / 'covarion' / 'CT_fixTopo_covarion.summary.tree')
    bt_summary_nexus = NexusReader.from_file(BASE_DIR / 'BT_fixTopo' / 'covarion' / 'BT_fixTopo_covarion.summary.tree')

    # Load the log files into a pandas data-frame
    ct = pd.read_csv(BASE_DIR / 'CT_fixTopo' / 'covarion' / 'CT_fixTopo_covarion.log', sep='\t')
    bt = pd.read_csv(BASE_DIR / 'BT_fixTopo' / 'covarion' / 'BT_fixTopo_covarion.log', sep='\t')
    bt_nL = pd.read_csv(BASE_DIR / 'BT_fixTopo' / 'covarion_noLoans' / 'BT_fixTopo_covarion_noLoans.log', sep='\t')
    # ct['mode'] = 'with contacTrees'
    # bt['mode'] = 'w/o contactTrees'
    # bt_nL['mode'] = 'w/o contacTrees\nloanwords removed'
    ct['mode'] = 'CT'
    bt['mode'] = 'noCT'
    bt_nL['mode'] = 'noCT-filtered'
    df = pd.concat([bt, bt_nL, ct])


    # # Before plotting anything we do a quick calculation of the marginal likelhoods and the Bayes factor for BT vs CT
    # ESS = 300
    #
    # MU_CEDGES = 0.25
    # print(ct.columns)
    # print(ct['convCount'])
    # print((ct.convCount == 0))
    #
    # import scipy.stats as stats
    # p0_post = (ct.convCount == 0).mean()
    # p0_prio = stats.poisson(MU_CEDGES).pmf(0)
    #
    # print(p0_post)
    # print(p0_prio)
    #
    # c0_post = p0_post*ESS + 1
    # c1_post = (1-p0_post)*ESS
    # c0_prio = p0_prio*ESS
    # c1_prio = (1-p0_prio)*ESS
    # print( f'{c1_post} * {c0_prio} / ({c1_prio} * {c0_post}) = {c1_post * c0_prio / (c0_post * c1_prio)}')
    #
    # sns.histplot(data=ct, x='convCount', stat='probability')
    #
    # posterior_mean = ct.convCount.mean()
    # posterior_variance = ct.convCount.var()
    # gamma_scale = posterior_variance / posterior_mean
    # gamma_shape = posterior_mean / gamma_scale
    # dgamma = stats.gamma(a=gamma_shape, scale=gamma_scale)
    # print(f'P[X<10] = {dgamma.cdf(10)}')
    # print(f'P[X<20] = {dgamma.cdf(20)}')
    #
    # x = np.linspace(30, 50, 1000)
    # plt.plot(x, dgamma.pdf(x))
    # plt.show()
    #
    # exit()

    # Set up the plot

    # fig, axes = plt.subplots(1, 4, figsize=(22, 8), gridspec_kw={'width_ratios':[3,3,2,2]})
    fig, axes = plt.subplot_mosaic(layout=[['tree', 'height'],
                                           ['tree', 'rate']],
                                   figsize=(18, 10),
                                   gridspec_kw={'width_ratios': [2.5, 1]})
    # gs = gridspec.GridSpec(1, 4, width_ratios=[2,2,1,1], figure=fig)
    # axes = [plt.subplot(gs[i]) for i in range(4)]

    # # Plot the BT tree
    # ax = axes[0]
    # ax.set_title('summary tree w/o contact', **subplot_title_args)
    # plot_language_tree(bt_summary_nexus, ax=ax, max_y=7300)

    # plot the CT tree
    ax = axes['tree']
    ax.set_title('CT: summary tree with contact edges', **subplot_title_args)
    # plot_language_tree(ct_summary_nexus, ax=ax, max_y=7300)
    plot_language_tree(ct_summary_nexus, ax=ax, max_y=6050)

    # Plot the tree height posterior distribution
    ax = axes['height']
    ax.set_title('tree height', **subplot_title_args)
    sns.violinplot(data=df, x='mode', y='rootHeight', inner="quartile",
                   palette=COLOR_PALETTE, ax=ax).tick_params(labelsize=12)
    ax.set_xlim(-0.7, 2.5)
    ax.set_ylim(-750, 13000)
    add_grid_lines(ax, min_y=0, max_y=13000, dy=2000, zorder=0)
    declutter_violin_plot(ax)

    # Add hlines for comparison with previous studies
    ref_color = 'navy'
    ax.axhline(4200, color=ref_color, zorder=0)
    ax.text(2.5, 4200, 'Chang 2015', color=ref_color,
            fontsize=11, va="bottom", ha="left", rotation=60)
    ax.axhline(5700, color=ref_color, zorder=0)
    ax.text(2.5, 5700, 'Bouckaert 2012', color=ref_color,
            fontsize=11, va="bottom", ha="left", rotation=60,)


    # Plot the clock rate posterior distribution
    ax = axes['rate']
    # ax = axes[1, 1]
    ax.set_title('clock rate', **subplot_title_args)
    sns.violinplot(data=df, x='mode', y='clockRate.medium', inner="quartile",
                   palette=COLOR_PALETTE, ax=ax).tick_params(labelsize=12)
    ax.set_xlabel('')
    ax.set_xlim(-0.7, 2.5)
    ax.set_ylim(-0.000005, 0.00009)
    add_grid_lines(ax, min_y=0.0000, max_y=0.0001, dy=0.00002, zorder=0, format=format_scientific)
    declutter_violin_plot(ax)

    # Add sub-plot labels
    # fig.text(0.0, 1., 'a', **subplot_label_args)
    # fig.text(0.3, 1., 'b', **subplot_label_args)
    # fig.text(0.604, 1., 'c', **subplot_label_args)
    # fig.text(0.816, 1., 'd', **subplot_label_args)
    fig.text(0.005, 1, 'a) ', **subplot_label_args)
    fig.text(0.7, 1, 'b)', **subplot_label_args)
    fig.text(0.7, 0.5, 'c)', **subplot_label_args)

    fig.tight_layout(pad=0.2, h_pad=4, w_pad=4)
    plt.show()
