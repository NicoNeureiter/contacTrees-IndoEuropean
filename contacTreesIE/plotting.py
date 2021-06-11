#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import colorsys
from typing import Callable

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from newick import Node

from contacTreesIE.newick_util import get_age

PI = np.pi
TAU = 2 * np.pi


def get_location(node):
    return node.location


def plot_edge(parent, child, ax=None, **kwargs_arrow):
    if ax is None:
        ax = plt.gca()
    color = kwargs_arrow.pop('color', 'k')
    alpha = kwargs_arrow.pop('alpha', 1.)
    return ax.plot([parent.x, child.x], [parent.y, child.y], c=color, alpha=alpha,
                   solid_capstyle='round', **kwargs_arrow)


def plot_tree(tree: Node, color='k', alpha=1.0, ax=None, lw=None, **plot_args):
    for node in tree.walk():
        if node.ancestor is None:
            continue

        plot_edge(node.ancestor, node, color=color, alpha=alpha, ax=ax, lw=lw, **plot_args)


def n_leafs(tree: Node, n=0):
    if tree.is_leaf:
        n = 1
    else:
        for c in tree.descendants:
            n += n_leafs(c)

    return n


def assign_node_coordinates(
        node: Node,
        left: float = 0,
        children_sort_key: Callable = None,
):
    """

    Args:
        node (Node):
        left:
        children_sort_key:
        y_parent (float): the vertical coordinate of the parent node.

    """
    if children_sort_key is None:
        children_sort_key = lambda c: len(c.get_leaves())

    node.y = get_age(node)

    if node.is_leaf:
        node.x = left

    else:
        children_sorted = sorted(node.descendants, key=children_sort_key)
        # children_sorted = node.descendants

        x_children = []
        for i, c in enumerate(children_sorted):
            assign_node_coordinates(
                node=c,
                left=left,
                children_sort_key=children_sort_key,
            )
            x_children.append(c. x)
            left += n_leafs(c)

        node.x = np.mean(x_children)


def plot_tree_topology(
        node: Node,
        node_plotter: Callable = None,
        annotate_leafs: bool = False,
        ax: plt.Axes = None,
        **plot_kwargs
):
    """

    Args:
        node (Node):
        node_plotter:
        annotate_leafs (bool): Whether or not to show tip node-labels
        ax:
        **plot_kwargs:

    Returns:
        float: the x coordinate of tree in the plot.
        float: the y coordinate of tree in the plot.
    """
    if ax is None:
        ax = plt.gca()

    if node_plotter is not None:
        node_plotter(node, node.x, node.y, ax=ax)

    if not ('c' in plot_kwargs or 'color' in plot_kwargs):
        plot_kwargs['color'] = 'k'

    if node.descendants:
        for c in node.descendants:
            ax.plot([node.x, c.x, c.x], [node.y, node.y, c.y], **plot_kwargs)
            plot_tree_topology(
                node=c,
                node_plotter=node_plotter,
                annotate_leafs=annotate_leafs,
                ax=ax,
                **plot_kwargs
            )

    if annotate_leafs and node.is_leaf:
        plt.text(node.x, node.y - 0.015, node.name,
                 horizontalalignment='center', verticalalignment='top',
                 fontsize=10, rotation=90)


def plot_network_topology(network, **plot_tree_args):
    tree = network.tree
    get_name = lambda n: n.name
    plot_tree_topology(tree, annotate_leafs=True, children_sort_key=get_name, **plot_tree_args)
    for e in network.contact_edges:
        b_start = e.branch_start
        b_end = e.branch_end

        plt.arrow(b_start.x, e.height, b_end.x-b_start.x, 0, width=0.001, head_length=0.05,
                  color='r', length_includes_head=True, zorder=10)


def colorline(ax, x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0), linewidth=3, alpha=1.0, **kwargs):
    """
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    from: https://nbviewer.jupyter.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    """

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])

    z = np.asarray(z)

    def make_segments(x, y):
        """
        Create list of line segments from x and y coordinates, in the correct format for LineCollection:
        an array of the form   numlines x (points per line) x 2 (x and y) array
        """

        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        return segments

    segments = make_segments(x, y)
    lc = LineCollection(segments, colors=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=1, zorder=1, **kwargs)

    # ax = plt.gca()
    ax.add_collection(lc)

    return lc


def color_tree(tree, hue=0.5, saturation=0.8, value=0.8, rate=1.):
    tree.color = np.array(colorsys.hsv_to_rgb(hue, saturation, value))
    sign = np.random.permutation([-1, 1])

    for i, child in enumerate(tree.descendants):
        child_hue = np.random.normal(hue, 2*rate*child.length) % 1
        # s = (i-0.5) * 2
        # child_hue = (hue + sign[i]*rate*child.length) % 1
        child_sat = np.clip(np.random.normal(saturation, rate*child.length), 0.6, 0.9)
        # child_sat = np.random.uniform(0.7, 0.8)
        child_value = np.clip(np.random.normal(value, rate*child.length), 0.5, 0.9)
        # child_value = np.random.uniform(0.6, 0.9)
        color_tree(child, hue=child_hue, saturation=child_sat, value=child_value, rate=rate*0.8)
