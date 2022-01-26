#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import colorsys
from typing import Callable, Optional

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from newick import Node

from contacTreesIE.newick_util import get_age, get_root, get_height

# Custom colors
GRAY_0 = (0.95, 0.95, 0.95)
GRAY_1 = (0.85, 0.85, 0.85)
GRAY_75 = (0.75, 0.75, 0.75)
GRAY_2 = (0.65, 0.65, 0.65)
GRAY_3 = (0.5, 0.5, 0.5)
CEDGE_COLOR = (0.85, 0.65, 0)
DONOR_COLOR = (0., 0.6, 0.55)
RECEIVER_COLOR = (0.8, 0.3, 0)


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
        children_sort_key: Optional[Callable] = None,
        flip_order_at_nodes = (),
        clade_label_getter: Optional[Callable] = None
):
    """

    Args:
        node: The root node of the tree for which the node coordinates are assigned.
        left: The x-coordinate of the left-most node.
        children_sort_key: The key-function for sorting the tree topology (which clade goes left or right at a split).
        flip_order_at_nodes: Names of the nodes where the order defined by children_sort_key should be reversed.
    """
    if children_sort_key is None:
        children_sort_key = lambda c: len(c.get_leaves())
    if clade_label_getter is None:
        clade_label_getter = lambda _: None

    node.y = get_age(node)

    if node.is_leaf:
        node.x = left
        return left

    else:
        children_sorted = (sorted(node.descendants, key=children_sort_key))
        if node.name in flip_order_at_nodes:
            children_sorted = children_sorted[::-1]

        x_children = []

        previous_clade = None
        for i, c in enumerate(children_sorted):
            clade = clade_label_getter(c)
            if (i > 0) and (clade != previous_clade):
                    left += 0.5

            previous_clade = clade

            left = assign_node_coordinates(
                node=c,
                left=left,
                children_sort_key=children_sort_key,
                flip_order_at_nodes=flip_order_at_nodes,
                clade_label_getter=clade_label_getter,
            )
            x_children.append(c. x)

            if i < len(children_sorted) - 1:
                left += 1

        node.x = np.mean(x_children)

        return left

def plot_tree_topology(
        node: Node,
        node_plotter: Callable = None,
        annotate_leafs: bool = False,
        annotate_internal_nodes: bool = False,
        leaf_label_args: dict = None,
        ax: plt.Axes = None,
        clade_color_getter: Optional[Callable] = None,
        **plot_kwargs
):
    """

    Args:
        node (Node): The root node of the tree to be plotted
        node_plotter (function): A function handler defining how a node is plotted.
        annotate_leafs (bool): Whether or not to show tip node labels
        annotate_internal_nodes (bool): Whether or not to show internal node labels
        ax: The Matplotlib axis used for plotting the tree topology.
        **plot_kwargs: Plotting parameters.

    Returns:
        float: the x coordinate of tree in the plot.
        float: the y coordinate of tree in the plot.
    """
    leaf_label_args = leaf_label_args or {}
    ax = ax or plt.gca()

    # if not ('c' in plot_kwargs or 'color' in plot_kwargs):
    color = plot_kwargs.pop('color', 'k')
    color = plot_kwargs.pop('c', color)
    if clade_color_getter is None:
        clade_color_getter = lambda _: color

    if node_plotter is not None:
        node_plotter(node, node.x, node.y, ax=ax)

    if node.descendants:
        for c in node.descendants:
            c: Node
            yoffset = 0 if c.is_leaf else 15
            ax.plot([node.x, c.x, c.x], [node.y, node.y, c.y + yoffset],
                    c=clade_color_getter(node),
                    **plot_kwargs)
            plot_tree_topology(
                node=c,
                node_plotter=node_plotter,
                annotate_leafs=annotate_leafs,
                annotate_internal_nodes=annotate_internal_nodes,
                leaf_label_args=leaf_label_args,
                clade_color_getter=clade_color_getter,
                ax=ax,
                **plot_kwargs
            )

    label_offset = 0.01 * get_height(get_root(node))
    if annotate_leafs and node.is_leaf:
        leaf_label_args.setdefault('fontsize', 10)
        leaf_label_args.setdefault('horizontalalignment', 'center')
        leaf_label_args.setdefault('verticalalignment', 'top')

        leaf_labels_at_0 = True
        if leaf_labels_at_0:
            y = 0
            if node.y > 0.000001:
                plt.plot([node.x, node.x], [0, node.y], color='lightgray', ls='dotted', zorder=0)
                print('...', node.name)
        else:
            y = node.y

        x_offset = leaf_label_args.pop('x_offset', 0.0)
        ax.text(node.x + x_offset, y - label_offset, node.name, **leaf_label_args)
        leaf_label_args['x_offset'] = x_offset

    if annotate_internal_nodes and not node.is_leaf:
        ax.text(node.x, node.y - label_offset, node.name,
                 horizontalalignment='center', verticalalignment='top',
                 fontsize=20)



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
