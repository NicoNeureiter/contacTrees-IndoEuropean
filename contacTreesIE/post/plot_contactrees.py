#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from pathlib import Path
from typing import List, Dict
from operator import attrgetter
from collections import defaultdict
from enum import IntEnum

import numpy as np
import pandas as pd
from newick import Node
from nexus import NexusReader

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from contacTreesIE.newick_util import get_age, ContactEdge, get_actual_leaves, collect_contactedges, drop_contactedges, \
    translate_node_names, remove_dead_end

GRAY_0 = (0.95, 0.95, 0.95)
GRAY_1 = (0.85, 0.85, 0.85)
GRAY_2 = (0.65, 0.65, 0.65)
GRAY_3 = (0.5, 0.5, 0.5)
CEDGE_COLOR = (0.85, 0.65, 0)
DONOR_COLOR = (0., 0.6, 0.55)
RECEIVER_COLOR = (0.8, 0.3, 0)

class ZOrder(IntEnum):

    GRID_LINE = 1
    CONTACT_EDGE = 2
    LANGUAGE_TREE_FRAME = 3
    CONTACT_EDGE_STUB = 4
    LANGUAGE_TREE = 5
    WORD_TREE = 6
    CONTACT_EDGE_HIGHLIGHTED = 10

# TODO plot each summary tree with transparent contact edge samples (densitree style)
# TODO same for single word trees to visualize
# TODO [not plotting] simply summarize where each loan word came from in each language


# class ContacTree(object):
#
#     """Class encapsulating python-newick trees and adds conversion edges to represent an
#     contacTree.
#
#     Attributes:
#         root (Node): ...
#         contact_edges (List[ContactEdge]): ...
#     """
#
#     def __init__(self,
#                  node: Node,
#                  contactedges: dict = None):
#         self.root = node
#         if contactedges is None:
#             self.contactedges = {}
#         else:
#             self.contactedges = contactedges


def read_contactedges(nexus: NexusReader, burnin=0.1, return_n_samples=False) -> List[ContactEdge]:
    translate = nexus.trees.translators
    tree_handler = nexus.trees
    n_burnin_samples = int(burnin * tree_handler.ntrees)
    trees: List[Node] = [t.newick_tree for t in tree_handler.trees[n_burnin_samples::4]]
    n_samples = len(trees)
    print(n_samples)

    contactedges = []
    for tree in trees:
        translate_node_names(tree, translate)
        contactedges += collect_contactedges(tree)

    if return_n_samples:
        return contactedges, n_samples
    else:
        return contactedges


def read_blockset(nexus: NexusReader):
    assert 'contactrees' in nexus.blocks
    return nexus.contactrees.block[1].strip(' \t;.,').split()[1:]


def read_tree(nexus: NexusReader) -> Node:
    translate = nexus.trees.translators
    trees = [t.newick_tree for t in nexus.trees.trees]

    assert len(trees) == 1
    tree = trees[0]
    translate_node_names(tree, translate)

    drop_contactedges(tree)

    return tree


def place_contactedges_in_tree(
        tree: Node,
        contactedges: List[ContactEdge],
        filter_height: bool = False
):
    leaves = sorted(get_actual_leaves(tree), key=attrgetter('name'))

    def compress_clade(node):
        clade = get_actual_leaves(node)
        return ''.join([str(int(l in clade)) for l in leaves])

    clade_map = {
        compress_clade(node): node for node in tree.walk()
    }

    for cedge in contactedges.copy():
        valid_donor = cedge.donor_clade in clade_map
        valid_receiver = cedge.receiver_clade in clade_map
        if not (valid_donor and valid_receiver):
            contactedges.remove(cedge)
            continue

        cedge.donor_node = clade_map[cedge.donor_clade]
        cedge.receiver_node = clade_map[cedge.receiver_clade]

        # if filter_height:
        #     if not (cedge.donor_node.height < cedge.height < cedge.donor_node.ancestor.height):
        #         contactedges.remove(cedge)
        #     elif not (cedge.receiver_node.height < cedge.height < cedge.receiver_node.ancestor.height):
        #         contactedges.remove(cedge)


def plot_contactedges(contactedges, word=None, arrow=False, **plot_kwargs):
    zorder = plot_kwargs.pop('zorder', ZOrder.CONTACT_EDGE)
    plot_kwargs['color'] = plot_kwargs.pop('c', 'k')
    width = plot_kwargs.pop('lw', 1.0)
    head_width = 3. * width

    if word is not None:
        contactedges = filter(lambda ce: word in ce.affected_blocks, contactedges)

    for cedge in contactedges:
        x1 = cedge.donor_node.x
        x2 = cedge.receiver_node.x
        y = cedge.height
        direction = np.sign(x2 - x1)

        if arrow:
            plt.arrow(x=x1, y=y, dx=(x2-x1), dy=0, width=width,
                      length_includes_head=True, head_width=head_width, head_length=0.1,
                      ec=(0., 0., 0.), lw=0.1,
                      zorder=ZOrder.CONTACT_EDGE_STUB, **plot_kwargs)

        else:
            plt.plot([x1, x2], [y, y], lw=width, zorder=zorder, **plot_kwargs)
            plt.plot([x2 - direction*0.5, x2 - direction*0.05], [y, y], lw=width, zorder=ZOrder.CONTACT_EDGE_STUB, **plot_kwargs)


def plot_contactedge_arrow(cedge, **plot_kwargs):
    zorder = plot_kwargs.pop('zorder', ZOrder.CONTACT_EDGE_HIGHLIGHTED)
    plot_kwargs['color'] = plot_kwargs.pop('c', 'k')
    width = plot_kwargs.pop('lw', 1.0)
    head_width = 3. * width

    x1 = cedge.donor_node.x
    x2 = cedge.receiver_node.x
    y = cedge.height
    direction = np.sign(x2 - x1)

    plt.arrow(x=x1 + 0.3*direction, y=y, dx=(x2-x1) - 0.6*direction, dy=0, width=width,
              length_includes_head=True, head_width=head_width, head_length=0.2,
              ec=(1., 1., 1.), lw=0.1,
              zorder=zorder, **plot_kwargs)

def plot_densiedges(trees_nexus, summarytree_nexus):
    contactedges = read_contactedges(nexus=trees_nexus)
    tree = read_tree(nexus=summarytree_nexus)
    place_contactedges_in_tree(tree, contactedges)

    assign_node_coordinates(tree)
    plot_tree_topology(tree, annotate_leafs=True)
    plot_contactedges(contactedges, word='lake', alpha=0.2)


def plot_network(nexus, annotate_leafs=True):
    contactedges = read_contactedges(nexus)
    tree = read_tree(nexus)
    place_contactedges_in_tree(tree, contactedges)

    assign_node_coordinates(tree)
    # plot_contactedges(contactedges, lw=0.006*get_height(tree), c=GRAY_1)
    plot_contactedges(contactedges, lw=3, c=GRAY_2)
    plot_tree_topology(tree, annotate_leafs=annotate_leafs, lw=7.5, color='white', zorder=ZOrder.LANGUAGE_TREE_FRAME)
    plot_tree_topology(tree, annotate_leafs=annotate_leafs, lw=6, color=GRAY_2, zorder=ZOrder.LANGUAGE_TREE)

    # Plot internal nodes
    x = [node.x for node in tree.walk() if not node.is_leaf] + [c.donor_node.x for c in contactedges]
    y = [node.y for node in tree.walk() if not node.is_leaf] + [c.height for c in contactedges]
    plt.scatter(x, y, s=100, color=GRAY_2, zorder=ZOrder.LANGUAGE_TREE)

    add_grid_lines(plt.gca(), max_y=get_age(tree))


def plot_wordtree(nexus, word):
    contactedges = read_contactedges(nexus)
    tree = read_tree(nexus)
    place_contactedges_in_tree(tree, contactedges)
    assign_node_coordinates(tree)
    for node in tree.walk():
        node.overshadowed = False
        node.original_ancestor = node.ancestor

    active_nodes = {node: node for node in tree.walk()}

    # events = sorted(contactedges + list(tree.walk()), key=get_height)
    events = sorted(contactedges, key=attrgetter('height'))
    for i, cedge in enumerate(events):
        donor: Node = active_nodes[cedge.donor_node]
        receiver: Node = active_nodes[cedge.receiver_node]
        new_parent: Node = donor.original_ancestor
        old_parent: Node = receiver.ancestor

        if word not in cedge.affected_blocks:
            continue
        if receiver.overshadowed:
            continue

        # 1. Remove reveiver_node from it's parent
        old_parent.descendants.remove(receiver)

        # 2. create a new node at the donor
        contactnode_length = get_age(new_parent) - cedge.height
        receiver.length = cedge.height - get_age(receiver)
        contactnode = Node(
            name=f'Contact_{i}',
            length=str(contactnode_length),
        )
        new_parent.add_descendant(contactnode)
        contactnode.add_descendant(receiver)
        if not donor.overshadowed:
            contactnode.add_descendant(donor)
            new_parent.descendants.remove(donor)
            donor.length = cedge.height - get_age(donor)

        contactnode.x = donor.x
        contactnode.y = get_age(contactnode)
        contactnode.original_ancestor = contactnode.ancestor
        contactnode.overshadowed = False

        receiver.overshadowed = True
        donor.overshadowed = False
        active_nodes[donor] = contactnode

        # If the old donor_ancestor is a dead-end, remove it
        remove_dead_end(old_parent)

    # plot_tree_topology(tree, color=GRAY_2, lw=1.5)
    plot_tree_topology(tree, color='black', lw=1, zorder=ZOrder.WORD_TREE)

    # Plot internal nodes
    internal_nodes = [node for node in tree.walk() if len(node.descendants) > 1]
    x = [node.x for node in internal_nodes]
    y = [node.y for node in internal_nodes]
    plt.scatter(x, y, s=30, color='k', zorder=100)


def clean_lexeme(lexeme):
    return lexeme.replace('  ', ' ').replace('  ', ' ').lower()


def plot_word_labels(
        nexus: NexusReader,
        words=None,
        word_support: Dict[str, float] = None,
        donor_clade: List[str] = None,
        receiver_clade: List[str] = None,
        header_interval: int = 20,
        no_tree: bool = False,
        ax: plt.Axes = None
):
    if ax is None:
        ax = plt.gca()

    tree = read_tree(nexus)
    drop_contactedges(tree)
    assign_node_coordinates(tree)

    # Load ielex dataset for lexemes
    df = pd.read_csv('resources/data-mittellatein-2021-09-30.csv', sep='\t', dtype=str)
    df = df.loc[~df.status.isin(['EXCLUDE', 'LOAN,EXCLUDE', 'WRONG'])]
    df['concept'] = df['cc_alias'].apply(lambda s: s.split('-')[0])

    for i, row in df.copy().iterrows():
        if row.language == 'Latin_M':
            row = row.copy()
            row.language = 'Latin_preserved'
            df = df.append(row, ignore_index=True)

    if words is None:
        words = np.unique(df.concept.to_numpy())

    if word_support:
        words = sorted(words, key=word_support.__getitem__, reverse=True)

    def get_word_in_language(language, concept):
        x = df.loc[(df.language == language) & (df.concept == concept)]
        if len(x) == 0:
            return '?'

        # Ensure unique entries per cognate class
        cc_done = set()
        lexemes = []
        for _, row in x.iterrows():
            if row.cc_id in cc_done:
                continue
            lexemes.append(str(row.lexeme))
            cc_done.add(row.cc_id)

        # assert len(x.lexeme) == len(set(x.lexeme.to_list()))
        # I wanted to assert that lexemes are unique, but this seems to be intended:
        # E.g. `karlmaðr` in Old Norse is in two cognate classes with the same meaning,
        # since `karl-` (male) and `-maðr` (man) have separate etymological origins.

        return clean_lexeme('\n'.join(lexemes))

    colors = plt.get_cmap('tab10')
    cognate_colors_by_concept = defaultdict(dict)

    def get_cognate_colors(language, concept):
        cognate_colors = cognate_colors_by_concept[concept]
        x = df.loc[(df.language == language) & (df.concept == concept)]
        x_colors = []

        # assert len(x.cc_id) == len(set(x.cc_id.to_list())), x
        # I wanted to assert that cognate IDs are unique (in a meaning class), but this
        # seems to be intended: E.g. `dʊir` and `dʒɔur` in Faroese have the same
        # meaning (animal) and etymological origin. This does not affect the analysis in
        # any way, but is probably still good to list as separate items in the data.

        for cc_id in x.cc_id:
            if cc_id not in cognate_colors:
                i_color = len(cognate_colors)
                cognate_colors[cc_id] = colors(i_color)

            c = cognate_colors[cc_id]
            if c not in x_colors:
                x_colors.append(c)

        return x_colors

    text_args = dict(
        horizontalalignment='center',
        verticalalignment='top',
        fontsize=7.5,
    )

    Y_TOP = -100.0
    DY = 300.0

    # Plot the cognate classes in color coding + lexeme for each languages and meaning
    for node in tree.get_leaves():
        if not no_tree:
            # Plot a visual connection from label to leaf, if leaf isn't contemporary
            if get_age(node) > 0.001:
                plt.plot(
                    [node.x, node.x], [0, node.y],
                    color='lightgray', ls='dotted',
                    zorder=0
                )

        y = Y_TOP + 0.5*DY
        for i, word in enumerate(words):
            if i % header_interval == 0:
                y -= 0.5*DY
                plt.text(x=node.x, y=y, s=node.name, weight='bold', **text_args)
                y -= DY

            plt.text(
                x=node.x,
                y=y,
                s=get_word_in_language(node.name, word),
                **text_args
            )
            c = get_cognate_colors(node.name, word)
            plt.scatter(
                x=node.x + np.linspace(-0.3, 0.3, len(c)+2)[1:-1],
                y=[y + 60] * len(c),
                c=c,
                s=20,
                marker='s'
            )

            y -= DY

        y_bot = y

        # Plot frame around donor + receiver clades
        is_donor = donor_clade and node.name in donor_clade
        is_receiver = receiver_clade and node.name in receiver_clade
        if is_donor or is_receiver:
            ax.add_patch(
                Rectangle(
                    xy=(node.x - 0.5, y_bot - 0.5*DY),
                    width=1,
                    height=Y_TOP - y_bot,
                    color=GRAY_0,
                    # alpha=0.5,
                    zorder=0,
                )
            )

    # Plot the meaning class labels
    y = Y_TOP + 0.5*DY
    for i, word in enumerate(words):
        if i % header_interval == 0:
            y -= 1.5*DY

        plt.text(
            x=-1.2,
            y=y,
            s=word,
            weight='bold',
            **text_args
        )

        if word_support:
            plt.text(
                x=-1.2,
                y=y - (0.3 * DY),
                s=f'{word_support[word]:.2f}',
                color='grey',
                **text_args
            )

        y -= DY



def plot_wordtrees(nexus: NexusReader, word_trees_directory=Path('./wordtrees/')):
    """Plot all word trees with the corresponding words and store the plots as PDF files
    in `word_trees_directory`."""

    # Create the word-trees directory
    os.makedirs(word_trees_directory, exist_ok=True)

    # Read the contact edges and attach them to nodes in the `nexus` tree.
    contactedges = read_contactedges(nexus)
    place_contactedges_in_tree(read_tree(nexus), contactedges)

    words_grouped = defaultdict(list)
    for w in read_blockset(nexus):
        g = ''.join(str(int(w in cedge.affected_blocks)) for cedge in contactedges)
        words_grouped[g].append(w)

    for g, words in words_grouped.items():
        print()
        fig_width = 40
        fig_height = 16 + (0.25 * len(words))
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        # if int(g) == 0:
        #     continue

        print('Plotting marginal tree', g)
        print(f'\tContact edges ({len(contactedges)}):')
        for i, cedge in enumerate(contactedges):
            if g[i] == '1':
                print('\t\t', ','.join(cedge.donor_node.get_leaf_names()),
                      '->', ','.join(cedge.receiver_node.get_leaf_names()))
        print(f'\tAffected words  ({len(words)}):')
        for word in words:
            print('\t\t', word)


        plot_network(nexus, annotate_leafs=False)
        plot_wordtree(nexus, word=words[0])
        plot_word_labels(nexus, words=words)
        plt.axis('off')
        plt.tight_layout(pad=0)
            # plt.grid(axis='y')
            # ax.set_xticklabels([])
            # ax.set_yticklabels([])
            # ax.set_frame_on(False)
            # ax.tick_params(tick1On=False)

        plt.savefig(word_trees_directory / f'wordtree_{words[0]}_et_al.pdf')
        # plt.show()


def plot_contactedge_with_data(summary_nexus: NexusReader,
                               samples_nexus: NexusReader,
                               edge_trees_directory=Path('./wordtrees/'),
                               burnin=0.1, summary_distance_threshold=4):
    """Plot the contacTree with one edge highlighted and the words borrowed at this edge.
    Store the plots as PDF files in `edge_trees_directory`."""
    # Create the word-trees directory
    os.makedirs(edge_trees_directory, exist_ok=True)

    # Read the contact edges and attach them to nodes in the `nexus` tree.
    tree = read_tree(summary_nexus)
    contactedges = read_contactedges(summary_nexus, burnin=burnin)
    place_contactedges_in_tree(tree, contactedges)

    cedge_samples, n_samples = read_contactedges(samples_nexus, burnin=burnin, return_n_samples=True)
    place_contactedges_in_tree(tree, cedge_samples)

    words_by_cedge = defaultdict(list)

    for cedge in contactedges:
        for w in cedge.affected_blocks:
            words_by_cedge[cedge].append(w)

    cedge_samples_by_receiver = defaultdict(list)
    for cedge in cedge_samples:
        cedge_samples_by_receiver[cedge.receiver_node.name].append(cedge)

    for cedge, words in [*words_by_cedge.items()]:
        donor_clade = cedge.donor_node.get_leaf_names()
        receiver_clade = cedge.receiver_node.get_leaf_names()

        print()
        fig_width = 6 + 0.8 * len(tree.get_leaf_names())
        fig_height = 16 + (0.25 * len(words))
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        print('Plotting marginal tree for edge:')
        print('\t\t', ','.join(cedge.donor_node.get_leaf_names()),
              '->', ','.join(cedge.receiver_node.get_leaf_names()))
        print(f'\tAffected words  ({len(words)}):')
        if len(words) > 80:
            continue
        for word in words:
            print('\t\t', word)

        assign_node_coordinates(tree)
        plot_contactedge_arrow(cedge, lw=40, c=CEDGE_COLOR)
        plot_tree_topology(tree, lw=5, color='white', zorder=ZOrder.LANGUAGE_TREE_FRAME)
        plot_tree_topology(tree, lw=4, color=GRAY_3, zorder=ZOrder.LANGUAGE_TREE)

        cedge_samples = cedge_samples_by_receiver[cedge.receiver_node.name]

        # Filter by receiver edge (should match the summary criterion)

        def node_distance(n1: Node, n2: Node):
            if n1 is n2: return 0
            if get_age(n1) < get_age(n2): return node_distance(n1.ancestor, n2) + 1
            else: return node_distance(n1, n2.ancestor) + 1

        def similar_donor(other):
            return node_distance(cedge.donor_node, other.donor_node) < summary_distance_threshold

        cedge_samples = list(filter(similar_donor, cedge_samples))

        donor_count = defaultdict(int)
        x_donor = []
        x_receiver = []
        y = []
        for c in cedge_samples:
            donor_count[c.donor_node] += 1
            x_donor.append(c.donor_node.x)
            x_receiver.append(c.receiver_node.x)
            y.append(c.height)

        plt.scatter(x_donor, y, s=50, color=DONOR_COLOR, alpha=0.3, lw=0.2, marker='_', zorder=ZOrder.WORD_TREE)
        plt.scatter(x_receiver, y, s=50, color=RECEIVER_COLOR, alpha=0.3, lw=0.2, marker='_', zorder=ZOrder.WORD_TREE)
        plt.scatter([], [], color='teal', marker='_', label='Donor end of a contact edge')
        plt.scatter([], [], color='magenta', marker='_', label='Receiver end of a contact edge')

        for node, count in donor_count.items():
            if node.ancestor is None: continue
            frequency = count / n_samples
            if frequency > 0.02:
                plt.text(
                    x=node.x + 0.2,
                    y=0.5 * (node.y + node.ancestor.y),
                    s='%.2f' % frequency,
                    color=DONOR_COLOR,
                    fontsize=12,
                    horizontalalignment='left',
                    verticalalignment='center',
                )

        # Plot internal nodes
        # x = [node.x for node in tree.walk() if not node.is_leaf]
        # y = [node.y for node in tree.walk() if not node.is_leaf]
        # plt.scatter(x, y, s=100, color=GRAY_3, zorder=ZOrder.LANGUAGE_TREE)

        add_grid_lines(plt.gca(), max_y=get_age(tree))

        plot_word_labels(nexus=summary_nexus, words=words, word_support=cedge.block_posterior,
                         donor_clade=donor_clade, receiver_clade=receiver_clade, ax=ax)
        plt.axis('off')
        plt.tight_layout()

        donor_clade_str = '|'.join(donor_clade)
        receiver_clade_str = '|'.join(receiver_clade)
        plt.savefig(edge_trees_directory / f'from_[{donor_clade_str}]_to_[{receiver_clade_str}].pdf')
        # plt.show()


def add_grid_lines(ax, max_y=5000.0, dy=1000.0):
    for y in np.arange(dy, max_y + 1, dy, dtype=int):
        ax.axhline(y, color='gray', lw=0.25, zorder=ZOrder.GRID_LINE)
        ax.text(
            -1.5, y - 20, str(y),
            horizontalalignment='center',
            verticalalignment='top',
            fontsize=10,
            color='lightgray'
        )

def plot_contact_to(receiver: str, summary_nexus: NexusReader, posterior_nexus: NexusReader):
    print(posterior_nexus.trees.translators)
    print(summary_nexus.trees.translators)
    for k, v in list(posterior_nexus.trees.translators.items()):
        posterior_nexus.trees.translators[k] = f'{k}_{v}'
        summary_nexus.trees.translators[k] = f'{k}_{v}'

    # Function to check whether a contact edge leads to a lineage above the specified ´receiver´ node.
    def to_receiver(cedge: ContactEdge):
        # return receiver in cedge.receiver_node.get_leaf_names()
        return cedge.receiver_node.name == receiver

    # Useful constants
    n_samples = len(posterior_nexus.trees.trees)

    # Read summary and posterior ctrees
    tree = read_tree(summary_nexus)
    contactedges = read_contactedges(posterior_nexus)
    place_contactedges_in_tree(tree, contactedges)
    contactedges = list(filter(to_receiver, contactedges))

    assign_node_coordinates(tree)
    plot_tree_topology(tree, annotate_leafs=True, lw=2, color='k', zorder=ZOrder.LANGUAGE_TREE)

    # # Plot internal nodes
    # x = [node.x for node in tree.walk() if not node.is_leaf] + [c.donor_node.x for c in contactedges]
    # y = [node.y for node in tree.walk() if not node.is_leaf] + [c.height for c in contactedges]
    # plt.scatter(x, y, s=100, color=GRAY_2, zorder=ZOrder.LANGUAGE_TREE)

    x_donor = [c.donor_node.x for c in contactedges]
    x_receiver = [c.receiver_node.x for c in contactedges]
    y = [c.height for c in contactedges]
    # plt.scatter(x_donor, y, s=50, color='teal', alpha=0.02, marker='_', zorder=ZOrder.WORD_TREE)
    # plt.scatter(x_receiver, y, s=50, color='magenta', alpha=0.02, marker='_', zorder=ZOrder.WORD_TREE)
    plt.scatter(x_donor, y, s=50, color='teal', alpha=1.0, lw=0.1, marker='_', zorder=ZOrder.WORD_TREE)
    plt.scatter(x_receiver, y, s=50, color='magenta', alpha=1.0, lw=0.1, marker='_', zorder=ZOrder.WORD_TREE)
    plt.scatter([], [], color='teal', marker='_', label='Donor end of a contact edge')
    plt.scatter([], [], color='magenta', marker='_', label='Receiver end of a contact edge')

    donor_count = defaultdict(int)
    for cedge in contactedges:
        donor_count[cedge.donor_node] += 1

    for node, count in donor_count.items():
        if node.ancestor is None: continue
        frequency = count / n_samples
        if frequency > 0.02:
            plt.text(
                x=node.x + 0.2,
                y=0.5 * (node.y + node.ancestor.y),
                s='%.2f' % frequency,
                color='teal',
                fontsize=12
            )

    add_grid_lines(plt.gca(), max_y=get_age(tree))
    plt.legend()

    # plot_word_labels(summary_nexus)

if __name__ == '__main__':
    from contacTreesIE.plotting import plot_tree_topology, assign_node_coordinates

    # fig, ax = plt.subplots(figsize=(30, 95))

    # TREES_PATH = Path('/home/nico/workspace/IE_case_study/small_runs/CT_fixTopo/noPreservedLatin/CT_fixTopo.trees')
    # TREES_PATH = Path('results/workstation/CT_fixedTopology_preservedLatin/CT_fixTopo_pL.trees')
    # SUMMARYTREE_PATH = Path('results/workstation/CT_fixedTopology_preservedLatin/CT_fixTopo_pL.summary.tree')
    # TREES_PATH = Path('results/workstation/CT_fixTopo/covarion/CT_fixTopo_covarion_1.trees')
    # SUMMARYTREE_PATH = Path('results/workstation/CT_fixTopo/covarion/CT_fixTopo_covarion_1.summary.tree')
    # TREES_PATH = Path('small_runs/CT_fixTopo/covarion/CT_fixTopo_covarion_1.trees')
    # SUMMARYTREE_PATH = Path('small_runs/CT_fixTopo/covarion/CT_fixTopo_covarion_1.summary.tree')
    TREES_PATH = Path('results/fixTopo_covarion.trees')
    SUMMARYTREE_PATH = Path('results/fixTopo_covarion.summary.tree')
    # samples_nexus = NexusReader.from_file(TREES_PATH)
    summary_nexus = NexusReader.from_file(SUMMARYTREE_PATH)

    fig, ax = plt.subplots(figsize=(60, 90))

    # plot_network(summary_nexus, annotate_leafs=False)
    # plot_wordtree(summary_nexus, word='lake')
    plot_word_labels(
        summary_nexus,
        header_interval=6,
        no_tree=True,
        ax=ax
    )
    # plot_word_labels(summary_nexus)

    # plot_contact_to('2_Danish', summary_nexus, posterior_nexus)
    # plot_contact_to('2_English', summary_nexus, posterior_nexus)
    # plot_contact_to('3_German', summary_nexus, posterior_nexus)

    # plot_contactedge_with_data(summary_nexus=summary_nexus,
    #                            samples_nexus=samples_nexus,
    #                            edge_trees_directory=Path('./loanwords_per_contact_edge/'))


    plt.axis('off')
    plt.tight_layout(pad=0.01)
    plt.savefig('wordlist.pdf')

    # plot_wordtrees(summary_nexus, word_trees_directory=Path('./wordtrees_tmp/'))
