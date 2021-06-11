#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, \
    unicode_literals, annotations

import os
from pathlib import Path
from typing import List
from operator import attrgetter
from collections import defaultdict
from enum import IntEnum

import numpy as np
import pandas as pd
from newick import Node
from nexus import NexusReader

import matplotlib.pyplot as plt

from contacTreesIE.newick_util import get_age, ContactEdge, get_actual_leaves, collect_contactedges, drop_contactedges, \
    translate_node_names, remove_dead_end



class ZOrder(IntEnum):

    GRID_LINE = 1
    CONTACT_EDGE = 2
    LANGUAGE_TREE_FRAME = 3
    CONTACT_EDGE_STUB = 4
    LANGUAGE_TREE = 5
    WORD_TREE = 6

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


def read_contactedges(nexus: NexusReader) -> List[ContactEdge]:
    translate = nexus.trees.translators
    trees: List[Node] = [t.newick_tree for t in nexus.trees.trees]

    contactedges = []
    for tree in trees:
        translate_node_names(tree, translate)
        contactedges += collect_contactedges(tree)

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
        contactedges: List[ContactEdge]
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


def plot_contactedges(contactedges, word=None, **plot_kwargs):
    plot_kwargs['color'] = plot_kwargs.pop('c', 'k')
    # width = plot_kwargs.pop('lw', 1.0)
    # head_width = 2. * width

    if word is not None:
        contactedges = filter(lambda ce: word in ce.affected_blocks, contactedges)

    for cedge in contactedges:
        x1 = cedge.donor_node.x
        x2 = cedge.receiver_node.x
        y = cedge.height
        direction = np.sign(x2 - x1)
        plt.plot([x1, x2], [y, y], zorder=ZOrder.CONTACT_EDGE, **plot_kwargs)
        plt.plot([x2 - direction*0.5, x2], [y, y], zorder=ZOrder.CONTACT_EDGE_STUB, **plot_kwargs)
        # plt.arrow(x=x1, y=y, dx=(x2-x1), dy=0, width=width,
        #           length_includes_head=True, head_width=head_width, head_length=0.3,
        #           ec=(0.65, 0.65, 0.65), **plot_kwargs)


def plot_densiedges(trees_nexus, summarytree_nexus):
    contactedges = read_contactedges(nexus=trees_nexus)
    tree = read_tree(nexus=summarytree_nexus)
    place_contactedges_in_tree(tree, contactedges)

    assign_node_coordinates(tree)
    plot_tree_topology(tree, annotate_leafs=True)
    plot_contactedges(contactedges, word='lake', alpha=0.2)


GRAY_2 = (0.65, 0.65, 0.65)

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
    # return lexeme.replace('  ', ' ').replace('  ', ' ').replace(' ', '_').lower()
    return lexeme.replace('  ', ' ').replace('  ', ' ').lower()


def plot_word_labels(nexus: NexusReader, words=None):
    tree = read_tree(nexus)
    drop_contactedges(tree)
    assign_node_coordinates(tree)

    # Load ielex dataset for lexemes
    df = pd.read_csv('resources/ielex-130421-ag-cc.tsv', sep='\t', dtype=str)
    df['concept'] = df['cc_alias'].apply(lambda s: s.split('-')[0])

    def get_word_in_language(language, concept):
        x = df.loc[(df.language == language) & (df.concept == concept)]
        if len(x) == 0:
            return '?'

        # assert len(x.lexeme) == len(set(x.lexeme.to_list()))
        # I wanted to assert that lexemes are unique, but this seems to be intended:
        # E.g. `karlmaðr` in Old Norse is in two cognate classes with the same meaning,
        # since `karl-` (male) and `-maðr` (man) have separate etymological origins.

        lexeme = '\n'.join(x.lexeme.to_numpy())
        lexeme = clean_lexeme(lexeme)

        return lexeme

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
        fontsize=8,
    )

    Y0 = -100.0
    DY = 250.0
    for node in tree.get_leaves():
        # Plot a visual connection from label to leaf, if leaf isn't contemporary
        if get_age(node) > 0.001:
            plt.plot(
                [node.x, node.x], [0, node.y],
                color='lightgray', ls='dotted',
                zorder=0
            )
        plt.text(x=node.x, y=Y0, s=node.name, weight='bold', **text_args)
        for i, word in enumerate(words):
            plt.text(
                x=node.x,
                y=Y0 - (1 + i) * DY,
                s=get_word_in_language(node.name, word),
                **text_args
            )
            c = get_cognate_colors(node.name, word)
            plt.scatter(
                x=node.x + np.linspace(-0.3, 0.3, len(c)+2)[1:-1],
                y=[Y0 - (1 + i) * DY + 50] * len(c),
                c=c,
                s=15,
                marker='s'
            )

    for i, word in enumerate(words):
        plt.text(
            x=-1.2,
            y=Y0 - (1 + i) * DY,
            s=word,
            weight='bold',
            **text_args
        )


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
        fig, ax = plt.subplots(figsize=(40, 16))

        if int(g) == 0:
            continue

        print('Plotting marginal tree', g)
        print('Contact edges:')
        for i, cedge in enumerate(contactedges):
            if g[i] == '1':
                print('\t', ','.join(cedge.donor_node.get_leaf_names()),
                      '->', ','.join(cedge.receiver_node.get_leaf_names()))
        print(f'Affected words  ({len(words)}):')
        for word in words:
            print('\t', word)


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


if __name__ == '__main__':
    from contacTreesIE.plotting import plot_tree_topology, assign_node_coordinates

    TREES_PATH = Path('ie.trees')
    SUMMARYTREE_PATH = Path('ie.summary.tree')
    posterior_nexus = NexusReader.from_file(TREES_PATH)
    summary_nexus = NexusReader.from_file(SUMMARYTREE_PATH)

    # plot_network(summary_nexus, annotate_leafs=False)
    # plot_wordtree(summary_nexus, word='lake')
    # plot_word_labels(summary_nexus, words=['lake', 'animal'])
    #
    # plt.axis('off')
    # plt.tight_layout(pad=0)
    # plt.show()

    plot_wordtrees(summary_nexus)
