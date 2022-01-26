from collections import deque
from typing import List
from itertools import product

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from nexus import NexusReader
from newick import Node

from contacTreesIE.compile_beast_xmls import read_dataset
from contacTreesIE.plotting.plot_contact_edges import read_contactedges, place_contactedges_in_tree
from contacTreesIE.newick_util import collect_contactedges, drop_contactedges, ContactEdge


class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'


def plot_precision_recall(loans: pd.DataFrame):
    fig, axes = plt.subplots(1, 3)

    for lang in loans.LANGUAGE.unique():
        x = loans[loans.LANGUAGE == lang]

        x.insert(0, 'idx', np.arange(len(x)))
        is_zero = x.p == 0.0
        if any(is_zero):
            cutoff = x[is_zero].iloc[0].idx
            print(lang, cutoff)
        else:
            print(lang, 'no nulls')

        tru_pos = np.cumsum(x.truth)
        tru = np.sum(x.truth)
        pos = np.arange(len(x.truth)) + 1
        recall = tru_pos / tru
        precision = tru_pos / pos

        if tru >= 8:
            axes[0].plot(pos, recall, label=lang, lw=1)
            axes[1].plot(pos, precision, label=lang, lw=1)
            axes[2].plot(recall, precision, label=lang, lw=1)

    for title, ax in zip(['Recall', 'Precision', 'Recall -- Precision'], axes):
        ax.legend()
        ax.set_ylim(0, 1)
        ax.set_title(title)

    axes[0].plot([0, 206], [0, 1], 'lightgrey', zorder=0)

    axes[0].set_xlim(0, 206)
    axes[1].set_xlim(0, 206)
    axes[2].set_xlim(0, 1)
    plt.tight_layout(pad=0.05, w_pad=0.0)
    plt.show()


def plot_receiver_operator_curve(loans: pd.DataFrame, show_example_threshold=True,
                                 color_cycle=None, ax=None, min_pos=10,
                                 show_label=True):
    if ax is None:
        ax = plt.gca()
    if color_cycle is None:
        color_cycle = deque(plt.rcParams['axes.prop_cycle'].by_key()['color'])

    for lang in loans.LANGUAGE.unique():
        x = loans[loans.LANGUAGE == lang]
        x.insert(0, 'idx', np.arange(len(x)))

        pos = np.sum(x.truth)
        neg = np.sum(~x.truth)

        if pos <= min_pos:
            continue

        color = color_cycle.popleft()

        tpr = []
        fpr = []
        first = True
        for i, p_thresh in reversed([*enumerate(x.p)]):
            positive_predictions = x.iloc[:i]
            tp = np.sum(positive_predictions.truth)
            fp = np.sum(~positive_predictions.truth)
            tpr.append(tp / pos)
            fpr.append(fp / neg)
            if first and p_thresh >= 0.33:
                print(lang, pos)
                print('%.2f:' % p_thresh, tpr[-1], fpr[-1])
                print()
                first = False

                if show_example_threshold:
                    plt.scatter([fpr[-1]], [tpr[-1]], color=color)
                    is_left = lang in ['Breton_ST', 'Welsh_N']
                    d = 1 if is_left else -1
                    t = plt.annotate(f'FPR = {100*fpr[-1]:.1f}%\nTPR = {100*tpr[-1]:.1f}%', (fpr[-1], tpr[-1]),
                                     va='center', ha='left' if is_left else 'right',
                                     color=color, weight='bold', xytext=(d*14, -d*8), textcoords='offset points')
                    t.set_bbox(dict(facecolor=(1.0, 1.0, 1.0, 1),
                                    edgecolor=color,
                                    lw=1.5))

        plt.plot(fpr, tpr, color=color, label=(lang if show_label else None))

    # ax.set_title('receiver-operator curve', fontsize=18, fontweight="bold",)
    ax.set_xlabel('false positive rate', fontsize=16)
    ax.set_ylabel('true positive rate', fontsize=16)
    ax.plot([0, 1], [0, 1], 'lightgrey', zorder=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend()


def get_IE_loans():
    data = read_dataset('resources/data-mittellatein-2021-09-30.csv')
    data = data.rename({'concept': 'WORD', 'language': 'LANGUAGE'}, axis=1)
    data = data.astype({'lexeme': str})

    loan_lexemes = data.loc[data.status.isin(['LOAN'])]
    loan_lexemes.rename({'comments 1': 'donor'}, axis=1, inplace=True)
    loan_lexemes.drop(['cc_alias', 'cc_id', 'lexeme_id', 'status', 'comments 2', 'comments 3'], axis=1, inplace=True)
    loan_lexemes = loan_lexemes.groupby(['LANGUAGE', 'WORD'])[['lexeme', 'donor']].apply(lambda col: col.apply('/'.join, axis=0)).reset_index()
    loan_lexemes = loan_lexemes.set_index(['LANGUAGE', 'WORD'])
    loan_lexemes.to_csv('ielex_loanwords.csv', index=True)

    data = data.drop(['cc_alias', 'cc_id', 'lexeme_id', 'status', 'comments 1', 'comments 2', 'comments 3'], axis=1)
    data = data.groupby(['LANGUAGE', 'WORD'])['lexeme'].apply(' | '.join).reset_index()
    data = data.set_index(['LANGUAGE', 'WORD'])
    data.lexeme[loan_lexemes.index] = loan_lexemes.lexeme
    data.lexeme.fillna('?', inplace=True)

    # Read known and predicted loanwords
    loans_known = pd.read_csv('loans.csv', index_col=0)
    loans_reconstructed = pd.read_csv('results/fix_clock_stdev/CT_fixTopo/covarion/loanwords.log', index_col=0)
    # loans_reconstructed = pd.read_csv('results/test/covarion_1/loanwords.log', index_col=0)

    # Copy 'WORD' from index into a separate column
    loans_known['WORD'] = loans_known.index
    loans_reconstructed['WORD'] = loans_reconstructed.index

    # Change data-frame format from 'words as rows, languages as columns' to
    # '(word, language) as rows, cell values in separate column'
    loans_known = pd.melt(loans_known, id_vars='WORD', var_name='LANGUAGE', value_name='truth')
    loans_reconstructed = pd.melt(loans_reconstructed, id_vars='WORD', var_name='LANGUAGE', value_name='p')

    # Join grount-truth and reconstruction into one data-frame
    loans = loans_known.join(loans_reconstructed.set_index(['LANGUAGE', 'WORD']), on=['LANGUAGE', 'WORD'])
    # Treat NAs as 0.0
    loans.p = loans.p.fillna(0.0)
    # Sorting by 'p' helps with computing FPs, FNs, etc for different thresholds.
    loans = loans.sort_values('p', ascending=False)
    # Join the original 'data' data-frame (for printing include lexemes)
    loans = loans.join(data, on=['LANGUAGE', 'WORD'], how='left', rsuffix='_')

    return loans


def read_trees(nexus: NexusReader) -> List[Node]:
    trees = [t.newick_tree for t in nexus.trees.trees]
    return trees


def read_contactedges(nexus: NexusReader, burnin=0.1,
                      tree_idx=None,
                      block_posterior_threshold=0.0) -> List[ContactEdge]:
    tree_handler = nexus.trees

    if tree_idx:
        tree = tree_handler.trees[tree_idx].newick_tree
        return collect_contactedges(tree, block_posterior_threshold=block_posterior_threshold)

    n_burnin_samples = int(burnin * tree_handler.ntrees)
    trees = [t.newick_tree for t in tree_handler.trees[n_burnin_samples::10]]
    contactedges = []
    for tree in trees:
        contactedges += collect_contactedges(tree, block_posterior_threshold=block_posterior_threshold)

    return contactedges

BLOCKS = [f'block{i}' for i in range(100)]
def get_loanwords(tree, cedges: List[ContactEdge]):
    df = pd.DataFrame(index=BLOCKS, columns=tree.get_leaf_names())
    df = df.fillna(0).astype(int)

    for cedge in cedges:
        leaves = cedge.receiver_node.get_leaf_names()
        for w_, lang in product(cedge.affected_blocks, leaves):
            word = w_.split('.')[0]
            df.loc[word, lang] = 1

    print(df)
    return df


def get_simulation_loans():
    # Load simulated trees and compute place
    truth_nexus = NexusReader.from_file('results/simulation_study/sampling/SAMPLING.trees')
    tree_1 = read_trees(truth_nexus)[5]
    drop_contactedges(tree_1)
    cedges_1 = read_contactedges(truth_nexus, tree_idx=5)
    place_contactedges_in_tree(tree_1, cedges_1)
    loans_known = get_loanwords(tree_1, cedges_1)

    # nexus_s1 = NexusReader.from_file('results/simulation_study/contactrees/test.001.trees')
    # tree_1_samples = read_trees(nexus_s1)
    loans_reconstructed = pd.read_csv('results/simulation_study/contactrees/loanwords.log', index_col=0)

    # Copy 'WORD' from index into a separate column
    loans_known['WORD'] = loans_known.index
    loans_reconstructed['WORD'] = loans_reconstructed.index

    # Change data-frame format from 'words as rows, languages as columns' to
    # '(word, language) as rows, cell values in separate column'
    loans_known = pd.melt(loans_known, id_vars='WORD', var_name='LANGUAGE', value_name='truth')
    loans_reconstructed = pd.melt(loans_reconstructed, id_vars='WORD', var_name='LANGUAGE', value_name='p')

    # Join grount-truth and reconstruction into one data-frame
    loans = loans_known.join(loans_reconstructed.set_index(['LANGUAGE', 'WORD']), on=['LANGUAGE', 'WORD'])
    # Treat NAs as 0.0
    loans.p = loans.p.fillna(0.0)
    # Sorting by 'p' helps with computing FPs, FNs, etc for different thresholds.
    loans = loans.sort_values('p', ascending=False)

    return loans


if __name__ == '__main__':
    loans = get_IE_loans()
    loans_simu = get_simulation_loans()

    fig, ax = plt.subplots(1, 1, figsize=(13, 7))

    # plot_precision_recall(loans)
    plot_receiver_operator_curve(loans, ax=ax, min_pos=7)

    plt.tight_layout(pad=0.05, w_pad=0.0)
    plt.show()

    #
    # # Write significant false positives/negatives to file:
    #
    # false_pos = loans[~loans.truth & (loans.p > 0.8)]
    # false_neg = loans[loans.truth & (loans.p < 0.2)]
    # false_pos['formatted'] = false_pos.apply(lambda row: '{p:.2f}: {lexeme} ({WORD})'.format(**row.to_dict()), axis=1)
    # false_neg['formatted'] = false_neg.apply(lambda row: '{p:.2f}: {lexeme} ({WORD})'.format(**row.to_dict()), axis=1)
    # with open('false_positives.txt', 'w') as fn_file:
    #     for lang, fp in false_pos.groupby('LANGUAGE'):
    #         fn_file.write(lang + '\n    ')
    #         fn_file.write('\n    '.join(fp.formatted) + '\n')
    # with open('false_negatives.txt', 'w') as fn_file:
    #     for lang, fn in false_neg.groupby('LANGUAGE'):
    #         fn_file.write(lang + '\n    ')
    #         fn_file.write('\n    '.join(fn.formatted) + '\n')
