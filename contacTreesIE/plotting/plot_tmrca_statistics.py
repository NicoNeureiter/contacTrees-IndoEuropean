from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sb

from newick import Node
from nexus import NexusReader

from contacTreesIE.plotting import color_tree, plot_tree_topology


def get_node_depths(root:Node, depth: float, depth_dict: dict):
    depth_dict[root] = depth
    for child in root.descendants:
        get_node_depths(child, depth + child.length, depth_dict)
    return depth_dict


def mrca(n1: Node, n2: Node, depths: dict):
    if n1 is n2:
        return n1

    if (depths[n1] > depths[n2]):
        return mrca(n1.ancestor, n2, depths)
    else:
        return mrca(n1, n2.ancestor, depths)


def plot_tmrca_statistics(trees_path: Path, label=None, axes=None):
    if axes is None:
        _, axes = plt.subplots(2)

    nex = NexusReader.from_file(trees_path)
    name_to_idx = {name: str(idx) for idx, name in enumerate(nex.taxa.taxa)}
    i_en = name_to_idx['English']
    i_de = name_to_idx['German']
    i_fr = name_to_idx['French']

    idx_1 = i_de
    idx_2 = i_en

    tmrcas = []
    for tree_line in nex.trees[100:104:2]:
    # for tree_line in nex.trees[1:5]:
        tree = tree_line.newick_tree
        depths = get_node_depths(tree, 0., {})
        max_depth = max(depths.values())
        node_1 = tree.get_node(idx_1)
        node_2 = tree.get_node(idx_2)

        color_tree(tree)
        plot_tree_topology(tree, ax=axes[1])

        tmrcas.append(
            max_depth - depths[mrca(node_1, node_2, depths)]
        )
        # print(tmrcas)
        # exit()

    print(len(tmrcas))
    sb.kdeplot(tmrcas, label=label, bw=0.2, ax=axes[0])
    # plt.hist(tmrcas, bins=50, label=label)


if __name__ == '__main__':
    PATH_CT = Path('results/topology_fixed/ielex_contactrees.trees')
    PATH_BT = Path('results/topology_fixed/ielex_basictrees.trees')
    PATH_WORD_TREES = Path('results/topology_fixed/wordtrees/')
    # PATH_WORD_TREES = Path('../contactrees/IE_run_20210128/')

    fig, axes = plt.subplots(2)
    plot_tmrca_statistics(PATH_BT, label='basic tree', axes=axes)
    plot_tmrca_statistics(PATH_CT, label='contactree', axes=axes)
    plt.legend()
    plt.show()

    # for word in ['all', 'animal', 'bark', 'lake', 'count']:
    for word in ['all', 'one', 'two', 'lake', 'count', 'animal', 'because']:
        plot_tmrca_statistics(PATH_WORD_TREES / f'ielex_contactrees.{word}.trees', label=word)

    plt.legend()
    plt.show()
