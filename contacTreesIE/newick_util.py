from operator import attrgetter
from typing import Union, List, Dict, Optional

from newick import Node


class ContactEdge(object):

    """Basic class to bundle information about a contact edge.

    Attributes:
        donor_clade (str): String representation of the clade where the edge starts
                           (the donor language of the loan words).
        receiver_clade (str): String representation of the clade where the edge ends
                              (the receiver language of the loan words).
        height (float): The height (units of time before present) of the contact event
        affected_blocks (list): ...
        block_posterior (dict): ...
        donor_node (Node): ...
        receiver_node (Node): ...

    """

    def __init__(self,
                 donor_clade: str = None,
                 receiver_clade: str = None,
                 height: float = None,
                 affected_blocks: List[str] = None,
                 block_posterior: Dict[str, float] = None):
        self.donor_clade = donor_clade
        self.receiver_clade = receiver_clade
        self.height = height
        self.affected_blocks = affected_blocks
        self.block_posterior = block_posterior

        self.donor_node = None
        self.receiver_node = None

    def __repr__(self):
        return f'ContactEdge({self.donor_clade}->{self.receiver_clade}:{self.height})'


def get_root(node: Node) -> Node:
    if node.ancestor is None:
        return node
    else:
        return get_root(node.ancestor)


def get_sibling(node: Node) -> Optional[Node]:
        if node.ancestor is None or len(node.ancestor.descendants) != 2:
            return None

        me_and_my_sibling = node.ancestor.descendants
        if node is me_and_my_sibling[0]:
            return me_and_my_sibling[1]
        else:
            return me_and_my_sibling[0]


def get_depth(node: Node) -> float:
    if node.ancestor is None:
        return 0.0
    else:
        return get_depth(node.ancestor) + node.length


def get_height(node: Union[ContactEdge, Node]) -> float:
    if isinstance(node, ContactEdge):
        return node.height

    if node.is_leaf:
        return 0.0
    else:
        return max(get_height(child) + child.length for child in node.descendants)


def get_age(node: Union[ContactEdge, Node]) -> float:
    if isinstance(node, ContactEdge):
        return node.height

    if node.ancestor is None:
        return get_height(node)
    else:
        return get_age(node.ancestor) - node.length


def get_node_by_name(tree: Node, name: str, ignore_duplicates: bool = False) -> Node:
    matches = get_all_nodes_named(tree, name)

    assert len(matches) >= 1, f'No node named `{name}` found in tree.'

    if not ignore_duplicates:
        assert len(matches) == 1, f'Multiple nodes named `{name}` found in tree.'

    return matches[0]


def get_all_nodes_named(tree: Node, name: str):
    matches = []
    for node in tree.walk():
        if node.name == name:
            matches.append(node)

    return matches

def parse_node_comment(comment: str) -> dict:
    # Drop the leading `&`
    assert comment.startswith('&')
    comment = comment[1:]

    # Read the attributes as a dict (using pythons dict() notation)
    attrs = eval(f'dict({comment})')

    # Manually parse the `affectedBlocks` string
    attrs['affectedBlocks'] = set(attrs['affectedBlocks'][1:-1].split(','))

    if 'blockPosterior' in attrs and isinstance(attrs['blockPosterior'], str):
        block_post_reformated = attrs['blockPosterior'].replace(':', '\':')\
                                                       .replace(',', ',\'')\
                                                       .replace('{', '{\'')
        attrs['blockPosterior'] = eval(block_post_reformated)
    else:
        attrs['blockPosterior'] = {b: 1. for b in attrs['affectedBlocks']}

    return attrs


def get_actual_leaves(node: Node) -> List[Node]:
    return [n for n in node.get_leaves() if not n.name.startswith('#')]


def collect_contactedges(tree: Node, block_posterior_threshold: float = 0.5) -> List:
    contactedges = {}

    leaves = sorted(get_actual_leaves(tree), key=attrgetter('name'))

    def compress_clade(node):
        clade = get_actual_leaves(node)
        return ''.join([str(int(l in clade)) for l in leaves])

    done = set()

    for node in tree.walk():
        if node.name is None:
            assert not node.is_leaf
            continue

        # if node.name and node.name.startswith('#'):
        if node.name.startswith('#'):
            if node.name not in contactedges:
                contactedges[node.name] = ContactEdge()
            cedge = contactedges[node.name]

            # The donor node contains information about the contactedge in the comment
            # The receiver node has no comments
            if node.comment and '&conv' in node.comment:
                assert node.is_leaf
                cedge.donor_clade = compress_clade(node.ancestor)
                node_attrs = parse_node_comment(node.comment)
                cedge.height = get_age(node)
                cedge.affected_blocks = [
                    b for b in node_attrs['affectedBlocks']
                    if node_attrs['blockPosterior'][b] > block_posterior_threshold
                ]
                cedge.block_posterior = {b: node_attrs['blockPosterior'][b] for b in node_attrs['affectedBlocks']}

            else:
                assert not node.is_leaf
                cedge.receiver_clade = compress_clade(node)

            cedge_hash = (cedge.donor_clade, cedge.receiver_clade)
            if None not in cedge_hash:
                # Check whether this is contact edge was already added
                if cedge_hash in done:
                    contactedges.pop(node.name)
                else:
                    done.add(cedge_hash)

    return list(contactedges.values())


def drop_contactedges(tree: Node):
    fake_tips = [n for n in tree.get_leaves() if n.name.startswith('#')]
    tree.prune(fake_tips)
    tree.remove_redundant_nodes()


def translate_node_names(tree: Node, translate: dict):
    def rename_node(node: Node):
        node.name = translate.get(node.name, node.name)
    tree.visit(rename_node)


def is_single_child(node: Node) -> bool:
    assert node in node.ancestor.descendants
    return len(node.ancestor.descendants) == 1


def remove_dead_end(node: Node):
    while node.ancestor and is_single_child(node):
        node.ancestor.descendants.remove(node)
        node = node.ancestor
