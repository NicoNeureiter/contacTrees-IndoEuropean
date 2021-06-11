#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
from collections import defaultdict
from pathlib import Path
from typing import List

import pandas as pd
from newick import Node, parse_node

from contacTreesIE.newick_util import get_sibling, translate_node_names, get_age
from contacTreesIE.pre.language_lists import *
from contacTreesIE.pre import xml_snippets
from contacTreesIE.pre.starting_trees import CHANG_MEDIUM_TREE
from contacTreesIE.pre.starting_trees import CHANG_MEDIUM_TRANSLATE

XML_TEMPLATE_PATH = 'resources/ie_template.xml'
LOANS_CSV_PATH = 'loans.csv'


DATASET = CHANG_MEDIUM
INCLUDED_CLADES = [*CELTIC, *GERMANIC, *ROMANCE]
INCLUDED_LANGUAGES = [l for l in DATASET if l in INCLUDED_CLADES]
RENAME = RENAME_CHANG
translate = CHANG_MEDIUM_TRANSLATE
STARTING_TREE = CHANG_MEDIUM_TREE

TIP_DATES = defaultdict(float)
TIP_DATES.update({
    "Old_Irish": 1250.0,
    "Cornish": 300.0,
    "Umbrian": 2200.0,
    "Oscan": 2200.0,
    "Latin": 2050.0,
    "Gothic": 1650.0,
    "Old_Norse": 775.0,
    "Old_High_German": 1050.0,
    "Old_English": 1000.0,
    "Old_Prussian": 500.0,
    "Old_Church_Slavonic": 1000.0,
    "Ancient_Greek": 2400.0,
    "Classical_Armenian": 1450.0,
    "Avestan": 2500.0,
    "Old_Persian": 2450.0,
    "Vedic_Sanskrit": 3000.0,
    "Tocharian_A": 1375.0,
    "Tocharian_B": 1350.0,
    "Hittite": 3450.0,
    "Luvian": 3350.0,
    "Lycian": 2400.0
})

SAMPLED_ANCESTORS = [
    'Latin',
    'Old_Irish',
    'Old_High_German',
    'Old_English',
    'Old_Norse'
]


class Clade(object):

    clade_template = '''\
        <distribution id="{clade}.prior" spec="MRCAPrior" tree="@acg" monophyletic="true">
          <taxonset id="{clade}" spec="TaxonSet">
{taxa}
          </taxonset>
{age_prior}        </distribution>\n\n'''

    taxon_template = '            <taxon idref="{taxon}" />'

    def __init__(self, name, members, age_prior=''):
        self.name = name
        self.members = [m for m in members
                        if isinstance(m, Clade) or (m in INCLUDED_LANGUAGES)]
        # self.members = members
        self.age_prior = age_prior

    def get_languages(self) -> List[str]:
        languages = []
        for m in self.members:
            if isinstance(m, Clade):
                languages += m.get_languages()
            else:
                assert isinstance(m, str)
                languages.append(m)
        return languages

    def to_xml(self):
        taxa_str = '\n'.join(
            self.taxon_template.format(taxon=l) for l in self.get_languages()
        )

        age_prior = ''
        if self.age_prior:
            age_prior = (' ' * 10) + self.age_prior + '\n'

        return self.clade_template.format(clade=self.name, taxa=taxa_str, age_prior=age_prior)


post_OI = Clade(
    name='Post_OI',
    members=['Irish_B', 'Scots_Gaelic']
)
old_irish_irish_scots = Clade(
    name='Old_Irish_Irish_Scots',
    members=['Old_Irish', post_OI]
)
brythonic = Clade(
    name='Brythonic',
    members=['Welsh_N', 'Breton_ST', 'Cornish'],
    age_prior='<distr id="Brythonic.tmrca" spec="Normal" mean="1550" sigma="300"/>'
)
celtic = Clade(
    name='Celtic',
    members=[old_irish_irish_scots, brythonic],
    age_prior='<distr id="Celtic.tmrca" spec="beast.math.distributions.Uniform" lower="1700" upper="5000.0"/>'
)
french_iberian = Clade(
    name='FrenchIberian',
    members=['Provencal', 'French', 'Walloon', 'Spanish', 'Portuguese_ST', 'Catalan'],
    age_prior='<distr id="FrenchIberian.tmrca" spec="Normal" mean="1400" sigma="100"/>'
)
post_latin = Clade(
    name='Post_Latin',
    members=[french_iberian, 'Sardinian_N', 'Sardinian_C', 'Rumanian_List', 'Vlach', 'Italian', 'Ladin', 'Friulian', 'Romansh']
)
latin_romance = Clade(
    name='Latin_Romance',
    members=['Latin', post_latin],
    age_prior='<distr id="Latin_Romance.tmrca" spec="Normal" mean="2100" sigma="200"/>'
)
italic = Clade(
    name='Italic',
    members=['Umbrian', 'Oscan', latin_romance]
)
post_ON = Clade(
    name='Post_ON',
    members=['Icelandic_ST', 'Faroese', 'Norwegian']
)
old_norse_icelandic_faroese = Clade(
    name='Old_Norse_Icelandic_Faroese',
    members=['Old_Norse', post_ON]
)
north_germanic = Clade(
    name='NorthGermanic',
    members=['Swedish_Up', 'Swedish_Vl', 'Swedish', 'Danish', old_norse_icelandic_faroese]
)
post_OHG = Clade(
    name='Post_OHG',
    members=['German', 'Luxembourgish', 'Schwyzerdutsch']
)
high_german = Clade(
    name='High_German',
    members=['Old_High_German', post_OHG]
)
english_clade = Clade(
    name='English_Clade',
    members=['Old_English', 'English']
)
west_germanic = Clade(
    name='WestGermanic',
    members=['Dutch_List', 'Flemish', 'Frisian', 'Afrikaans', english_clade, high_german],
    age_prior='<distr id="WestGermanic.tmrca" spec="Normal" mean="1550" sigma="30"/>'
)
north_west_germanic = Clade(
    name='NorthWestGermanic',
    members=[north_germanic, west_germanic],
    age_prior='<distr id="NorthWestGermanic.tmrca" spec="Normal" mean="1875" sigma="67"/>'
)
germanic = Clade(
    name='Germanic',
    members=['Gothic', north_west_germanic]
)
# lithuanian_latvian = Clade(
#     name='LithuanianLatvian',
#     members=['Lithuanian_ST', 'Latvian'],
#     age_prior='<distr id="LithuanianLatvian.tmrca" spec="Normal" mean="1350" sigma="30"/>'
# )
# baltic = Clade(
#     name='Baltic',
#     members=['Old_Prussian', 'Lithuanian_ST', 'Latvian']
# )
# south_slavic = Clade(
#     name='South_Slavic',
#     members=['Macedonian', 'Bulgarian', 'Serbocroatian', 'Old_Church_Slavonic']
# )
# slavic = Clade(
#     name='Slavic',
#     members=['Slovenian', 'Lower_Sorbian', 'Upper_Sorbian', 'Czech', 'Slovak', 'Czech_E', 'Ukrainian', 'Byelorussian', 'Polish', 'Russian', 'Macedonian', 'Bulgarian', 'Serbocroatian', 'Old_Church_Slavonic'],
#     age_prior='<distr id="Slavic.tmrca" spec="LogNormalDistributionModel" offset="1200.0" M="300" S="0.6" meanInRealSpace="true"/>'
# )
# balto_slavic = Clade(
#     name='Balto-Slavic',
#     members=['Old_Prussian', 'Lithuanian_ST', 'Latvian', 'Slovenian', 'Lower_Sorbian', 'Upper_Sorbian', 'Czech', 'Slovak', 'Czech_E', 'Ukrainian', 'Byelorussian', 'Polish', 'Russian', 'Macedonian', 'Bulgarian', 'Serbocroatian', 'Old_Church_Slavonic'],
#     age_prior='<distr id="Balto-Slavic.tmrca" spec="beast.math.distributions.Uniform" lower="2200" upper="3400"/>'
# )
# modern_greek = Clade(
#     name='Modern_Greek',
#     members=['Greek_Mod', 'Greek_Ml']
# )
# greek = Clade(
#     name='Greek',
#     members=['Ancient_Greek', 'Greek_Mod', 'Greek_Ml']
# )
# eastern_iranian = Clade(
#     name='Eastern_Iranian',
#     members=['Wakhi', 'Avestan', 'Digor_Ossetic', 'Iron_Ossetic']
# )
# albanian = Clade(
#     name='Albanian',
#     members=['Albanian_Top', 'Albanian_G', 'Albanian_K', 'Albanian_C']
# )

CLADES = [
    post_OI,
    old_irish_irish_scots,
    brythonic,
    celtic,
    french_iberian,
    post_latin,
    latin_romance,
    italic,
    post_ON,
    old_norse_icelandic_faroese,
    north_germanic,
    post_OHG,
    high_german,
    english_clade,
    west_germanic,
    north_west_germanic,
    germanic,
]


def read_dataset(csv_path):
    """Read the cognate data-set from a TSV file and run some basic parsing."""
    ielex = pd.read_csv(csv_path, sep='\t', dtype=str)
    ielex['concept'] = ielex['cc_alias'].map(lambda s: s.split('-')[0].strip())
    return ielex


def drop_noncoded(ielex):
    """Drop rows which do not have a cognate ID."""
    return ielex.loc[~pd.isna(ielex.cc_id)]


def parse_loanwords(ielex: pd.DataFrame):
    loans = pd.DataFrame(data=0, dtype=int,
                         columns=ielex.language.unique(),
                         index=ielex.concept.unique())

    for concept, concepts_grouped in ielex.groupby('concept'):
        concepts_grouped.set_index('language')
        for i_cc, (_, cc_grouped) in enumerate(concepts_grouped.groupby('cc_id')):
            for _, row in cc_grouped.iterrows():
                if row.status in ('LOAN', 'LOAN,EXCLUDE'):
                    loans.loc[concept, row.language] = 1
    return loans


def parse_data_matrix(ielex: pd.DataFrame, exclude_loans=False):
    """Parse the pandas dataframe, containing the IELEX data into nested dictionary, which
    is easier to use for compiling the BEAST XML.

    Args:
        ielex (pd.DataFrame): the Pandas DataFrame containing the raw IELex dataset.
        exclude_loans (boolean): Whether to exclude loan words from the analysis.

    Returns:
        dict: The parsed data set which is nested dictionary:
                {language (str) ->  {concept (str) -> data (string of 0/1/?)}}
    """
    languages = ielex.language.unique()
    data = {lang: {} for lang in languages}
    for concept, concepts_grouped in ielex.groupby('concept'):
        concepts_grouped.set_index('language')
        n_cc = len(concepts_grouped.cc_id.unique())
        concept_data = defaultdict(list)
        for i_cc, (_, cc_grouped) in enumerate(concepts_grouped.groupby('cc_id')):
            for _, row in cc_grouped.iterrows():
                if row.status in ('EXCLUDE', 'LOAN,EXCLUDE'):
                    continue
                if row.status == 'LOAN':
                    if exclude_loans:
                        continue
                concept_data[row.language].append(i_cc)

        for lang in languages:
            if lang in concept_data:
                data[lang][concept] = ''.join(['1' if i_cc in concept_data[lang] else '0'
                                               for i_cc in range(n_cc)])
            else:
                data[lang][concept] = ''.join(['?' for _ in range(n_cc)])

    return data


def encode_binary_array(a):
    """Compact bit-string encoding for binary arrays."""
    return ''.join(['%i' % x for x in a])


def compile_ielex_xml(data_raw: dict, ascertainment_correction=True, min_coverage=0.0,
                      fixed_topolgy=False, fixed_node_heights=False,
                      use_contactrees=True, expected_conversions=0.15,
                      exclude_loans=False):
    """Compile the IELEX data together with hardcoded settings into a BEAST XML.

    Args:
        data_raw (pd.DataFrame):
        ascertainment_correction (bool):
        min_coverage (float):
        use_contactrees (bool):
        fixed_topolgy (bool):
        fixed_node_heights (bool):

    Returns:
        str: The compiled BEAST2 XML file.
    """
    data = parse_data_matrix(data_raw, exclude_loans=exclude_loans)

    language_alignments = defaultdict(str)
    concept_ranges = {}
    concept_n_sites = {}
    coverage = defaultdict(int)

    n_concepts = 0
    for lang, data_l in data.items():
        i = 0
        for concept, data_l_c in data_l.items():
            is_nan = data_l_c.startswith('?')
            if not is_nan:
                coverage[lang] += 1

            if ascertainment_correction:
                if is_nan:
                    data_l_c = '?' + data_l_c
                else:
                    data_l_c = '0' + data_l_c

            language_alignments[lang] += data_l_c

            if concept not in concept_ranges:
                i_next = i + len(data_l_c)
                concept_ranges[concept] = (i, i_next-1)
                concept_n_sites[concept] = len(data_l_c)
                i = i_next
                n_concepts += 1


    # site_counts = [*concept_n_sites.values()]
    # plt.bar(*np.unique(site_counts, return_counts=True))
    # plt.show()
    # exit()

    # Filter languages with insufficient data
    for lang in list(language_alignments.keys()):
        if lang not in INCLUDED_LANGUAGES:
            # print(lang.ljust(20) + '%.1f%%' % (100 * coverage[lang] / n_concepts))
            language_alignments.pop(lang)
        # else:
        #     print(' '*30 + lang.ljust(20) + '%.1f%%' % (100 * coverage[lang] / n_concepts))

    # Fill in alignments for each language
    alignments = ''
    for lang, alignment in language_alignments.items():
        if coverage[lang] < min_coverage * n_concepts:
            continue
        alignments += xml_snippets.ALIGNMENT.format(tax_name=lang, data=alignment)

    # Fill in filtered alignments for each block
    filtered_alignments = ''
    for concept, (start, end) in concept_ranges.items():
        filtered_alignments += xml_snippets.FILTERED_ALIGNMENT.format(
            concept=concept,
            start=start+1,
            end=end+1,
            excludeto=int(ascertainment_correction)
        )

    # Compile list of tip-dates
    format_date = lambda l: f'          {l} = {TIP_DATES[l]}'
    dates = ',\n'.join(map(format_date, language_alignments.keys()))

    # Compile MCRA priors
    mrca_priors = ''.join(clade.to_xml() for clade in CLADES)

    # Compile likelihoods
    # with slow or fast clock depending on the number of sites per concept
    likelihoods = ''
    for concept, n_sites in concept_n_sites.items():
        # TODO: maybe reconsider the bins
        # With 3 bins I would do slow 1...5, medium 6...9, fast 10...\inf
        slow_or_fast = 'slow' if n_sites <= 6 else 'fast'
        if use_contactrees:
            likelihoods += xml_snippets.CONTACTREES_LIKELIHOOD.format(
                concept=concept,
                slow_or_fast=slow_or_fast
            )
        else:
            likelihoods += xml_snippets.BASICTREES_LIKELIHOOD.format(
                concept=concept,
                slow_or_fast=slow_or_fast
            )
            expected_conversions = 0.0

    # Compile operators
    operators = '\n'
    if fixed_node_heights:
        assert fixed_topolgy
    if not fixed_topolgy:
        operators += xml_snippets.TOPOLOGY_OPERATORS
    if not fixed_node_heights:
        operators += xml_snippets.NODE_HEIGHTS_OPERATORS
    if use_contactrees:
        operators += xml_snippets.CONTACT_OPERATORS

    # Prepare the starting tree
    starting_tree = fix_tree(STARTING_TREE)

    # Add word-tree loggers if required
    if use_contactrees:
        # word_tree_loggers = xml_snippets.WORD_TREE_LOGGERS
        word_tree_loggers = ''
    else:
        word_tree_loggers = ''

    # Load the BEAST2 XML template
    with open(XML_TEMPLATE_PATH, 'r') as xml_template_file:
        xml_template = xml_template_file.read()

    # Put everything together...
    xml_filled = xml_template.format(
        concepts=','.join(concept_ranges.keys()),
        languages=','.join(language_alignments.keys()),
        alignments=alignments,
        filtered_alignments=filtered_alignments,
        concept_clocks=...,
        dates=dates,
        mrca_priors=mrca_priors,
        likelihood=likelihoods,
        operators=operators,
        starting_tree=starting_tree,
        expected_conversions=str(expected_conversions),
        word_tree_loggers=word_tree_loggers
    )

    return xml_filled


def fix_tree(old_newick):
    old_newick_no_attr = drop_attributes(old_newick)

    #####################################
    #
    # tree: Tree = Tree.from_newick(old_newick_no_attr)
    # tree.rename_nodes(translate)
    # tree.rename_nodes(RENAME_CHANG)
    #
    # for name in INCLUDED_LANGUAGES:
    #     assert name in [n.name for n in tree.iter_descendants()], name
    # tree.remove_nodes_except_names(INCLUDED_LANGUAGES + [''])
    #
    # for node in tree.iter_leafs():
    #     err = node.get_age() - TIP_DATES[node.name]
    #
    #     node._length += err
    #     if node.length < 0:
    #         node.parent._length += node.length - 1.0
    #         if node.sibling:
    #             node.sibling._length -= node.length - 1.0
    #         node._length = 1.0
    #
    # new_newick = tree.to_newick()
    #
    #########################

    tree: Node = parse_node(old_newick_no_attr.strip(' ;'))
    translate_node_names(tree, translate)
    translate_node_names(tree, RENAME_CHANG)

    for name in INCLUDED_LANGUAGES:
        assert name in tree.get_leaf_names(), name

    tree.prune_by_names(INCLUDED_LANGUAGES, inverse=True)
    tree.remove_redundant_nodes()
    tree.length = 0.0
    tree = parse_node(tree.newick)

    for node in tree.get_leaves():
        err = get_age(node) - TIP_DATES[node.name]

        node.length += err
        if node.length < 0:
            node.ancestor.length += node.length - 1.0
            sibling = get_sibling(node)
            if sibling is not None:
                sibling.length -= node.length - 1.0
            node.length = 1.0

    return tree.newick


def drop_attributes(newick):
    while '[' in newick:
        before_attr, _, rest = newick.partition('[')
        attr_str, _, after_attr = rest.partition(']')
        newick = before_attr + after_attr
    return newick


def filter_languages(df):
    include = df.language.isin(INCLUDED_LANGUAGES)
    return df[include]


if __name__ == '__main__':
    DATA_PATH = Path('resources/ielex-130421-ag-cc.tsv')
    ielex_df = read_dataset(DATA_PATH)
    ielex_df = drop_noncoded(ielex_df)
    print(ielex_df.columns.to_numpy())
    # for lang in ielex_df.language.unique():
    #     # if lang in [*CELTIC, *GERMANIC, *ROMANCE]:
    #     if lang.lower().startswith('irish'):
    #         print(str.ljust(lang, 20), ielex_df.concept[ielex_df.language == lang].to_numpy())

    # A = ielex_df[ielex_df.language == 'English']
    # B = ielex_df[ielex_df.language == 'German']
    # print(np.mean([a in B.cc_id.to_list() for a in A.cc_id.to_list()]))
    # # print([a in B for a in A])
    # for i, a in A.iterrows():
    #     if a['cc_id'] in B.cc_id.to_list(): continue
    #     concept = a['concept']
    #     a_lexeme = a['lexeme']
    #     b_lexemes = B[B.concept == concept].lexeme.to_list()
    #     print(concept.ljust(15), a_lexeme.ljust(15), b_lexemes)
    #
    # exit()
    ielex_df = filter_languages(ielex_df)
    ielex_df.to_csv('resources/ielex-subset.tsv', sep='\t', index=False)

    # for l in INCLUDED_LANGUAGES:
    #     if l not in TIP_DATES:
    #         print(l)
    #
    # print(len(INCLUDED_LANGUAGES))
    # exit()

    N_RUNS = 2

    RUN_CONFIGURATIONS = {
        'CT_full': {'use_contactrees': True,
                    'fixed_topolgy': False,
                    'fixed_node_heights': False,
                    'exclude_loans': False},

        'BT_full': {'use_contactrees': False,
                    'fixed_topolgy': False,
                    'fixed_node_heights': False,
                    'exclude_loans': False},

        'CT_fixedTopology': {'use_contactrees': True,
                              'fixed_topolgy': True,
                              'fixed_node_heights': False,
                              'exclude_loans': False},

        'BT_fixedTopology': {'use_contactrees': False,
                              'fixed_topolgy': True,
                              'fixed_node_heights': False,
                              'exclude_loans': False},

        'CT_fixedTree': {'use_contactrees': True,
                          'fixed_topolgy': True,
                          'fixed_node_heights': True,
                          'exclude_loans': False},

        'BT_fixedTree': {'use_contactrees': False,
                          'fixed_topolgy': True,
                          'fixed_node_heights': True,
                          'exclude_loans': False},

        # basicTrees excluding loans:

        'BT_full_noLoans': {'use_contactrees': False,
                            'fixed_topolgy': False,
                            'fixed_node_heights': False,
                            'exclude_loans': True},

        'BT_fixedTopology_noLoans': {'use_contactrees': False,
                                     'fixed_topolgy': True,
                                     'fixed_node_heights': False,
                                     'exclude_loans': True},

        'BT_fixedTree_noLoans': {'use_contactrees': False,
                                 'fixed_topolgy': False,
                                 'fixed_node_heights': False,
                                 'exclude_loans': True},

    }

    for run_name, kwargs in RUN_CONFIGURATIONS.items():
        # Create folder structure
        run_directory = Path(f'runs/{run_name}')
        if not os.path.exists(run_directory):
            os.mkdir(run_directory)

        # Compile the BEAST-XML from the data and the run-configuration
        ie_xml_str = compile_ielex_xml(ielex_df, **kwargs)

        # Write the XML to ´N_RUNS different files´
        for i in range(1, N_RUNS+1):
            # Create folder structure
            # if kwargs['use_contactrees']:
            #     if not os.path.exists(run_directory / f'wordtrees_run{i}/'):
            #         os.mkdir(run_directory / f'wordtrees_run{i}/')

            with open(run_directory / f'ie_run{i}.xml', 'w') as ie_xml_file:
                ie_xml_file.write(ie_xml_str)

        # Create a shell script to start all runs in parallel
        with open(run_directory / 'start_runs.sh', 'w') as ie_xml_file:
            lines = []
            for i in range(1, N_RUNS + 1):
                lines.append(f'beast -threads 4 -overwrite ie_run{i}.xml &')
            ie_xml_file.write('\n'.join(lines))

    # loans = parse_loanwords(ielex_df)
    # loans.to_csv(LOANS_CSV_PATH)
