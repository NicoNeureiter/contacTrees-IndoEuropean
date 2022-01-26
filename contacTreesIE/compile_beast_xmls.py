#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
from collections import defaultdict
from pathlib import Path
from typing import List, Union
from copy import deepcopy
from enum import Enum

import pandas as pd
from newick import Node, parse_node

from contacTreesIE.newick_util import get_sibling, translate_node_names, get_age, get_node_by_name
from contacTreesIE.preprocessing.language_lists import *
from contacTreesIE.preprocessing import xml_snippets
from contacTreesIE.preprocessing.xml_snippets import Samplers
from contacTreesIE.preprocessing.starting_trees import CHANG_MEDIUM_TREE
from contacTreesIE.preprocessing.starting_trees import CHANG_MEDIUM_TRANSLATE

XML_TEMPLATE_PATH = 'resources/ie_template.xml'
LOANS_CSV_PATH = 'loans.csv'
ZOMBIE_LATIN = 'Latin_preserved'
MEDIEVAL_LATIN = 'Latin_M'

DATASET = CHANG_MEDIUM
# DATASET = TINY_SET
INCLUDED_CLADES = CELTIC + GERMANIC + ROMANCE
INCLUDED_LANGUAGES = [l for l in DATASET if l in INCLUDED_CLADES] + [MEDIEVAL_LATIN]
RENAME = RENAME_CHANG
translate = CHANG_MEDIUM_TRANSLATE
STARTING_TREE = CHANG_MEDIUM_TREE

SA_PRIOR_SIGMA = 100.0

TIP_DATES = defaultdict(float)
TIP_DATES.update({
    "Old_Irish": 1250.0,
    "Cornish": 300.0,
    "Umbrian": 2200.0,
    "Oscan": 2200.0,
    "Latin": 2050.0,
    MEDIEVAL_LATIN: 1000.0,
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


class CalibrationModes(Enum):
    BOUCKAERT = 1
    CHANG = 2


CALIBRATION_MODE = CalibrationModes.BOUCKAERT


class Clade(object):

    clade_template = '''\
        <distribution id="{clade}.prior" spec="MRCAPrior" tree="@acg" monophyletic="true">
          <taxonset id="{clade}" spec="TaxonSet">
{taxa}
          </taxonset>
{age_prior}        </distribution>\n\n'''

    taxon_template = '            <taxon idref="{taxon}" />'

    def __init__(self,
                 name: str,
                 members: list,
                 age_prior: Union[str, dict] = ''):
        self.name = name
        self.members = [m for m in members
                        if isinstance(m, Clade) or (m in INCLUDED_LANGUAGES)]
        # self.members = members
        if isinstance(age_prior, dict):
            self.age_prior = age_prior[CALIBRATION_MODE]
        else:
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
    members=['Irish_B', 'Scots_Gaelic'],
    age_prior={
        CalibrationModes.CHANG: '<distr id="Post_OI.tmrca" spec="beast.math.distributions.Uniform" lower="1050" upper="10000"/>',
        CalibrationModes.BOUCKAERT: '',
    }
)
old_irish_irish_scots = Clade(
    name='Old_Irish_Irish_Scots',
    members=['Old_Irish', post_OI],
    age_prior=f'<distr id="OldIrish_SA.tmrca" spec="Normal" offset="1250" sigma="{SA_PRIOR_SIGMA}"/>'
)
brythonic = Clade(
    name='Brythonic',
    members=['Welsh_N', 'Breton_ST', 'Cornish'],
    age_prior={
        CalibrationModes.CHANG: '<distr id="Brythonic.tmrca" spec="beast.math.distributions.Uniform" lower="1250" upper="10000"/>',
        CalibrationModes.BOUCKAERT: '<distr id="Brythonic.tmrca" spec="Normal" mean="1550" sigma="100"/>',
    }
)
celtic = Clade(
    name='Celtic',
    members=[old_irish_irish_scots, brythonic],
    age_prior={
        CalibrationModes.CHANG: '',
        CalibrationModes.BOUCKAERT: '<distr id="Celtic.tmrca" spec="beast.math.distributions.Uniform" lower="1700" upper="10000.0"/>',
    }
)
french_iberian = Clade(
    name='FrenchIberian',
    members=['Provencal', 'French', 'Walloon', 'Spanish', 'Portuguese_ST', 'Catalan'],
    age_prior={
        CalibrationModes.CHANG: '',
        CalibrationModes.BOUCKAERT: '<distr id="FrenchIberian.tmrca" spec="Normal" mean="1400" sigma="100"/>',
    }
)

post_latin = Clade(
    name='Post_Latin',
    members=[french_iberian, 'Sardinian_N', 'Sardinian_C', 'Rumanian_List', 'Vlach', 'Italian', 'Ladin', 'Friulian', 'Romansh']
)

latin_romance = Clade(
    name='Latin_Romance',
    members=['Latin', post_latin],
    age_prior={
        CalibrationModes.CHANG: '<distr id="Latin_Romance.tmrca" spec="beast.math.distributions.Uniform" lower="1750" upper="10000"/>',
        CalibrationModes.BOUCKAERT: '<distr id="Latin_Romance.tmrca" spec="Normal" mean="2050" sigma="150"/>',
    }
)
"""COMMENT: For the tMRCA of Romance the lower bound is given by the sample date of Latin.
We add an informative prior as an upper bound to ensure that the branch above latin is not
 too long (we know it was at least close to a sampled ancestor of romance languages)."""

# italic = Clade(
#     name='Italic',
#     members=['Umbrian', 'Oscan', latin_romance]
# )
post_ON = Clade(
    name='Post_ON',
    members=['Icelandic_ST', 'Faroese', 'Norwegian']
)
old_norse_icelandic_faroese = Clade(
    name='Old_Norse_Icelandic_Faroese',
    members=['Old_Norse', post_ON],
    age_prior=f'<distr id="OldNorse_SA.tmrca" spec="Normal" offset="775" sigma="{SA_PRIOR_SIGMA}"/>'
)
north_germanic = Clade(
    name='NorthGermanic',
    members=['Swedish_Up', 'Swedish_Vl', 'Swedish', 'Danish', old_norse_icelandic_faroese],
    age_prior={
        CalibrationModes.CHANG: '<distr id="NorthGermanic.tmrca" spec="beast.math.distributions.Uniform" lower="1500" upper="10000"/>',
        CalibrationModes.BOUCKAERT: '',
    },
)
post_OHG = Clade(
    name='Post_OHG',
    members=['German', 'Luxembourgish', 'Schwyzerdutsch']
)
high_german = Clade(
    name='High_German',
    members=['Old_High_German', post_OHG],
    age_prior=f'<distr id="High_German_SA.tmrca" spec="Normal" mean="1050" sigma="{SA_PRIOR_SIGMA}"/>',
)
english_clade = Clade(
    name='English_Clade',
    members=['Old_English', 'English'],
    age_prior=f'<distr id="English_SA.tmrca" spec="Normal" mean="1000" sigma="{SA_PRIOR_SIGMA}"/>'
)
west_germanic = Clade(
    name='WestGermanic',
    members=['Dutch_List', 'Flemish', 'Frisian', 'Afrikaans', english_clade, high_german],
    age_prior={
        CalibrationModes.CHANG: '',
        CalibrationModes.BOUCKAERT: '<distr id="WestGermanic.tmrca" spec="Normal" mean="1550" sigma="100"/>',
    },
)
north_west_germanic = Clade(
    name='NorthWestGermanic',
    members=[north_germanic, west_germanic],
    age_prior={
        CalibrationModes.CHANG: '',
        CalibrationModes.BOUCKAERT: '<distr id="NorthWestGermanic.tmrca" spec="Normal" mean="1875" sigma="150"/>',
    },
)
germanic = Clade(
    name='Germanic',
    members=['Gothic', north_west_germanic],
    age_prior={
        CalibrationModes.CHANG: '<distr id="Germanic.tmrca" spec="beast.math.distributions.Uniform" lower="2250" upper="10000"/>',
        CalibrationModes.BOUCKAERT: '',
    },
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
    # italic,
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
    try:
        ielex['concept'] = ielex['cc_alias'].map(lambda s: s.split('-')[0].strip())
    except AttributeError as e:
        for cc_alias in ielex.cc_alias:
            if not isinstance(cc_alias, str) or ('-' not in cc_alias):
                print('Invalid cc_alias:', cc_alias)
        raise e

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
                if row.status in ('LOAN',):  #, 'LOAN,EXCLUDE'):
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

    # ´data´ is a nested dict:
    # {language: {concept: "absence/presence of each cognate for this concept as a binary string"}}
    data = {lang: {} for lang in languages}

    for concept, concepts_grouped in ielex.groupby('concept'):
        concepts_grouped.set_index('language')
        n_cc = len(concepts_grouped.cc_id.unique())
        concept_data = defaultdict(list)

        # Collect all cognates for the current concept
        for i_cc, (_, cc_grouped) in enumerate(concepts_grouped.groupby('cc_id')):
            for _, row in cc_grouped.iterrows():
                if row.status in ('EXCLUDE', 'LOAN,EXCLUDE', 'WRONG'):
                    continue
                elif row.status == 'LOAN':
                    if exclude_loans:
                        continue
                elif not pd.isna(row.status):
                    raise ValueError(f'Unknown status "{row.status}" in {row.language}:{row.cc_alias}')
                concept_data[row.language].append(i_cc)

        # Assemble cognate absence/presence for the current concept and each language in a binary string
        for lang in languages:
            if lang in concept_data:
                data[lang][concept] = ''.join(['1' if i_cc in concept_data[lang] else '0'
                                               for i_cc in range(n_cc)])
                assert '1' in data[lang][concept]
            else:
                data[lang][concept] = ''.join(['?' for _ in range(n_cc)])

    return data


def encode_binary_array(a):
    """Compact bit-string encoding for binary arrays."""
    return ''.join(['%i' % x for x in a])


def compile_ielex_xml(data_raw: dict, ascertainment_correction=True, min_coverage=0.0,
                      fixed_topolgy=False, fixed_node_heights=False,
                      use_contactrees=True, expected_conversions=0.25,
                      exclude_loans=False, add_zombie_latin=True,
                      add_medieval_latin=True,
                      use_covarion=True, sample_acg_prior_params=False,
                      clock_stdev_prior=0.04, fix_clock_stdev=True,
                      sampler=Samplers.MCMC, chain_length=20000000):
    """Compile the IELEX data together with hardcoded settings into a BEAST XML.

    Args:
        data_raw (pd.DataFrame): The pandas data-frame containing the raw IELex data.
        ascertainment_correction (bool): whether to apply ascertainment correction.
        min_coverage (float): the minimum fraction of cognates a language need to be included.
        use_contactrees (bool): whether use contacTrees for the reconstruction.
        fixed_topolgy (bool): whether to fix the tree topology.
        fixed_node_heights (bool): whether to fix the height of internal nodes in the tree.
        expected_conversions (float): The expected number of conversions in the prior.
        exclude_loans (bool): whether to exclude loans from the data or not.
        add_zombie_latin (bool): whether to add a preserved copy of the Latin taxon as
            a contemporary language to allow for recent borrowing of Latin words.

    Returns:
        str: The compiled BEAST2 XML file.
    """
    clades = deepcopy(CLADES)
    data = parse_data_matrix(data_raw, exclude_loans=exclude_loans)

    if add_medieval_latin:
        INCLUDED_LANGUAGES.append(MEDIEVAL_LATIN)

        for clade in clades:
            if 'Latin' in clade.members:
                clade.members.append(MEDIEVAL_LATIN)

        clades.append(
            Clade(
                name='Latin_descendants',
                members=['Latin', MEDIEVAL_LATIN],
                age_prior='<distr id="Latin_SA.tmrca" spec="Normal" offset="2050" sigma="20"/>',
            )
        )

    if add_zombie_latin:
        INCLUDED_LANGUAGES.append(ZOMBIE_LATIN)

        # Use the most recent form of latin as data for zombie latin
        if add_medieval_latin:
            data[ZOMBIE_LATIN] = data[MEDIEVAL_LATIN]
        else:
            data[ZOMBIE_LATIN] = data['Latin']

        for clade in clades:
            if 'Latin' in clade.members:
                clade.members.append(ZOMBIE_LATIN)

    if add_medieval_latin and add_zombie_latin:
        clades.append(
            Clade(
                name='Latin_M_descendants',
                members=[MEDIEVAL_LATIN, ZOMBIE_LATIN],
                age_prior='<distr id="Latin_M_SA.tmrca" spec="Normal" offset="1000" sigma="10"/>',
            )
        )



    language_alignments = defaultdict(str)
    concept_ranges = {}
    concept_n_sites = {}
    coverage = defaultdict(int)

    n_concepts = 0
    words_per_language = defaultdict(int)
    for lang, data_l in data.items():
        i = 0
        for concept, data_l_c in data_l.items():
            is_nan = data_l_c.startswith('?')
            if not is_nan:
                coverage[lang] += 1
                words_per_language[lang] += sum(map(int, data_l_c))
            elif lang == 'English':
                print(concept)

            # if is_nan and lang == 'Afrikaans':
            #     print(concept)

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

    for c, lang in sorted([(c, l) for l, c in coverage.items()]):
        print(lang.ljust(15, ' '), str(c).ljust(5, ' '), words_per_language[lang], ' ', words_per_language[lang]-c)
    # exit()

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
            print('Coverage too low:', lang, coverage[lang])
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
    format_date = lambda l: f'        {l} = {TIP_DATES[l]}'
    dates = ',\n'.join(map(format_date, language_alignments.keys()))

    # Compile MCRA priors
    mrca_priors = ''.join(clade.to_xml() for clade in clades)

    # Compile likelihoods
    # with slow or fast clock depending on the number of sites per concept
    likelihoods = ''
    for concept, n_sites in concept_n_sites.items():
        # TODO: maybe reconsider the bins
        # With 3 bins I would do slow 1...5, medium 6...9, fast 10...\inf
        if n_sites <= 5:
            site_category = 'slow'
        elif n_sites <= 9:
            site_category = 'medium'
        else:
            site_category = 'fast'

        if use_contactrees:
            likelihoods += xml_snippets.CONTACTREES_LIKELIHOOD.format(
                concept=concept,
                site_cat=site_category
            )
        else:
            likelihoods += xml_snippets.BASICTREES_LIKELIHOOD.format(
                concept=concept,
                site_cat=site_category
            )
            expected_conversions = 0.0

    frozen_taxa = xml_snippets.FROZEN_TAXA if add_zombie_latin else ''

    if use_covarion:
        substitution_model = xml_snippets.COVARION_MODEL
        substitution_model_prior = xml_snippets.COVARION_PRIORS
        data_type = xml_snippets.COVARION_DATA_TYPE
    else:
        substitution_model = xml_snippets.CTMC_MODEL
        substitution_model_prior = xml_snippets.CTMC_PRIORS
        data_type = xml_snippets.CTMC_DATA_TYPE

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
    starting_tree = fix_tree(STARTING_TREE, add_zombie_latin, add_medieval_latin)

    # Add word-tree loggers if required
    if use_contactrees:
        # word_tree_loggers = xml_snippets.WORD_TREE_LOGGERS
        word_tree_loggers = ''
    else:
        word_tree_loggers = ''

    # Prepare the `run` tag
    run_tag = xml_snippets.SAMPLER_TAG[sampler].format(chain_length=chain_length)

    # Load the BEAST2 XML template
    with open(XML_TEMPLATE_PATH, 'r') as xml_template_file:
        xml_template = xml_template_file.read()

    if fix_clock_stdev:
        clock_stdev_operator = ''
    else:
        clock_stdev_operator = xml_snippets.CLOCK_STDEV_OPERATOR

    # Put everything together...
    xml_filled = xml_template.format(
        concepts=','.join(concept_ranges.keys()),
        languages=','.join(language_alignments.keys()),
        alignments=alignments,
        filtered_alignments=filtered_alignments,
        dates=dates,
        mrca_priors=mrca_priors,
        likelihood=likelihoods,
        clock_rate_prior=xml_snippets.CLOCK_RATE_PRIOR_FLAT,
        clock_stdev_prior=clock_stdev_prior,
        subst_model=substitution_model,
        subst_model_prior=substitution_model_prior,
        data_type=data_type,
        operators=operators,
        clock_stdev_operator=clock_stdev_operator,
        starting_tree=starting_tree,
        expected_conversions=str(expected_conversions),
        frozen_taxa=frozen_taxa,
        word_tree_loggers=word_tree_loggers,
        run_tag=run_tag,
    )

    return xml_filled


def fix_tree(old_newick: str,
             add_zombie_latin: bool,
             add_medieval_latin: bool) -> str:
    old_newick_no_attr = drop_attributes(old_newick)

    tree = parse_node(old_newick_no_attr.strip(' ;'))
    translate_node_names(tree, translate)
    translate_node_names(tree, RENAME_CHANG)

    for name in INCLUDED_LANGUAGES:
        if name in [ZOMBIE_LATIN, MEDIEVAL_LATIN]: continue
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

    if add_medieval_latin:
        latin: Node = get_node_by_name(tree, 'Latin')
        old_parent: Node = latin.ancestor
        medieval_latin: Node = Node.create(
            name=MEDIEVAL_LATIN,
            length='%.8f' % (get_age(latin) + 1.0 - TIP_DATES[MEDIEVAL_LATIN])
        )
        new_parent: Node = Node.create(
            name='',
            length='%.8f' % (latin.length - 1.0),
            descendants=[latin, medieval_latin]
        )
        latin.length = 1.0
        old_parent.descendants.remove(latin)
        old_parent.add_descendant(new_parent)


    if add_zombie_latin:
        if add_medieval_latin:
            parent = get_node_by_name(tree, MEDIEVAL_LATIN)
        else:
            parent = get_node_by_name(tree, 'Latin')

        old_parent: Node = parent.ancestor
        zombie_latin: Node = Node.create(
            name=ZOMBIE_LATIN,
            length='%.8f' % (get_age(parent) + 1.0)
        )
        new_parent: Node = Node.create(
            name='',
            length='%.8f' % (parent.length - 1.0),
            descendants=[parent, zombie_latin]
        )
        parent.length = 1.0
        old_parent.descendants.remove(parent)
        old_parent.add_descendant(new_parent)

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
    # DATA_PATH = Path('resources/ielex-130421-ag-cc.tsv')
    DATA_PATH = Path('resources/data-mittellatein-2021-09-30.csv')
    ielex_df = read_dataset(DATA_PATH)
    ielex_df = drop_noncoded(ielex_df)
    print(ielex_df.columns.to_numpy())
    print(ielex_df.language.unique())
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
    path_base, _, path_ext = str(DATA_PATH).rpartition('.')
    subset_path = path_base + '-subset.tsv'
    ielex_df = filter_languages(ielex_df)
    ielex_df.to_csv(subset_path, sep='\t', index=False)

    # for l in INCLUDED_LANGUAGES:
    #     if l not in TIP_DATES:
    #         print(l)
    #
    # print(len(INCLUDED_LANGUAGES))
    # exit()

    N_RUNS = 5

    RUN_CONFIGURATIONS = {

        'BT_fixTopo/covarion': {
            'sampler': Samplers.MCMC,
            'chain_length': 20000000,
            'use_contactrees': False,
            'fixed_topolgy': True,
            'fixed_node_heights': False,
            'exclude_loans': False,
            'add_zombie_latin': True,
            'use_covarion': True,
        },

        'BT_fixTopo/covarion_noLoans': {
            'sampler': Samplers.MCMC,
            'chain_length': 20000000,
            'use_contactrees': False,
            'fixed_topolgy': True,
            'fixed_node_heights': False,
            'exclude_loans': True,
            'add_zombie_latin': True,
            'use_covarion': True,
        },

        'BT_full': {
            'sampler': Samplers.MCMC,
            'chain_length': 20000000,
            'use_contactrees': False,
            'fixed_topolgy': False,
            'fixed_node_heights': False,
            'exclude_loans': False,
            'add_zombie_latin': True,
            'use_covarion': True
        },

        'CT_fixTopo/covarion': {
            'sampler': Samplers.MC3,
            'chain_length': 20000000,
            'use_contactrees': True,
            'fixed_topolgy': True,
            'fixed_node_heights': False,
            'exclude_loans': False,
            'add_zombie_latin': True,
            'use_covarion': True,
        },

        'CT_full': {
            'sampler': Samplers.MC3,
            'chain_length': 25000000,
            'use_contactrees': True,
            'fixed_topolgy': False,
            'fixed_node_heights': False,
            'exclude_loans': False,
            'add_zombie_latin': True,
            'use_covarion': True
        },

    }

    for run_name, kwargs in RUN_CONFIGURATIONS.items():

        # Create folder structure
        run_directory = Path(f'runs/fix_clock_stdev/{run_name}')
        os.makedirs(run_directory, exist_ok=True)

        # Compile the BEAST-XML from the data and the run-configuration
        ie_xml_str = compile_ielex_xml(ielex_df, **kwargs)

        # Write the XML to ´N_RUNS different files´
        fname_base = run_name.replace("/", "_")
        for i in range(1, N_RUNS+1):
            # Create folder structure
            # if kwargs['use_contactrees']:
            #     if not os.path.exists(run_directory / f'wordtrees_run{i}/'):
            #         os.mkdir(run_directory / f'wordtrees_run{i}/')

            fname = f'{fname_base}_{i}.xml'
            with open(run_directory / fname, 'w') as ie_xml_file:
                ie_xml_file.write(ie_xml_str)

        # Create a shell script to start all runs in parallel
        with open(run_directory / 'start_runs.sh', 'w') as ie_xml_file:
            lines = []
            for i in range(1, N_RUNS + 1):
                threads = xml_snippets.SAMPLER_THREADS[kwargs['sampler']]
                lines.append(f'beast -threads {threads} -overwrite {fname_base}_{i}.xml > {fname_base}_{i}.screenlog &')
            ie_xml_file.write('\n'.join(lines))

    loans = parse_loanwords(ielex_df)
    loans.to_csv(LOANS_CSV_PATH)
