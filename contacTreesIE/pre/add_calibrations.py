import argparse
from pathlib import Path
import xml.etree.ElementTree as ET

from tqdm import tqdm
import pycldf
import pyglottolog
from cldfcatalog import Catalog


def normal(mean, std):
    return dict(
        tag="Normal", name="distr", offset="0.0", mean=f"{mean:}", sigma=f"{std:}"
    )


def until(lower, upper):
    # Assume lower to upper is the 2 σ interval of a normal distribution
    mean, std = (lower + upper) / 2, (upper - lower) / 4
    return normal(mean, std)


IE_CALIBRATIONS = [
    {"languages": {"Hittite"}, "d": normal(3400, 100)},
    {"languages": {"Vedic_Sanskrit"}, "d": normal(3250, 250)},
    {"languages": {"Avestan"}, "d": normal(2500, 50)},
    {"languages": {"Ancient_Greek"}, "d": normal(2450, 50)},
    {"languages": {"Latin"}, "d": normal(2150, 50)},
    {"languages": {"Gothic"}, "d": normal(1650, 25)},
    {"languages": {"Old_High_German"}, "d": normal(1150, 50)},
    {"languages": {"Old_English"}, "d": normal(1000, 50)},
    {"languages": {"Old_Norse"}, "d": normal(800, 50)},
    {"languages": {"Classical_Armenian"}, "d": normal(1550, 50)},
    {"languages": {"Tocharian_B"}, "d": normal(1350, 150)},
    {"languages": {"Old_Irish"}, "d": normal(1200, 100)},
    {"languages": {"Cornish"}, "d": normal(300, 100)},
    {"languages": {"Old_Church_Slavonic"}, "d": normal(1000, 50)},
    {
        "languages": {
            "Gothic",
            "Old_Norse",
            "Icelandic_ST",
            "Faroese",
            "Norwegian",
            "Swedish",
            "Danish",
            "Old_English",
            "English",
            "Frisian",
            "Old_High_German",
            "German",
            "Luxembourgish",
            "Schwyzerdutsch",
            "Dutch_List",
            "Flemish",
            "Afrikaans",
        },
        "name": "Germanic",
        "d": {"tag": "Uniform", "name": "distr", "lower": "2250", "upper": "20000"},
    },
    {
        "languages": {
            "Latin",
            "Sardinian_N",
            "Sardinian_C",
            "Rumanian_List",
            "Catalan",
            "Portuguese_ST",
            "Spanish",
            "French",
            "Provencal",
            "Walloon",
            "Ladin",
            "Romansh",
            "Friulian",
            "Italian",
        },
        "name": "Romance",
        "d": {"tag": "Uniform", "name": "distr", "lower": "1750", "upper": "20000"},
    },
    {
        "languages": {
            "Old_Norse",
            "Icelandic_ST",
            "Faroese",
            "Norwegian",
            "Swedish",
            "Danish",
        },
        "name": "Scandinavian",
        "d": {"tag": "Uniform", "name": "distr", "lower": "1500", "upper": "20000"},
    },
    {
        "languages": {
            "Czech",
            "Slovak",
            "Polish",
            "Upper_Sorbian",
            "Ukrainian",
            "Byelorussian",
            "Russian",
            "Slovenian",
            "Macedonian",
            "Bulgarian",
            "Serbian",
            "Old_Church_Slavonic",
        },
        "name": "Slavic",
        "d": {"tag": "Uniform", "name": "distr", "lower": "1500", "upper": "20000"},
    },
    {
        "languages": {
            "Lithuanian_ST",
            "Latvian",
        },
        "name": "East_Baltic",
        "d": {"tag": "Uniform", "name": "distr", "lower": "1300", "upper": "20000"},
    },
    {
        "languages": {
            "Welsh_N",
            "Breton_ST",
            "Cornish",
        },
        "name": "British_Celtic",
        "d": {"tag": "Uniform", "name": "distr", "lower": "1250", "upper": "20000"},
    },
    {
        "languages": {
            "Irish_B",
            "Scots_Gaelic",
        },
        "name": "Modern_Irish-Scots_Gaelic",
        "d": {"tag": "Uniform", "name": "distr", "lower": "1050", "upper": "20000"},
    },
    {
        "languages": {
            "Tadzik",
            "Persian",
        },
        "name": "Persian-Tajik",
        "d": {"tag": "Uniform", "name": "distr", "lower": "750", "upper": "20000"},
    },
]


FBD_REPLACEMENTS = {
        "TipDatesRandomWalker": "beast.evolution.operators.SampledNodeDateRandomWalker",
}


def calibration(run, prior, trait, all_languages, d, languages=None, glottolog_clade=None,
                mean=0.0, name=None, replacements=None, monophyletic=False):
    if languages is None:
        languages = []
    if replacements is None:
        replacements = {}

    if glottolog_clade is not None:
        languages = {
            l for l, lineage in all_languages.items() if glottolog_clade in lineage
        }
    if name is None:
        if glottolog_clade is None:
            name = '_'.join(languages)
        else:
            name = glottolog_clade

    if mean == 0.0:
        mean = d.get("mean", mean)

    tag = d.pop('tag')

    if len(languages) == 0:
        return
    elif len(languages) == 1:
        language = list(languages)[0]
        if glottolog_clade:
            mrcaprior = ET.SubElement(
                prior,
                "distribution",
                id=f"{language:}_originateMRCA",
                monophyletic="true" if monophyletic else "false",
                spec="beast.math.distributions.MRCAPrior",
                tree="@tree",
                useOriginate="true",
            )
            taxonset = ET.SubElement(
                mrcaprior, "taxonset", id=f"tx_{language:}", spec="TaxonSet"
            )
            ET.SubElement(taxonset, "taxon", idref=f"{language:}")
            ET.SubElement(mrcaprior, tag, **d)
        else:
            if not trait.text or not trait.text.strip():
                trait.text = f"\n{language:} = {mean:}"
            else:
                trait.text += f",\n{language:} = {mean:}"

            op = ET.SubElement(
                run,
                "operator",
                id=f"TipDatesandomWalker:{language:}",
                spec=replacements.get("TipDatesRandomWalker", "beast.evolution.operators.TipDatesRandomWalker"),
                windowSize="1",
                tree="@tree",
                weight="3.0",
            )
            ET.SubElement(op, "taxonset", idref=f"{language:}_tip")

            mrcaprior = ET.SubElement(
                prior,
                "distribution",
                id=f"{language:}_tipMRCA",
                monophyletic="true" if monophyletic else "false",
                spec="beast.math.distributions.MRCAPrior",
                tree="@tree",
                tipsonly="true",
            )
            taxonset = ET.SubElement(
                mrcaprior, "taxonset", id=f"{language:}_tip", spec="TaxonSet"
            )
            ET.SubElement(taxonset, "taxon", idref=f"{language:}")
            ET.SubElement(mrcaprior, tag, **d)
    else:
        mrcaprior = ET.SubElement(
            prior,
            "distribution",
            id=f"{name}_tipMRCA",
            monophyletic="true" if monophyletic else "false",
            spec="beast.math.distributions.MRCAPrior",
            tree="@tree",
        )
        taxonset = ET.SubElement(
            mrcaprior, "taxonset", id=f"{name}", spec="TaxonSet"
        )
        plate = ET.SubElement(
            taxonset, "plate", range=",".join(sorted(languages)), var="language"
        )
        ET.SubElement(plate, "taxon", idref="$(language)")
        ET.SubElement(mrcaprior, tag, **d)
    mrcaprior.tail = "\n"


def add_ie_calibrations(calibrations, first_writing):
    parser = argparse.ArgumentParser(
        description="Export a CLDF dataset (or similar) to bioinformatics alignments"
    )
    parser.add_argument(
        "--family",
        "-f",
        help="""Only include languages within this Glottolog clade""",
    )
    parser.add_argument(
        "--first-writing",
        "-w",
        type=float,
        help="The date (BP) when writing started in the region, and thus a notable chance of lanugages of languages being sampled begins. (Default: Don't modify this parameter in the template.)"
    )
    parser.add_argument(
        "--subset",
        type=argparse.FileType('r'),
        help="A file (or '-' for stdin) containing one language to be included per line"
    )
    parser.add_argument(
        "--output-file",
        "-o",
        type=Path,
        help="""File to write output to. (If output file exists, add tags in there.) Default: Write to stdout""",
    )
    parser.add_argument(
        "--sampled-ancestors",
        "-s",
        action="store_true",
        default=False,
        help="Work for a sampled ancestor tree, which needs variant operators."
    )
    parser.add_argument(
        "--metadata",
        "-m",
        type=Path,
        default="raw_cldf/cldf-metadata.json",
        help="""Metadata file, for language list""",
    )
    args = parser.parse_args()

    ds = pycldf.Wordlist.from_metadata(args.metadata)

    if args.subset:
        langs = {l.strip() for l in args.subset}

    with Catalog.from_config("glottolog", tag="v4.3") as glottolog_repo:
        glottolog = pyglottolog.Glottolog(glottolog_repo.dir)
        languages = {}
        for language in tqdm(
            ds["LanguageTable"], total=ds["LanguageTable"].common_props["dc:extent"]
        ):
            id = language[ds.column_names.languages.id]
            if args.subset and id not in langs:
                continue
            try:
                languoid = glottolog.languoid(language["Glottocode"])
            except TypeError:
                continue
            ancestors = [languoid.id] + [a.id for a in languoid.ancestors]
            if args.family and args.family not in ancestors:
                continue
            languages[id] = ancestors

    if args.output_file is None:
        root = ET.fromstring(
            """
        <beast><tree /><run><distribution id="posterior" spec="util.CompoundDistribution"><distribution id="prior" spec="util.CompoundDistribution" /></distribution></run></beast>
        """
        )
        et = ET.ElementTree(root)
    elif args.output_file.exists():
        et = ET.parse(args.output_file)
        root = et.getroot()
    else:
        root = ET.fromstring(
            """
        <beast><tree /><run><distribution id="posterior" spec="util.CompoundDistribution"><distribution id="prior" spec="util.CompoundDistribution" /></distribution></run></beast>
        """
        )
        et = ET.ElementTree(root)

    prior = list(root.iter("distribution"))[1]
    assert prior.attrib["id"] == "prior"

    run = list(root.iter("run"))[0]

    traits = list(root.iter("trait"))
    if not traits:
        tree = list(root.iter("tree"))[0]
        trait = ET.SubElement(
            tree,
            "trait",
            id="datetrait",
            spec="beast.evolution.tree.TraitSet",
            taxa="@taxa",
            traitname="date-backward",
        )
    else:
        trait = traits[0]
        assert trait.attrib["traitname"] == "date-backward"

    for lang in languages:
        print(lang)

    for c in calibrations:
        calibration(run, prior, trait, languages, replacements=FBD_REPLACEMENTS if args.sampled_ancestors else {}, **c)

    if args.first_writing is not None:
        run = [tag for tag in root.iter() if tag.attrib.get("id") == "SamplingChangeTime"][0]
        run.text = f"0. {args.first_writing:f}"

    et.write(args.output_file, encoding="unicode")
    return et, root, run, prior, trait, languages
