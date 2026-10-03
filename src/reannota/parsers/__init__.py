"""Parsers for annotation file formats (EggNOG, InterPro, etc.)."""

from reannota.parsers.bgc import build_gff_rows, load_antismash_gbk, load_gecco
from reannota.parsers.eggnog import (
    build_egg_dictionary_clean,
    egg_dict_to_tsv,
    egg_gff_to_dataframe,
)
from reannota.parsers.interpro import ipr_dictotsv, ipr_termfinder
from reannota.parsers.pseudofinder import parse_pseudogff_to_dict

__all__ = [
    "build_egg_dictionary_clean",
    "egg_dict_to_tsv",
    "egg_gff_to_dataframe",
    "ipr_termfinder",
    "ipr_dictotsv",
    "build_gff_rows",
    "parse_pseudogff_to_dict",
    "load_antismash_gbk",
    "load_gecco",

]
