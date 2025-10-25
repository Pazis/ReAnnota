"""Parsers for annotation file formats (EggNOG, InterPro, etc.)."""

from ReAnnota.parsers.antismash import build_gff_rows, load_regions_json
from ReAnnota.parsers.eggnog import (
    build_egg_dictionary_clean,
    egg_dict_to_tsv,
    egg_gff_to_dataframe,
)
from ReAnnota.parsers.gecco import (
    build_gecco_gff_rows,
    combine_gbk_files,
    extract_cluster_genes,
    read_file_paths,
)
from ReAnnota.parsers.interpro import ipr_dictotsv, ipr_termfinder
from ReAnnota.parsers.pseudofinder import parse_pseudogff_to_dict

__all__ = [
    "build_egg_dictionary_clean",
    "egg_dict_to_tsv",
    "egg_gff_to_dataframe",
    "ipr_termfinder",
    "ipr_dictotsv",
    "load_regions_json",
    "build_gff_rows",
    "parse_pseudogff_to_dict",
    "build_gecco_gff_rows",
    "extract_cluster_genes",
    "read_file_paths",
    "combine_gbk_files"

]
