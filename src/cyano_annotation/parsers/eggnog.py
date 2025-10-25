"""EggNOG annotation file parser."""

import csv
import logging
import re

import pandas as pd

logger = logging.getLogger("cyano_annotation")

E_VALUE_THRESHOLD = 1e-5


def egg_gff_to_dataframe(input_file):
    """
    Read a tab-delimited EggNOG/GFF file and return a cleaned DataFrame with selected columns.

    Parameters:
    input_file : str
    Path to the input tab-delimited file.

    Returns:
    pd.DataFrame
    DataFrame containing only the relevant columns for downstream processing.
    """

    # Define all column names expected in the EggNOG/GFF file
    col_names = [
        "query",
        "seed_ortholog",
        "evalue",
        "score",
        "eggNOG_OGs",
        "max_annot_lvl",
        "COG_category",
        "Description",
        "Preferred_name",
        "GOs",
        "EC",
        "KEGG_ko",
        "KEGG_Pathway",
        "KEGG_Module",
        "KEGG_Reaction",
        "KEGG_rclass",
        "BRITE",
        "KEGG_TC",
        "CAZy",
        "BiGG_Reaction",
        "PFAMs",
    ]

    # Read the file into a DataFrame
    df = pd.read_csv(input_file, sep="\t", header=None, names=col_names, comment="#", skiprows=4)

    # Keep only relevant columns for annotation
    keep_cols = [
        "query",
        "evalue",
        "eggNOG_OGs",
        "COG_category",
        "Description",
        "GOs",
        "PFAMs",
        "KEGG_ko",
    ]

    clean_df = df[keep_cols]
    return clean_df


def build_egg_dictionary_clean(input_file):
    """
    Build a dictionary of cleaned annotations from the EggNOG/GFF file, filter by e-value and ignore generic terms.

    Parameters:
    input_file : str
        Path to the input EggNOG/GFF file.

    Returns:
    dict
        Dictionary keyed by Query_ID with values containing lists of:
        - GOs
        - PFAM
        - COG_category
        - COG_ref
        - Description
    """
    logger.debug(f"Reading EggNOG file: {input_file}")
    clean_df = egg_gff_to_dataframe(input_file)
    logger.debug(f"Loaded {len(clean_df)} rows from EggNOG file")

    dictionary = {row["query"]: row.drop("query").to_dict() for _, row in clean_df.iterrows()}
    dictionary_clean = {}

    for query, values in dictionary.items():
        query_id = query
        description = str(values["Description"])
        evalue = float(values["evalue"])
        eggNOG_OGs = str(values["eggNOG_OGs"])
        COG_category = str(values["COG_category"])
        GOs = str(values["GOs"])
        PFAMs = str(values["PFAMs"])
        KEGGs = str(values["KEGG_ko"])
        outsiders = [
            "DUF",
            "Domain",
            "Protein of",
            "Protein conscerved",
            "multi",
            "Alternative",
            "restriction",
            "alpha",
            "domain",
        ]

        if evalue < E_VALUE_THRESHOLD:  # <--Filter by e-value

            # Collect GO terms if valid
            go_terms = []
            if "-" not in GOs:
                go_terms.append(GOs)

            # Collect PFAM terms, ignoring generic descriptions
            pfam_terms = []
            if "-" not in PFAMs:
                if not any(PFAMs.startswith(x) for x in outsiders):
                    pfam_terms.append(PFAMs)

            # Collect COG categories
            COG_terms = []
            if "-" not in COG_category:
                COG_terms.append(COG_category)

            # Collect description if not generic
            Desc = []
            if "-" not in description:
                if not any(description.startswith(x) for x in outsiders):
                    Desc.append(description)

            # Extract COG reference from eggNOG_OGs
            EGG_COG = []
            if "COG" in str(eggNOG_OGs):
                match = re.search(r"(COG[^@]+)@", str(eggNOG_OGs))
                if match:
                    EGG_COG.append(match.group(1))

            KEGG = []
            if KEGGs != "-" and KEGGs.strip():
                # Clean each KO
                cleaned = [k.strip().replace("ko:K", "KO").upper() for k in KEGGs.split(",")]
                # Append each cleaned KO to the KEGG list
                for k in cleaned:
                    KEGG.append(k)

        # Only keep if there are clean hits
        if go_terms or pfam_terms or Desc or KEGG or EGG_COG:
            dictionary_clean[query_id] = {
                "GOs": go_terms,
                "PFAM": pfam_terms,
                "COG_category": COG_terms,
                "COG_ref": [f"COG:{cog}" for cog in EGG_COG],
                "KEGGs": KEGG,
                "Description": Desc,
            }

    logger.info(
        f"EggNOG processing: {len(dictionary_clean)}/{len(dictionary)} entries passed "
        f"filtering (e-value < {E_VALUE_THRESHOLD})"
    )
    return dictionary_clean


def egg_dict_to_tsv(dictionary, output_file):
    """
    Write the cleaned annotation dictionary to a tab-delimited TSV file.

    Parameters:
    dictionary : dict
        Annotation dictionary from build_dictionary_clean().
    output_file : str
        Path to the output TSV file.
    """
    logger.debug(f"Writing EggNOG results to TSV: {output_file}")
    with open(output_file, "w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")

        writer.writerow(
            ["Query_ID", "GOs", "PFAM", "COG_category", "COG_ref", "KEGG_ko", "Description"]
        )

        for query_id, vals in dictionary.items():
            go_terms = vals.get("GOs", []) or [""]
            pfam_terms = vals.get("PFAM", []) or [""]
            desc = ",".join(vals.get("Description", []))
            cog_terms = vals.get("COG_category", []) or [""]
            cog_ref_terms = vals.get("COG_ref", []) or [""]
            kegg_terms = vals.get("KEGGs", []) or [""]
            max_len = max(len(go_terms), len(pfam_terms), len(cog_terms), len(cog_ref_terms))

            for i in range(max_len):
                go_val = go_terms[i] if i < len(go_terms) else ""
                pfam_val = pfam_terms[i] if i < len(pfam_terms) else ""
                cog_val = cog_terms[i] if i < len(cog_terms) else ""
                cog_ref_val = cog_ref_terms[i] if i < len(cog_ref_terms) else ""
                kegg_val = kegg_terms[i] if i < len(kegg_terms) else ""
                writer.writerow([query_id, go_val, pfam_val, cog_val, cog_ref_val, kegg_val, desc])

    return output_file
