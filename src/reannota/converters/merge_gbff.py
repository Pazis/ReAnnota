"""Merge annotation data into GenBank files."""

import logging

import pandas as pd
from Bio import SeqIO

from reannota.parsers import parse_pseudogff_to_dict

logger = logging.getLogger("ReAnnota")

DESCRIPTION_LENGTH_MAX = 50


def merge_csv_to_gbff(egg_file, interpro_file, gbff_in, gbff_out , pseudofile=None ):
    """
    Merge functional annotations from EggNOG and InterPro CSV files into a GenBank (GBFF) file.

    Parameters:
    -----------
    egg_file : str
        Path to the EggNOG annotation CSV file (tab-delimited) containing functional information.
    interpro_file : str
        Path to the InterPro annotation CSV file (tab-delimited) containing InterPro and GO annotations.
    gbff_in : str
        Path to the input GenBank (.gbff) file to be annotated.
    gbff_out : str
        Path to the output GenBank (.gbff) file with merged annotations.


    1. Read the EggNOG and InterPro CSV files into pandas DataFrames.
    2. Convert the DataFrames into dictionaries keyed by "Query_ID"
    3. Parse the input GenBank file into SeqRecord objects using Bio.SeqIO from biopython.
    4. Loop over each CDS feature in the GenBank sequences.
        a. Add InterPro annotations:
            - Update 'product' with InterPro description if available.
            - Add GO and InterPro terms to 'db_xref'.
        b. Add EggNOG annotations:
            - Updates 'product' or 'note' with EggNOG description.
            - Add COG references, COG categories, GO terms, and PFAM information to 'db_xref' or 'note'.
    5. Write the updated sequences to a new _enchanced GenBank output file.
    """
    logger.debug(f"Merging annotations from EggNOG and InterPro into GBFF: {gbff_in}")

    # Read EggNOG and InterPro CSV files as pandas DataFrames
    if egg_file is not None:
        logger.debug(f"Reading EggNOG annotations from: {egg_file}")
        egg_df = pd.read_csv(egg_file, sep="\t")
        logger.debug(f"Loaded {len(egg_df)} EggNOG annotation entries")
    else:
        egg_df = pd.DataFrame()  # empty DataFrame if missing

    if interpro_file is not None:
        logger.debug(f"Reading InterPro annotations from: {interpro_file}")
        ipr_df = pd.read_csv(interpro_file, sep="\t")
        logger.debug(f"Loaded {len(ipr_df)} InterPro annotation entries")
    else:
        ipr_df = pd.DataFrame()  # empty DataFrame if missing

    if pseudofile is not None :
        logger.debug(f"Reading Pseudogenes from: {pseudofile}")
        pseudo_dict = parse_pseudogff_to_dict(pseudofile)
        logger.debug(f"Loaded {len(pseudo_dict)} Pseudogenes")

    # Convert InterPro DataFrame to a dictionary keyed by Query_ID
    ipr_dictionary = {row["Query_ID"]: row.to_dict() for _, row in ipr_df.iterrows()}

    # Convert EggNOG DataFrame to a dictionary keyed by Query_ID
    egg_dictionary = {row["Query_ID"]: row.to_dict() for _, row in egg_df.iterrows()}

    gbff_read = list(SeqIO.parse(gbff_in, "genbank"))
    logger.debug(f"Loaded {len(gbff_read)} sequences from GBFF file")

    # Counters for statistics
    cds_count = 0
    ipr_annotated = 0
    egg_annotated = 0
    pseudogenes = 0

    # Loop through each sequence in the GenBank file
    for sequence in gbff_read:
        # Loop through each feature in the sequence
        for feature in sequence.features:
            if feature.type == "CDS":  # <-- Filter and focus only on CDS
                cds_count += 1
                locus_tag = feature.qualifiers.get("locus_tag", [""])[0]

                # ---------------Adding Interpro annotation to Gbff-----------------#

                if locus_tag in ipr_dictionary:
                    ipr_annotated += 1
                    values = ipr_dictionary[locus_tag]
                    Description_ipr = values.get("Description", "")
                    GO_terms = str(values.get("GO_terms", ""))
                    Interpro_terms = str(values.get("Interpro_terms", ""))

                    # Update 'product' field if InterPro description is available
                    if (
                        Description_ipr is not None
                        and str(Description_ipr).lower() != "nan"
                        and str(Description_ipr).strip() != ""
                    ):
                        feature.qualifiers["product"] = [str(Description_ipr)]

                    feature.qualifiers.setdefault("db_xref", [])  # <--Ensure db_xref field exists
                    db_xrefs = feature.qualifiers["db_xref"]

                    # Add GO terms if not already present
                    if (
                        GO_terms
                        and GO_terms != "nan"
                        and not any(x.startswith("GO:") for x in db_xrefs)
                    ):
                        feature.qualifiers["db_xref"].append(GO_terms)

                    # Add InterPro terms if not already present
                    if (
                        Interpro_terms
                        and Interpro_terms != "nan"
                        and not any(x.startswith("Interpro:") for x in db_xrefs)
                    ):
                        feature.qualifiers["db_xref"].append(Interpro_terms)

                # --------------Adding eggnog annotation to gbff and checking overlaps------#

                if locus_tag in egg_dictionary:
                    egg_annotated += 1
                    values = egg_dictionary[locus_tag]
                    Description = str(values.get("Description", ""))
                    COG_ref = str(values.get("COG_ref", ""))
                    COG_category = str(values.get("COG_category", ""))
                    GOs = str(values.get("GOs", ""))
                    PFAMs = str(values.get("PFAM", ""))
                    KEGGs = str(values.get("KEGG_ko", ""))

                    # Update 'product' or add to 'note' depending on length and existing product
                    feature.qualifiers.setdefault("note", [])
                    if Description and Description != "nan":
                        if len(Description) > DESCRIPTION_LENGTH_MAX:
                            feature.qualifiers["note"].append(
                                f"Description={values['Description']}"
                            )
                        elif feature.qualifiers["product"] == ["hypothetical protein"]:
                            feature.qualifiers["product"] = [Description]

                    feature.qualifiers.setdefault("db_xref", [])  # <--Ensure db_xref field exists
                    db_xrefs = feature.qualifiers["db_xref"]

                    # Add COG reference if not already present
                    if (
                        COG_ref
                        and COG_ref != "nan"
                        and not any(x.startswith("COG:") for x in db_xrefs)
                    ):
                        feature.qualifiers["db_xref"].append(COG_ref)

                    # Add COG category if not already present
                    if (
                        COG_category
                        and COG_category != "nan"
                        and not any(x.startswith("COG:") for x in db_xrefs)
                    ):
                        feature.qualifiers["db_xref"].append(f"COG:{COG_category}")

                    # Add GO terms if not already present
                    if GOs and GOs != "nan" and not any(x.startswith("GO:") for x in db_xrefs):
                        feature.qualifiers["db_xref"].append(GOs)

                    if (
                        KEGGs
                        and KEGGs != "nan"
                        and not any(x.startswith("KEGG:") for x in db_xrefs)
                    ):
                        feature.qualifiers["db_xref"].append(f"KEGG:{KEGGs}")

                    # Add PFAM information as a note
                    if (
                        PFAMs
                        and PFAMs != "nan"
                        and not any(x.startswith("PFAM:") for x in db_xrefs)
                    ):
                        feature.qualifiers["note"].append(f"Eggnog PFAM comment:{PFAMs}")

                    #--------------Add Pseudogenes-------------------
                if pseudofile != None:
                    if locus_tag in pseudo_dict:
                        pseudogenes += 1
                        values = pseudo_dict[locus_tag]
                        Description_pseudo = values.get("attributes", "").replace("note=", "")

                        feature.qualifiers.setdefault("note", [])

                        if feature.qualifiers["product"] == ["hypothetical protein"]:
                            feature.qualifiers["product"] = ["Pseudogene"]
                            feature.qualifiers["note"].append(Description_pseudo)
                        elif feature.qualifiers["product"] != ["hypothetical protein"]:
                            feature.qualifiers["note"].append(Description_pseudo)


    # Write the updated sequences to the output GenBank file
    logger.debug(f"Writing annotated GBFF to: {gbff_out}")
    with open(gbff_out, "w") as out_handle:
        SeqIO.write(gbff_read, out_handle, "genbank")

    logger.info(
        f"Merged annotations: {cds_count} CDS features processed, "
        f"{ipr_annotated} annotated from InterPro, {egg_annotated} annotated from EggNOG"
    )
    return gbff_out
