"""InterPro annotation file parser."""

import csv
import logging

logger = logging.getLogger("cyano_annotation")


def ipr_termfinder(input_file):
    """
    Parses the intepro GFF3-like file to extract GO terms, InterPro terms, and functional descriptions.

    Parameters:
    -----------
    input_file : str
        Path to the input file, typically a GFF3 file or similar tab-delimited format.

    Returns:
    --------
    dict
        A dictionary keyed by Query_ID, with values containing:
        - "GO": list of GO terms
        - "InterPro": list of InterPro terms
        - "Description": list of functional descriptions
    """
    logger.debug(f"Reading InterPro file: {input_file}")
    dictionary = {}
    line_count = 0

    with open(input_file) as in_handle:
        for line in in_handle:
            line_count += 1
            if line.startswith("#"):  # <-- Skip comment lines
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:  # <--Skip incomplete lines
                continue

            query_id = parts[0]  # Feature_ID
            attributes = parts[8]  # Column 9 contains all the attributes

            # Terms that should be ignored in the description
            outsiders = [
                "Protein of unknown function",
                "Domain of unknown function",
                "Region",
                "Signal",
            ]

            # Extract GO terms
            go_terms = []
            if '"GO:' in attributes:
                for terms in attributes.split(";"):
                    if '"GO:' in terms:
                        for values in terms.replace("=", ",").split(","):
                            values = values.strip()
                            if values.startswith('"GO:') and values not in go_terms:
                                go_terms.append(values)

            # Extract InterPro terms
            ipr_terms = []
            if '"InterPro:' in attributes:
                for terms in attributes.split(";"):
                    if '"InterPro:' in terms:
                        for values in terms.replace("=", ",").split(","):
                            values = values.strip()
                            if values.startswith('"InterPro:') and values not in ipr_terms:
                                ipr_terms.append(values)

            # Extract functional descriptions (ignoring generic/unknown terms)
            description = []
            for terms in attributes.split(";"):
                if terms.startswith("signature_desc="):
                    value = terms.split("=", 1)[1].strip()
                    if not any(value.startswith(x) for x in outsiders):
                        description.append(value)

            if go_terms or ipr_terms:
                dictionary[query_id] = {
                    "GO": go_terms,
                    "InterPro": ipr_terms,
                    "Description": description,
                }

    logger.info(
        f"InterPro processing: {len(dictionary)} entries with annotations found "
        f"from {line_count} lines"
    )
    return dictionary


def ipr_dictotsv(dictionary, output_file):
    """
    Write the annotation dictionary into a tab-delimited TSV file.

    dictionary : dict
    Dictionary returned by termfinder(), containing GO, InterPro, and Description for each Query_ID.
    output_file : str
    Path to the output TSV file.
    """
    logger.debug(f"Writing InterPro results to TSV: {output_file}")
    with open(output_file, "w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")

        # Write header
        writer.writerow(["Query_ID", "GO_terms", "Interpro_terms", "Description"])

        for query_id, vals in dictionary.items():
            go_terms = vals["GO"] or [""]  # <--Use empty string if no GO terms
            ipr_terms = vals["InterPro"] or [""]  # <-- Use empty string if no InterPro terms
            desc = ",".join(vals["Description"])  # <--Combine descriptions into a single string

            max_len = max(len(go_terms), len(ipr_terms))

            # Write each GO/IPR term on a separate line, reusing description
            for i in range(max_len):
                go_val = go_terms[i] if i < len(go_terms) else ""
                ipr_val = ipr_terms[i] if i < len(ipr_terms) else ""
                writer.writerow([query_id, go_val, ipr_val, desc])

    return output_file
