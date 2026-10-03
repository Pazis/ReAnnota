"""InterPro annotation file parser."""

import csv
import logging
import re
from collections import defaultdict

logger = logging.getLogger("ReAnnota")

def ipr_termfinder(input_file):
    """
    Parses the InterPro GFF3-like file to extract GO terms, InterPro terms, 
    and functional descriptions using a hierarchical database algorithm.

    Parameters:
    -----------
    input_file : str
        Path to the input file (GFF3 format).

    Returns:
    --------
    dict
        A dictionary keyed by Query_ID containing accumulated, non-redundant annotations.
    """
    logger.debug(f"Reading InterPro file: {input_file}")
    
    # Using defaultdict prevents the "Overwrite Bug"
    dictionary = defaultdict(lambda: {
        "GO": set(),
        "InterPro": set(),
        "Desc_Tier1": set(), # Whole-protein databases
        "Desc_Tier2": set(), # High-quality domain databases
        "Desc_Tier3": set()  # Structural/Generic databases
    })
    
    line_count = 0

    # Terms that carry no biological meaning
    outsiders = [
        "Protein of unknown function",
        "Domain of unknown function",
        "Region",
        "Signal",
        "Uncharacterized protein",
        "hypothetical protein"
    ]

    # Define database hierarchy for descriptions
    tier1_dbs = {"TIGRFAM", "NCBIFAM", "PANTHER", "HAMAP", "PRINTS"}
    tier2_dbs = {"PFAM", "CDD", "SMART", "PROSITEPROFILES", "PROSITEPATTERNS"}

    with open(input_file) as in_handle:
        for line in in_handle:
            line_count += 1
            if line.startswith("#"):
                continue
            
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            query_id = parts[0]
            analysis_db = parts[1].upper() # e.g., 'PFAM', 'NCBIFAM'
            attributes = parts[8]

            # 1. Extract GO terms using Regex (much safer than string splitting)
            # Looks for exactly "GO:" followed by 7 digits
            go_matches = re.findall(r'GO:\d{7}', attributes)
            dictionary[query_id]["GO"].update(go_matches)

            # 2. Extract InterPro terms using Regex
            # Looks for exactly "IPR" followed by 6 digits
            ipr_matches = re.findall(r'IPR\d{6}', attributes)
            dictionary[query_id]["InterPro"].update([f"InterPro:{ipr}" for ipr in ipr_matches])

            # 3. Extract and Categorize Descriptions
            desc_match = re.search(r'signature_desc=([^;]+)', attributes)
            if desc_match:
                desc = desc_match.group(1).strip()
                
                # Filter out useless descriptions
                if not any(desc.startswith(x) for x in outsiders):
                    
                    # Route the description to the correct biological Tier
                    if analysis_db in tier1_dbs:
                        dictionary[query_id]["Desc_Tier1"].add(desc)
                    elif analysis_db in tier2_dbs:
                        dictionary[query_id]["Desc_Tier2"].add(desc)
                    else:
                        dictionary[query_id]["Desc_Tier3"].add(desc)

    logger.info(f"InterPro processing: {len(dictionary)} unique proteins found from {line_count} lines")
    return dictionary


def ipr_dictotsv(dictionary, output_file):
    """
    Write the hierarchical annotation dictionary into a tab-delimited TSV file,
    resolving the best description for each protein.
    """
    logger.debug(f"Writing InterPro results to TSV: {output_file}")
    
    with open(output_file, "w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")

        # Note: Gene_name is intentionally left blank. InterPro 'Name=' is a signature ID (e.g. PF1234), not a gene name.
        writer.writerow(["Query_ID", "GO_terms", "Interpro_terms", "Description", "Gene_name"])

        for query_id, vals in dictionary.items():
            
            # Join non-redundant GO and IPR terms
            go_terms_string = ",".join(sorted(vals["GO"]))
            ipr_terms_string = ",".join(sorted(vals["InterPro"]))
            
            # --- THE BEST MATCH ALGORITHM ---
            # Attempt 1: Use Tier 1 (Whole Protein names like "ATP-dependent DNA helicase RecG")
            if vals["Desc_Tier1"]:
                # If there are multiple, pick the longest one (usually the most descriptive)
                best_desc = max(vals["Desc_Tier1"], key=len)
            
            # Attempt 2: Use Tier 2 (Specific Domains)
            elif vals["Desc_Tier2"]:
                # A protein might have multiple valid domains. Join them with a semicolon.
                # Example: "DEAD/DEAH box helicase; RecG wedge domain"
                best_desc = "; ".join(sorted(vals["Desc_Tier2"]))
            
            # Attempt 3: Fallback to Tier 3
            elif vals["Desc_Tier3"]:
                best_desc = "; ".join(sorted(vals["Desc_Tier3"]))
            
            # No description found
            else:
                best_desc = ""

            # Write row (Leaving Gene_name empty so it doesn't pollute the GenBank file)
            writer.writerow([query_id, go_terms_string, ipr_terms_string, best_desc, ""])

    return output_file