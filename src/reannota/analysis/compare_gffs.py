"""Compare two GFF files for annotation differences."""

import logging
import os

import pandas as pd

logger = logging.getLogger("ReAnnota")


def compare_gff(gffstart, gff_enhanced, output_csv):
    """
    Compare two GFF files and count annotation types.

    Parameters:
    -----------
    gffstart : str
        Path to the starting/original GFF file.
    gff_enhanced : str
        Path to the enhanced/updated GFF file.
    output_csv : str
        Path to output CSV file with comparison results.

    Returns:
    --------
    str
        Path to the output CSV file.
    """

    def count_features(gff_path):
        """Count annotation features in a GFF file."""
        counts = {"GO": 0, "InterPro": 0, "PFAM": 0, "KEGG": 0,"Pseudogene_candidates":0,"BGCs":0, "Hypothetical": 0}
        in_antismash=False
        with open(gff_path, "r") as f:
            for line in f:
                if line.startswith("## AntiSMASH predicted BGC features"):
                    in_antismash = True
                    continue
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue

                if not in_antismash:
                    attributes = parts[8]
                    if "GO:" in attributes:
                        counts["GO"] += 1
                    if "InterPro:" in attributes:
                        counts["InterPro"] += 1
                    if "PFAM:" in attributes:
                        counts["PFAM"] += 1
                    if "KEGG" in attributes:
                        counts["KEGG"] += 1
                    if "product=hypothetical" in attributes:
                        counts["Hypothetical"] += 1
                    if "Note=Pseudogene" in attributes:
                        counts["Pseudogene_candidates"] +=1
                else:
                    # AntiSMASH section: count by feature type in column 3
                    type = parts[2]
                    if "biosynthetic-gene-cluster" in type:
                        counts["BGCs"] += 1

        return counts

    # Count features for each file
    start_counts = count_features(gffstart)
    enhanced_counts = count_features(gff_enhanced)

    # Create DataFrame
    data = {
        "File": [os.path.basename(gffstart), os.path.basename(gff_enhanced)],
        "GO_entries": [start_counts["GO"], enhanced_counts["GO"]],
        "InterPro_entries": [start_counts["InterPro"], enhanced_counts["InterPro"]],
        "PFAM_entries": [start_counts["PFAM"], enhanced_counts["PFAM"]],
        "KEGG_entries": [start_counts["KEGG"], enhanced_counts["KEGG"]],
        "Pseudogene_candidates":[start_counts["Pseudogene_candidates"], enhanced_counts ["Pseudogene_candidates"]],
        "BGCs": [start_counts["BGCs"], enhanced_counts["BGCs"]],
        "Hypotheticals": [start_counts["Hypothetical"], enhanced_counts["Hypothetical"]],
    }
    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False, sep="\t")

    # Calculate differences and percentages
    diffs = {}
    for key in ["GO", "InterPro", "PFAM", "KEGG","Pseudogene_candidates","BGCs" ,"Hypothetical"]:
        start_val = start_counts[key]
        enh_val = enhanced_counts[key]
        diff = enh_val - start_val
        if start_val > 0:
            percent = (diff / start_val) * 100
        else:
            percent = float("inf") if enh_val > 0 else 0.0
        diffs[key] = (diff, percent)

    # === Logging summary ===
    logger.info("=== GFF Comparison Summary ===")
    logger.info(f"Input file: {os.path.basename(gffstart)}")
    logger.info(f"Enhanced file: {os.path.basename(gff_enhanced)}")
    for key, (diff, percent) in diffs.items():
        sign = "+" if diff >= 0 else "-"
        if percent == float("inf"):
            percent_str = "∞%"
        else:
            percent_str = f"{percent:+.2f}%"
        logger.info(f"{key:<15}: {sign}{abs(diff)} ({percent_str}) change")

    keys_for_total = ["GO", "InterPro", "PFAM", "KEGG", "Pseudogene_candidates","BGCs"]
    total_features_start = sum(start_counts[k] for k in keys_for_total)
    total_features_enhanced = sum(enhanced_counts[k] for k in keys_for_total)
    total_diff = total_features_enhanced - total_features_start
    total_percent = (total_diff / total_features_start * 100) if total_features_start > 0 else 0

    logger.info(
        f"Total features: {total_features_start} → {total_features_enhanced} ({total_percent:+.2f}%)"
    )
    logger.info(f"Comparison CSV saved: {output_csv}")
    logger.info("======================================")

    return output_csv
