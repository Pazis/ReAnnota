"""Convert GenBank (GBFF) files to GFF3 format."""

import logging
from collections import defaultdict

from Bio import SeqIO

from reannota.parsers import build_gff_rows, load_antismash_gbk, load_gecco

logger = logging.getLogger("ReAnnota")


def gbff_to_gff(gbff_in, gff_out, antismash_gbk=None, gecco_csv=None, gecco_summary=None, bakta_gff=None):
    """
    Convert GenBank (GBFF) to GFF3 and optionally merge antiSMASH predictions 
    and original CRISPR coordinates from Bakta GFF.

    Parameters:
    - gbff_in: input GenBank file
    - gff_out: output GFF3 file
    - antismash_gbk: optional antiSMASH genbank file
    - gecco_csv: optional GECCO csv file
    - gecco_summary: optional GECCO summary tsv file
    - bakta_gff: optional original Bakta GFF file to rescue missing CRISPR coordinates
    """
    logger.debug(f"Converting GBFF to GFF3: {gbff_in} → {gff_out}")

    feature_counts = {"CDS": 0, "tRNA": 0, "rRNA": 0, "tmRNA": 0, "ncRNA": 0, "repeat_region": 0, "CRISPR_arrays": 0}

    # Load antiSMASH BGCs if provided
    all_bgcs = []

    # Load AntiSMASH (expects .gbk now)
    if antismash_gbk and str(antismash_gbk) != "None":
        logging.info("antiSMASH file detected | Integrating antiSMASH output")
        try:
            as_bgcs = load_antismash_gbk(str(antismash_gbk))
            all_bgcs.extend(as_bgcs)
        except Exception as e:
            logger.error(f"Failed to load antiSMASH file: {e}")

    # Load GECCO (expects .tsv summary + .csv paths)
    if gecco_summary and str(gecco_summary) != "None":
        logging.info("GECCO files detected | Integrating GECCO output")
        try:
            ge_bgcs = load_gecco(str(gecco_summary), gbk_paths_csv=str(gecco_csv))
            all_bgcs.extend(ge_bgcs)
        except Exception as e:
            logger.error(f"Failed to load GECCO files: {e}")

    # Load CRISPR features from Bakta GFF
    crispr_by_seqid = defaultdict(list)
    if bakta_gff and str(bakta_gff) != "None":
        logging.info("Bakta GFF detected | Extracting CRISPR coordinates")
        try:
            with open(bakta_gff, "r") as bgff:
                for line in bgff:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    
                    if len(parts) >= 9:
                        seqid = parts[0]
                        feature_type = parts[2]
                        
                        # Catch 'CRISPR', 'crispr-repeat', and 'crispr-spacer'
                        if "CRISPR" in feature_type.upper():
                            start_coord = int(parts[3])
                            crispr_by_seqid[seqid].append((start_coord, line.strip()))
                            
                            # Only count the main array for the summary log
                            if feature_type.upper() == "CRISPR":
                                feature_counts["CRISPR_arrays"] += 1
                                
        except Exception as e:
            logger.error(f"Failed to parse Bakta GFF for CRISPRs: {e}")

    with open(gff_out, "w") as out:
        out.write("##gff-version 3\n")  # GFF3 header

        # Parse GenBank records one by one
        for record in SeqIO.parse(gbff_in, "genbank"):
            seqid = record.id
            contig_length = len(record.seq)

            # Write sequence-region line for Circos/Browsers
            out.write(f"##sequence-region {seqid} 1 {contig_length}\n")

            # Ensure an organism annotation exists
            if "organism" not in record.annotations:
                record.annotations["organism"] = "unknown_organism"

            contig_features = [] # Store features here to sort them later

            # Iterate through features of interest
            for feature in record.features:
                if feature.type in ["CDS", "tRNA", "rRNA", "tmRNA", "ncRNA", "repeat_region"]:
                    feature_counts[feature.type] += 1

                    # Convert feature location to 1-based coordinates (GFF standard)
                    start = int(feature.location.start) + 1
                    end = int(feature.location.end)
                    strand = "+" if feature.location.strand == 1 else "-"

                    # Collect feature attributes
                    attributes = []
                    if "locus_tag" in feature.qualifiers:
                        attributes.append(f"ID={feature.qualifiers['locus_tag'][0]}")
                    if "gene" in feature.qualifiers:
                        attributes.append(f"Name={feature.qualifiers['gene'][0]}")
                    if "product" in feature.qualifiers:
                        attributes.append(f"product={feature.qualifiers['product'][0]}")
                    if "note" in feature.qualifiers:
                        attributes.append(f"Note={','.join(feature.qualifiers['note'])}")
                    if "db_xref" in feature.qualifiers:
                        attributes.append(f"Dbxref={','.join(feature.qualifiers['db_xref'])}")
                    if "pseudogene" in feature.qualifiers:
                        attributes.append(f"pseudogene={','.join(feature.qualifiers['pseudogene'])}")

                    # Join attributes and format GFF3 line
                    attr_str = ";".join(attributes)
                    gff_line = f"{seqid}\tBiopython\t{feature.type}\t{start}\t{end}\t.\t{strand}\t.\t{attr_str}"
                    
                    # Note: Using a tuple with (start_coord, feature_type, line) 
                    # ensures parent 'CRISPR' features are written right before their repeats/spacers if they share a start coord
                    contig_features.append((start, feature.type, gff_line))

            # Add extracted CRISPR coordinates for this contig
            if seqid in crispr_by_seqid:
                for start_coord, gff_line in crispr_by_seqid[seqid]:
                    # Extract feature type from the saved line just for sorting purposes
                    ftype = gff_line.split("\t")[2]
                    contig_features.append((start_coord, ftype, gff_line))

            # Sort all features (GBFF + CRISPRs)
            # Sorting by [0] (start coordinate). 
            # If start coordinates tie, it sorts by [1] (feature_type string) which keeps things deterministic.
            contig_features.sort(key=lambda x: (x[0], x[1]))

            # Write the sorted features to the GFF out file
            for _, _, line in contig_features:
                out.write(line + "\n")

        # Write antiSMASH/GECCO features after the genomic features
        if all_bgcs:
            out.write("## Predicted BGC features (antiSMASH / GECCO) \n")
            for row in build_gff_rows(all_bgcs):
                out.write("\t".join(row) + "\n")

    logger.info(
        f"GFF3 conversion complete: {feature_counts['CDS']} CDS, "
        f"{feature_counts['tRNA']} tRNA, {feature_counts['rRNA']} rRNA, {feature_counts['ncRNA']} ncRNA, "
        f"{feature_counts['CRISPR_arrays']} CRISPR arrays integrated and written"
    )
    return gff_out