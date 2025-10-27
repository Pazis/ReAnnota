"""Convert GenBank (GBFF) files to GFF3 format."""

import logging

from Bio import SeqIO

from reannota.parsers import (
    build_gecco_gff_rows,
    build_gff_rows,
    extract_cluster_genes,
    load_regions_json,
)

logger = logging.getLogger("ReAnnota")


def gbff_to_gff(gbff_in, gff_out, antismash_json=None,gecco_gbk=None, antismash_version="7.1.0"):
    """
    Convert GenBank (GBFF) to GFF3 and optionally merge antiSMASH predictions.

    Parameters:
    - gbff_in: input GenBank file
    - gff_out: output GFF3 file
    - antismash_json: optional antiSMASH regions.json file
    - antismash_version: version string for antiSMASH features
    """
    logger.debug(f"Converting GBFF to GFF3: {gbff_in} → {gff_out}")

    feature_counts = {"CDS": 0, "tRNA": 0, "rRNA": 0 , "tmRNA": 0}

    # Load antiSMASH BGCs if provided
    antismash_rows = []
    if antismash_json != None:
        if antismash_json != "None":
            logging.info("antiSMASH file detected | Intergrating antiSMASH output")
            bgcs = load_regions_json(antismash_json)
            antismash_rows = list(build_gff_rows(bgcs, antismash_version))

    # Load Gecco BGCs if provided
    gecco_rows=[]
    if gecco_gbk != None:
        if gecco_gbk != "None":
            logging.info("Gecco file detected | Intergrating Gecco output")
            gecco_bgcs = extract_cluster_genes(gecco_gbk)
            gecco_rows = list(build_gecco_gff_rows(gecco_bgcs))



    with open(gff_out, "w") as out:
        out.write("##gff-version 3\n")  # GFF3 header

        # Parse GenBank records one by one
        for record in SeqIO.parse(gbff_in, "genbank"):
            seqid = record.id
            contig_length = len(record.seq)

            # Write sequence-region line for Circos
            out.write(f"##sequence-region {seqid} 1 {contig_length}\n")

            # Ensure an organism annotation exists
            if "organism" not in record.annotations:
                record.annotations["organism"] = "unknown_organism"

            # Iterate through features of interest
            for feature in record.features:
                if feature.type in ["CDS", "tRNA", "rRNA", "tmRNA"]:
                    feature_counts[feature.type] += 1

                    # Convert feature location to 1-based coordinates (GFF standard)
                    start = int(feature.location.start) + 1
                    end = int(feature.location.end)
                    strand = "+" if feature.location.strand == 1 else "-"

                    # Collect feature attributes
                    attributes = []
                    if "locus_tag" in feature.qualifiers:
                        attributes.append(f"ID={feature.qualifiers['locus_tag'][0]}")
                    if "product" in feature.qualifiers:
                        attributes.append(f"product={feature.qualifiers['product'][0]}")
                    if "note" and "pseudogene" in feature.qualifiers:
                        attributes.append(f"Note=Pseudogene candidate. Reason {','.join(feature.qualifiers['note'])}")
                    elif "note" in feature.qualifiers:
                        attributes.append(f"Note={','.join(feature.qualifiers['note'])}")
                    else:
                        attributes.append("Note=")
                    if "db_xref" in feature.qualifiers:
                        attributes.append(f"Dbxref={','.join(feature.qualifiers['db_xref'])}")

                    # Join attributes and format GFF3 line
                    attr_str = ";".join(attributes)
                    gff_line = f"{seqid}\tBiopython\t{feature.type}\t{start}\t{end}\t.\t{strand}\t.\t{attr_str}\n"
                    out.write(gff_line)

        # Write antiSMASH features after GenBank features
        if antismash_rows:
            out.write("## AntiSMASH predicted BGC features \n")
            for row in antismash_rows:
            # Convert row list to tab-separated line
                out.write("\t".join(map(str, row)) + "\n")

        # Write Gecco features after GenBank features
        if gecco_rows:
            out.write("## Gecco predicted BGC features \n")
            for row in gecco_rows:
            # Convert row list to tab-separated line
                out.write(row + "\n")

    logger.info(
        f"GFF3 conversion complete: {feature_counts['CDS']} CDS, "
        f"{feature_counts['tRNA']} tRNA, {feature_counts['rRNA']} rRNA features written"
    )
    return gff_out
