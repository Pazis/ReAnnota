"""Convert GenBank (GBFF) files to GFF3 format."""

import logging

from Bio import SeqIO

from reannota.parsers import build_gff_rows , load_antismash_gbk, load_gecco

logger = logging.getLogger("ReAnnota")


def gbff_to_gff(gbff_in, gff_out, antismash_gbk=None, gecco_csv=None , gecco_summary=None):
    """
    Convert GenBank (GBFF) to GFF3 and optionally merge antiSMASH predictions.

    Parameters:
    - gbff_in: input GenBank file
    - gff_out: output GFF3 file
    - antismash_json: optional antiSMASH regions.json file
    - antismash_version: version string for antiSMASH features
    """
    logger.debug(f"Converting GBFF to GFF3: {gbff_in} → {gff_out}")

    feature_counts = {"CDS": 0, "tRNA": 0, "rRNA": 0 , "tmRNA": 0 , "ncRNA":0 , "repeat_region":0}

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
            # We pass the CSV path to the gbk_paths_csv argument we added
            ge_bgcs = load_gecco(str(gecco_summary), gbk_paths_csv=str(gecco_csv))
            all_bgcs.extend(ge_bgcs)
        except Exception as e:
            logger.error(f"Failed to load GECCO files: {e}")



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
                if feature.type in ["CDS", "tRNA", "rRNA", "tmRNA" , "ncRNA" , "repeat_region"]:
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
                    gff_line = f"{seqid}\tBiopython\t{feature.type}\t{start}\t{end}\t.\t{strand}\t.\t{attr_str}\n"
                    out.write(gff_line)

        # Write antiSMASH features after GenBank features
        if all_bgcs:
            out.write("## Predicted BGC features (antiSMASH / GECCO) \n")
            # Your new script yields lists of strings, we join them with tabs
            for row in build_gff_rows(all_bgcs):
                out.write("\t".join(row) + "\n")

    logger.info(
        f"GFF3 conversion complete: {feature_counts['CDS']} CDS, "
        f"{feature_counts['tRNA']} tRNA, {feature_counts['rRNA']} rRNA, {feature_counts['ncRNA']} ncRNA features written"
    )
    return gff_out
