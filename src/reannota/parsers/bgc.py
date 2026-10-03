import csv
import os
from collections import namedtuple

from Bio import SeqIO

# Named tuples for structure
# Added 'tool' and 'product' to BGC to keep track of prediction metadata
BGC = namedtuple("BGC", "contig_name bgc_name start end product tool orfs")
ORF = namedtuple("ORF", "locus_tag start end strand type product")

def load_antismash_gbk(file_path):
    """
    Load BGCs and ORFs from an antiSMASH GenBank (.gbk) file.
    """
    bgcs = []

    # Iterate over contigs in the GBK
    for record in SeqIO.parse(file_path, "genbank"):
        cluster_idx = 1

        for feature in record.features:
            if feature.type == "protocluster":
                # 1. Extract BGC info
                bgc_start = int(feature.location.start) + 1
                bgc_end = int(feature.location.end)
                product_list = feature.qualifiers.get("product", ["unknown"])
                bgc_product = ";".join(product_list)
                bgc_name = f"{record.id}_bgc{cluster_idx}"

                # 2. Extract ORFs within this BGC
                orfs = []
                for sub_feat in record.features:
                    if sub_feat.type == "CDS" and \
                       sub_feat.location.start >= feature.location.start and \
                       sub_feat.location.end <= feature.location.end:

                        # Get standard info
                        locus_tag = sub_feat.qualifiers.get("locus_tag", ["unknown"])[0]
                        prod = sub_feat.qualifiers.get("product", ["hypothetical protein"])[0]
                        strand = sub_feat.location.strand # 1 or -1

                        orfs.append(ORF(
                            locus_tag=locus_tag,
                            start=int(sub_feat.location.start) + 1,
                            end=int(sub_feat.location.end),
                            strand=strand,
                            type="CDS",
                            product=prod
                        ))

                bgcs.append(BGC(
                    contig_name=record.id,
                    bgc_name=bgc_name,
                    start=bgc_start,
                    end=bgc_end,
                    product=bgc_product,
                    tool="antiSMASH",
                    orfs=orfs
                ))
                cluster_idx += 1

    return bgcs

def load_gecco(tsv_path, gbk_paths_csv=None):
    """
    Load BGCs from GECCO output with coordinate correction and Deduplication.
    """
    bgcs = []
    gbk_lookup = {}

    # --- STRATEGY 1: Load explicit paths from CSV ---
    if gbk_paths_csv and os.path.exists(gbk_paths_csv):
        with open(gbk_paths_csv, 'r') as f:
            for line in f:
                clean_line = line.strip().split(',')[0]
                if clean_line and clean_line.endswith(".gbk"):
                    filename = os.path.basename(clean_line)
                    bgc_id_key = filename.replace(".gbk", "")
                    gbk_lookup[bgc_id_key] = clean_line

    # --- STRATEGY 2: Guess directory based on TSV location ---
    tsv_dir = os.path.dirname(os.path.abspath(tsv_path))
    
    with open(tsv_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')

        for row in reader:
            contig_id = row['sequence_id']
            bgc_id = row['cluster_id']
            bgc_type = row['type']
            bgc_start_global = int(row['start']) 
            bgc_end_global = int(row['end'])

            gbk_path = None
            if bgc_id in gbk_lookup:
                gbk_path = gbk_lookup[bgc_id]
            if not gbk_path:
                 potential_path = os.path.join(tsv_dir, f"{bgc_id}.gbk")
                 if os.path.exists(potential_path):
                     gbk_path = potential_path
            
            orfs = []
            if gbk_path and os.path.exists(gbk_path):
                offset = bgc_start_global - 1
                
                # --- NEW: Deduplication Set ---
                # We will store (start, end, strand) tuples here to check for duplicates
                seen_coords = set()

                try:
                    for record in SeqIO.parse(gbk_path, "genbank"):
                        for feature in record.features:
                            if feature.type == "CDS":
                                # Calculate coordinates first
                                global_orf_start = offset + int(feature.location.start) + 1
                                global_orf_end = offset + int(feature.location.end)
                                strand = feature.location.strand
                                
                                # Create a unique key for this gene location
                                coord_key = (global_orf_start, global_orf_end, strand)

                                # --- DEDUPLICATION CHECK ---
                                # If we have already seen a gene at this exact spot, skip this one
                                if coord_key in seen_coords:
                                    continue
                                
                                # Add to seen list so we don't add it again
                                seen_coords.add(coord_key)

                                l_tag = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                                
                                # Better Product Logic: Check product, then function, then fallback
                                if "product" in feature.qualifiers:
                                    prod = feature.qualifiers["product"][0]
                                elif "function" in feature.qualifiers:
                                    prod = feature.qualifiers["function"][0]
                                else:
                                    prod = "hypothetical protein"

                                orfs.append(ORF(
                                    locus_tag=l_tag,
                                    start=global_orf_start, 
                                    end=global_orf_end,
                                    strand=strand,
                                    type="CDS",
                                    product=prod
                                ))
                except Exception as e:
                    print(f"Error parsing GBK {gbk_path}: {e}")
            else:
                print(f"MISSING GBK: Could not find file for cluster '{bgc_id}'")

            bgcs.append(BGC(
                contig_name=contig_id,
                bgc_name=bgc_id,
                start=bgc_start_global,
                end=bgc_end_global,
                product=bgc_type,
                tool="GECCO",
                orfs=orfs
            ))

    return bgcs

def build_gff_rows(bgcs):
    """
    Generate GFF3 rows.
    - Column 2 (Source): Distinguishes the tool (antiSMASH vs GECCO).
    - Column 3 (Type): Always 'biosynthetic-gene-cluster'.
    - Column 9 (Attributes): Uses prefixes (as: or gecco:) to ensure unique IDs.
    """
    for bgc in bgcs:
        
        # 1. Source (Column 2): This tells you WHICH tool found it
        source_val = bgc.tool  # "antiSMASH" or "GECCO"
        
        # 2. Type (Column 3): Standardized as requested
        feature_type = "biosynthetic-gene-cluster"

        # 3. ID Prefixes: Keep unique IDs so tools don't merge them
        if bgc.tool == "antiSMASH":
            id_prefix = "as"
        else:
            id_prefix = "gecco"

        # e.g. ID=as:contig1_bgc1
        unique_id = f"{id_prefix}:{bgc.bgc_name}"

        # -------------------------
        # Parent BGC Entry
        # -------------------------
        attributes_bgc = [
            f"ID={unique_id}",
            f"product={bgc.product}",
            f"tool={bgc.tool}"  # Extra tag for clarity
        ]
        
        yield [
            bgc.contig_name,
            source_val,
            feature_type,      # <--- Now "biosynthetic-gene-cluster"
            str(bgc.start),
            str(bgc.end),
            ".",
            ".",
            ".",
            ";".join(attributes_bgc)
        ]
        
        # -------------------------
        # Child ORF Entries
        # -------------------------
        for orf in bgc.orfs:
            # Prefix ORF IDs too (as:locus_tag or gecco:locus_tag)
            orf_unique_id = f"{id_prefix}:{orf.locus_tag}"
            
            attributes_orf = [
                f"ID={orf_unique_id}",
                f"Parent={unique_id}",
                f"product={orf.product}"
            ]
            
            strand_char = "+" if orf.strand == 1 else "-"
            
            yield [
                bgc.contig_name,
                source_val,
                "CDS",
                str(orf.start),
                str(orf.end),
                ".",
                strand_char,
                "0",
                ";".join(attributes_orf)
            ]


