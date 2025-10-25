import csv
from pathlib import Path

from Bio import SeqIO


def read_file_paths(csv_path):
    """Read all file paths from a CSV. Return None if empty."""
    file_paths = []

    if not Path(csv_path).exists():
        return None

    with open(csv_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if row:
                path = row[0].strip()
                if path:
                    file_paths.append(path)

    return file_paths if file_paths else None

def combine_gbk_files(csv_path, output_file=None):
    """
    Combine GBK files listed in a CSV. 
    Returns combined records as a list in memory.
    Optionally writes to output_file if provided.
    Returns None if no valid GBK files found.
    """
    gbk_files = read_file_paths(csv_path)
    if not gbk_files:
        print("⚠️ No GECCO GBK files found in CSV. Skipping GECCO processing.")
        return None

    records = []
    for f in gbk_files:
        f_path = Path(f)
        if not f_path.exists():
            print(f"Warning: GECCO .gbk file not found — {f}")
            continue
        print(f"Reading {f} ...")
        records.extend(list(SeqIO.parse(f, "genbank")))

    if not records:
        print("⚠️ No valid GBK records found. Skipping GECCO processing.")
        return None

    if output_file:
        print(f"Writing combined GBK to {output_file} ...")
        SeqIO.write(records, output_file, "genbank")
        print(f"✅ Done! Combined {len(records)} records into {output_file}")

    return records

def extract_cluster_genes(gecco_file):
    """
    Parse a GECCO GBFF file and extract all genes in each cluster
    along with the cluster type.

    Parameters
    ----------
    gbff_file : str
        Path to the GBFF file (GenBank format).

    Returns
    -------
    genes_list : list of dict
        Each dict contains:
        - cluster_id
        - cluster_type
        - gene_name
        - product
        - start
        - end
        - strand
    """
    genes_list = []

    for record in SeqIO.parse(gecco_file, "genbank"):
        # Extract cluster type from COMMENT
        comment = record.annotations.get("comment", "")
        cluster_type = "Unknown"
        if comment:
            for line in comment.splitlines():
                if "cluster_type" in line and "::" in line:
                    cluster_type = line.split("::")[1].strip()
                    break

        # 2. If still unknown, try scanning the raw record text (fallback)
        if cluster_type == "Unknown":
            # Sometimes features or comments may be stored as text
            raw_text = record.format("genbank")
            for line in raw_text.splitlines():
                if "cluster_type" in line and "::" in line:
                    cluster_type = line.split("::")[1].strip()
                    break

        # Extract genes and CDS
        for feature in record.features:
            if feature.type in ["CDS"]:
                start = int(feature.location.start) + 1  # GFF 1-based
                end = int(feature.location.end)
                strand = "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "."
                gene_name = feature.qualifiers.get("locus_tag", ["N/A"])[0]
                product = feature.qualifiers.get("product", ["N/A"])[0]

                genes_list.append({
                    "cluster_id": record.id,
                    "cluster_type": cluster_type,
                    "gene_name": gene_name,
                    "product": product,
                    "start": start,
                    "end": end,
                    "strand": strand
                })

    return genes_list

#TODO change gecco input from gbk to csv and make function that take the csv
#  with file path and merges all gbk files into one that will be used for parsing

def build_gecco_gff_rows(genes_list):
    """
    Convert a list of gene dictionaries (from extract_cluster_genes)
    into GFF3-formatted rows.

    Parameters
    ----------
    genes_list : list of dict
        Output from extract_cluster_genes.

    Returns
    -------
    gff_lines : list of str
        Each string is a GFF3-formatted line.
    """
    gff_lines = []

    for gene in genes_list:

        if gene['gene_name'].startswith("contig_"):
            continue
        else:
        # Prepare attributes
            attributes = f"ID={gene['gene_name']};ClusterType={gene['cluster_type']};Product={gene['product']}"

            # Build GFF3 line
            line = (
                f"{gene['cluster_id']}\t"     # seqid
                f"GECCO\t"                    # source
                f"CDS\t"                     # type
                f"{gene['start']}\t"          # start
                f"{gene['end']}\t"            # end
                f".\t"                        # score
                f"{gene['strand']}\t"         # strand
                f".\t"                        # phase
                f"{attributes}"               # attributes
            )
            gff_lines.append(line)

    return gff_lines

