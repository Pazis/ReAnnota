import csv
import re
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

        if cluster_type == "Unknown":
            raw_text = record.format("genbank")
            for line in raw_text.splitlines():
                if "cluster_type" in line and "::" in line:
                    cluster_type = line.split("::")[1].strip()
                    break

        for feature in record.features:
            if feature.type == "CDS":
                start = int(feature.location.start) + 1
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


def build_gecco_gff_rows(genes_list):
    """
    Convert a flat list of genes (from extract_cluster_genes)
    into GFF3-formatted rows with one biosynthetic cluster per record.
    """
    gff_lines = []

    # --- Group genes by cluster_id ---
    clusters = {}
    for g in genes_list:
        if g["gene_name"].startswith("contig_"):
            continue  # skip placeholder genes
        cid = g["cluster_id"]
        if cid not in clusters:
            clusters[cid] = {"cluster_type": g["cluster_type"], "genes": []}
        clusters[cid]["genes"].append(g)

    # --- Build GFF rows ---
    for i, (cluster_id, data) in enumerate(clusters.items(), start=1):
        match = re.search(r"(contig_\d+)", cluster_id)
        seqid = match.group(1) if match else cluster_id

        cluster_start = min(g["start"] for g in data["genes"])
        cluster_end = max(g["end"] for g in data["genes"])

        cluster_tag = f"{seqid}_bgc{i}"
        cluster_line = (
            f"{seqid}\t"
            f"GECCO\t"
            f"biosynthetic-gene-cluster\t"
            f"{cluster_start}\t{cluster_end}\t"
            f".\t.\t.\t"
            f"ID={cluster_tag}"
        )
        gff_lines.append(cluster_line)

        for gene in data["genes"]:
            # Clean seqid again for CDS lines
            gmatch = re.search(r"(contig_\d+)", gene['cluster_id'])
            gseqid = gmatch.group(1) if gmatch else gene['cluster_id']

            attributes = (
                f"ID={gene['gene_name']};"
                f"Parent={cluster_tag};"
                f"ClusterType={gene['cluster_type']};"
                f"Product={gene['product']}"
            )

            gene_line = (
                f"{gseqid}\t"      
                f"GECCO\t"
                f"CDS\t"
                f"{gene['start']}\t{gene['end']}\t"
                f".\t{gene['strand']}\t.\t"
                f"{attributes}"
            )
            gff_lines.append(gene_line)

    return gff_lines


