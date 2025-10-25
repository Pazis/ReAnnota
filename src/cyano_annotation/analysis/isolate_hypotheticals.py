"""Extract hypothetical proteins from GenBank files.

This module can be used as a standalone script or imported as a library.
"""

import argparse

from Bio import SeqIO


def extract_hypotheticals(input_file, output_file):
    """
    Extract hypothetical protein sequences from a GenBank file.

    Parameters:
    -----------
    input_file : str
        Path to input GenBank file.
    output_file : str
        Path to output FASTA file.

    Returns:
    --------
    int
        Number of hypothetical proteins extracted.
    """
    count = 0
    with open(output_file, "w") as out_f:
        # Parse the GenBank file sequence by sequence
        for x in SeqIO.parse(input_file, "genbank"):
            # Loop through all features in the sequence
            for feature in x.features:
                if feature.type == "CDS":
                    qualifiers = feature.qualifiers

                    # Only process features labeled as "hypothetical protein"
                    if "hypothetical protein" in qualifiers.get("product", [""])[0].lower():
                        # Ensure the feature has a protein translation
                        if "translation" in qualifiers:
                            prot_seq = qualifiers["translation"][0]
                            prot_id = qualifiers["locus_tag"][0]

                            # Compile relevant db_xref entries (COG, PFAM, EC, GO)
                            db_entries = ";".join(
                                [
                                    e
                                    for e in qualifiers.get("db_xref", [])
                                    if e.startswith(("COG:", "PFAM:", "EC:", "GO:"))
                                ]
                            )
                            header = f"{prot_id}" + (f" | {db_entries}" if db_entries else "")
                            out_f.write(f">{header}\n{prot_seq}\n")
                            count += 1
    return count


def main():
    """Command-line interface for extracting hypothetical proteins."""
    parser = argparse.ArgumentParser(description="Extract Hypothetical Proteins from Genbank file")
    parser.add_argument("-i", "--input", required=True, help="Input Genbank file")
    parser.add_argument("-o", "--output", required=True, help="Output fasta file")
    args = parser.parse_args()

    count = extract_hypotheticals(args.input, args.output)
    print(f"Extracted {count} hypothetical proteins to {args.output}")


if __name__ == "__main__":
    main()
