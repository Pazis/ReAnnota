# antismash_gff_utils.py

import json
import re
from collections import namedtuple

# Named tuples for structure
BGC = namedtuple("BGC", "contig_name bgc_name start end orfs")
ORF = namedtuple("ORF", "locus_tag start end strand type product")


def load_regions_json(file_path):
    """
    Load BGCs and ORFs from antiSMASH JSON or JS file.
    Supports:
      - JSON file (regions.json)
      - JS file (regions.js) without js2py
    """
    if file_path.endswith(".js"):
        with open(file_path) as f:
            js_text = f.read()

        # Assume format: var recordData = [...]
        match = re.search(r"recordData\s*=\s*(\[.*\]);", js_text, re.DOTALL)
        if not match:
            raise ValueError("Could not find recordData array in JS file")
        json_str = match.group(1)

        # Χρήση demjson3 για robust parsing
        import demjson3
        json_data = demjson3.decode(json_str)


    else:
        # Assume JSON
        with open(file_path) as f:
            json_data = json.load(f)

    # Parse JSON into BGC and ORF namedtuples
    bgcs = []
    for entry in json_data:
        regions = entry.get("regions", [])
        for region in regions:
            region_start = region["start"]
            region_end = region["end"]
            index = region["idx"]
            orfs = []
            for orf in region.get("orfs", []):
                orfs.append(ORF(
                    locus_tag=orf.get("locus_tag", ""),
                    start=orf.get("start", 0),
                    end=orf.get("end", 0),
                    strand=int(orf.get("strand", 1)),
                    type=orf.get("type", ""),
                    product=orf.get("product", "")
                ))
            bgcs.append(BGC(
                contig_name=entry.get("seq_id", "unknown_contig"),
                bgc_name=f"{entry.get('seq_id', 'contig')}_bgc{index}",
                start=region_start,
                end=region_end,
                orfs=orfs
            ))
    return bgcs


def build_gff_rows(bgcs, antismash_version):
    """Generate GFF3 rows for BGCs and ORFs"""
    for bgc in bgcs:
        # Parent BGC feature
        yield [
            bgc.contig_name,
            f"antiSMASH:{antismash_version}",
            "biosynthetic-gene-cluster",
            bgc.start,
            bgc.end,
            ".",
            ".",
            ".",
            f"ID={bgc.bgc_name}"
        ]
        for orf in bgc.orfs:
            attrs = [f"ID={orf.locus_tag}", f"Parent={bgc.bgc_name}"]
            if orf.product:
                attrs.append(f"product={orf.product}")
            if orf.type:
                attrs.append(f"function={orf.type}")
            yield [
                bgc.contig_name,
                f"antiSMASH:{antismash_version}",
                "CDS",
                orf.start,
                orf.end,
                ".",
                "+" if orf.strand == 1 else "-",
                ".",
                ";".join(attrs)
            ]
