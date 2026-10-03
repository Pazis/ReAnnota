import collections
import logging

import matplotlib.pyplot as plt
from pycirclize import Circos


def gff_to_circos_png(gff_file, output_png, title=None):
    """
    Generate a Circos plot with:
    - Outer Rings: Genes & BGCs
    - Inner Track: Pseudogene Locations (Barcode style)
    """
    logging.info(f"Creating pycirclize plot for {gff_file}")

    # --- 1. Parse Data ---
    seq_lengths = {}

    # Outer Rings (Boxes)
    feats_fwd = collections.defaultdict(list)
    feats_rev = collections.defaultdict(list)
    feats_as  = collections.defaultdict(list)
    feats_ge  = collections.defaultdict(list)

    # Pseudogene Locations (Exact Coordinates)
    pseudo_coords = collections.defaultdict(list)

    # Valid types for the box tracks
    valid_func_types = ["tRNA", "rRNA", "tmRNA", "ncRNA"]

    with open(gff_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 9: continue

            raw_seqid, source, ftype, start, end, _, strand, _, attr = parts

            # Clean SeqID
            seqid = raw_seqid.replace("contig_", "c").replace("contig", "c")
            start, end = int(start), int(end)
            seq_lengths[seqid] = max(seq_lengths.get(seqid, 0), end)

            # --- PARSE BGCs (Outer Rings) ---
            is_bgc_row = ftype in ["biosynthetic-gene-cluster", "antiSMASH_region", "GECCO_region"]
            if is_bgc_row:
                tool = "AS" if ("antiSMASH" in source or "as:" in attr) else "GECCO"
                if tool == "AS": feats_as[seqid].append((start, end))
                else: feats_ge[seqid].append((start, end))
                continue

            # --- PARSE GENES ---
            if ftype == "biosynthetic-gene-cluster": continue

            attr_lower = attr.lower()
            is_hypo = "hypothetical protein" in attr_lower
            is_pseudo = "pseudogene" in ftype.lower() or "pseudogene" in attr_lower

            # 1. STORE PSEUDOGENE LOCATIONS
            # We store the tuple (start, end) so we can draw a bar of the correct width
            if is_pseudo:
                pseudo_coords[seqid].append((start, end))

            # 2. POPULATE BOXES (Outer Tracks)
            color = None
            if is_pseudo:
                color = "#000000"
            elif ftype == "CDS":
                color = "#FFFFFF00" if is_hypo else ("#1f77b4" if strand == "+" else "#2ca02c")
            elif ftype in valid_func_types:
                color = "#ff7f0e"

            if color:
                feature = (start, end, color)
                if strand == "+": feats_fwd[seqid].append(feature)
                else: feats_rev[seqid].append(feature)

    if not seq_lengths: return

   # --- ΠΡΟΣΘΗΚΗ: ΕΛΕΓΧΟΣ ΑΡΙΘΜΟΥ CONTIGS ΚΑΙ ΔΥΝΑΜΙΚΟ ΚΕΝΟ ---
    num_contigs = len(seq_lengths)

    # Αν τα contigs είναι πάνω από 40, σταματάμε εδώ
    if num_contigs > 40:
        logging.warning(f"Skipping Circos plot for {gff_file}: Found {num_contigs} contigs (limit is 40).")
        return

    # Αν συνεχίσουμε (<=40), υπολογίζουμε το κατάλληλο κενό (space)
    dynamic_space = min(5.0, 180.0 / num_contigs)
    # -----------------------------------------------------------

    # --- 2. Initialize Circos ---
    sorted_sectors = dict(sorted(seq_lengths.items(), key=lambda item: item[1], reverse=True))

    # --- ΠΡΟΣΘΗΚΗ: Περνάμε τη μεταβλητή dynamic_space αντί για το σταθερό "5" ---
    circos = Circos(sorted_sectors, space=dynamic_space)

    for sector in circos.sectors:
        seqid = sector.name
        length = sector.size

        # --- Track 1: Axis ---
        sector.text(f"{seqid}", size=12, r=110)
        track_axis = sector.add_track((98, 100))
        track_axis.axis(fc="none", ec="#333333", lw=0.5)
        track_axis.xticks_by_interval(500000, label_orientation="vertical", line_kws=dict(ec="black", lw=0.5), text_kws=dict(size=8))

        # --- Track 2 & 3: Genes (Boxes) ---
        track_fwd = sector.add_track((93, 97))
        track_fwd.axis(fc="#f0f8ff", ec="#cccccc", lw=0.3)
        for s, e, c in feats_fwd[seqid]: track_fwd.rect(s, e, fc=c, ec="none")

        track_rev = sector.add_track((87, 91))
        track_rev.axis(fc="#f0fff0", ec="#cccccc", lw=0.3)
        for s, e, c in feats_rev[seqid]: track_rev.rect(s, e, fc=c, ec="none")

        # --- Track 4 & 5: BGCs ---
        track_as = sector.add_track((78, 85))
        if feats_as[seqid]:
            track_as.axis(fc="#eaeaea", ec="none")
            for s, e in feats_as[seqid]: track_as.rect(s, e, fc="#9467bd", ec="black", lw=0.2)
        else: track_as.axis(fc="none", ec="none")

        track_ge = sector.add_track((70, 77))
        if feats_ge[seqid]:
            track_ge.axis(fc="#f9f9f9", ec="none")
            for s, e in feats_ge[seqid]: track_ge.rect(s, e, fc="#e377c2", ec="black", lw=0.2)
        else: track_ge.axis(fc="none", ec="none")

        # --- Track 6: PSEUDOGENE LOCATIONS (Barcode Style) ---
        # Create a narrow track for the ticks
        track_pseudo = sector.add_track((60, 65))
        track_pseudo.axis(fc="#f5f5f5", ec="none") # Very light grey background

        if pseudo_coords[seqid]:
            for start, end in pseudo_coords[seqid]:
                # Draw a rectangle for each pseudogene
                # Since pseudogenes can be small, we ensure a minimum width for visibility if needed,
                # but drawing exact start/end is most accurate.
                track_pseudo.rect(start, end, fc="black", ec="none")

    # --- 3. Legend & Save ---
    if title: circos.text(title, size=16, r=120)
    fig = circos.plotfig()

    handles = [
        plt.Rectangle((0,0),1,1, color="#1f77b4", label="Functional CDS (+)"),
        plt.Rectangle((0,0),1,1, color="#2ca02c", label="Functional CDS (-)"),
        plt.Rectangle((0,0),1,1, color="#000000", label="Pseudogene"),
        plt.Rectangle((0,0),1,1, color="#ff7f0e", label="RNA"),
        plt.Rectangle((0,0),1,1, color="#9467bd", label="antiSMASH"),
        plt.Rectangle((0,0),1,1, color="#e377c2", label="GECCO"),
    ]
    plt.legend(handles=handles, bbox_to_anchor=(0.5, 0.5), loc="center", frameon=False)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    logging.info(f"pycirclise plot saved to {output_png}")
