import collections
import logging

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pycircos


def gff_to_circos_png(
    gff_file, output_png, feature_types=None, title=None):
    """
    Generate a Circos-like genome visualization using PyCircos,
    with support for scatterplots, lineplots, barplots, heatmaps, and links.

    Parameters:
        gff_file (str): Input GFF file.
        output_png (str): Path to save the output figure.
        feature_types (list, optional): Feature types to visualize.
        title (str, optional): Title for the figure.
        links_file (str, optional): CSV file with link information.
    """
    Garc = pycircos.Garc
    Gcircle = pycircos.Gcircle

    if feature_types is None:
        feature_types = ["CDS", "tRNA", "rRNA", "pseudogene", "BGC"]

    color_map = {
        "CDS": "#0C5CA7",
        "tRNA": "#22C55E",
        "rRNA": "#F59E0B",
        "pseudogene": "#EF4444",
        "BGC": "#641CA8",
        "other": "#9CA3AF"
    }

    logging.info(f"🧬 Creating Circos plot for {gff_file}")

    # --- Create figure and Gcircle ---
    circle = Gcircle(figsize=(12, 8))

    seq_lengths = {}
    features_per_seq = collections.defaultdict(list)

    # --- Parse GFF ---
    with open(gff_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attributes = parts
            start, end = int(start), int(end)
            seq_lengths[seqid] = max(seq_lengths.get(seqid, 0), end)
            if ftype not in feature_types:
                continue
            feature_name = "NA"
            for attr in attributes.split(";"):
                if attr.startswith(("Name=", "gene=", "ID=")):
                    feature_name = attr.split("=")[1]
                    break
            features_per_seq[seqid].append({
                "type": ftype, "start": start, "end": end,
                "strand": strand, "name": feature_name,
                "color": color_map.get(ftype, color_map["other"])
            })

    if not seq_lengths:
        raise ValueError(f"No valid sequences found in {gff_file}")

    # --- Καθορισμός χρωμάτων ανά χρωμόσωμα ---
    chromosomes = list(seq_lengths.keys())
    num_chroms = len(chromosomes)
    # Επιλογή palette (π.χ. tab20)
    palette = cm.get_cmap("tab20", num_chroms)
    chrom_colors = {chrom: palette(i) for i, chrom in enumerate(chromosomes)}

    # --- Χρήση χρωμάτων όταν δημιουργούμε τα arcs ---
    for seqid, length in seq_lengths.items():
        arc = Garc(
            arc_id=seqid,
            size=length,
            interspace=2,
            raxis_range=(950,1000),
            facecolor=chrom_colors[seqid],   # διαφορετικό χρώμα ανά χρωμόσωμα
            edgecolor="#030B1B",
            linewidth=1.2,
            label=seqid,
            label_visible=True,
            labelposition=60
        )
        circle.add_garc(arc)

    circle.set_garcs(-65, 245)

    for text in circle.ax.texts:
        text.set_fontsize(5)

    # --- Feature lineplots & scatterplots ---
    for seqid, feats in features_per_seq.items():
        for f in feats:
            circle.lineplot(seqid, data=[1, 1], positions=[f["start"], f["end"]],
                            raxis_range=(850, 940), linecolor=f["color"], linewidth=2)
        starts = [f["start"] for f in feats]
        ends = [f["end"] for f in feats]
        colors = [f["color"] for f in feats]
        midpoints = [(s + e) / 2 for s, e in zip(starts, ends)]
        circle.scatterplot(seqid, data=[1] * len(midpoints), positions=midpoints,
                           facecolor=colors, raxis_range=(850, 940),
                           edgecolor="black", linewidth=0.1)

    # --- CDS density heatmap ---
    for seqid, feats in features_per_seq.items():
        cds_positions = [(f["start"], f["end"]) for f in feats if f["type"] == "CDS"]
        if cds_positions:
            arc = circle.garc_dict[seqid]
            density = arc.calc_density(cds_positions, window_size=10000)
            circle.heatmap(seqid, density, raxis_range=(770, 830), cmap=plt.cm.Blues)


    # --- Legend ---
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='CDS', markerfacecolor=color_map["CDS"], markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', label='tRNA', markerfacecolor=color_map["tRNA"], markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', label='rRNA', markerfacecolor=color_map["rRNA"], markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', label='Pseudogene', markerfacecolor=color_map["pseudogene"], markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', label='BGC', markerfacecolor=color_map["BGC"], markersize=8),
    ]
    circle.ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.25, 1),
                     title="Feature Types", fontsize=9, frameon=False)


    # --- Final touches & save ---
    circle.ax.set_facecolor("white")
    if title:
        circle.ax.set_title(title, fontsize=16, pad=20, fontweight="bold")
    circle.figure.savefig(output_png, dpi=600)

