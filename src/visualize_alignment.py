import os
import json
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from Bio import SeqIO
import numpy as np
from pathlib import Path

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.linewidth": 1.0,
        "pdf.fonttype": 42,
        "figure.dpi": 100,
        "savefig.dpi": 300,
    }
)

# Color scheme for bases
BASE_COLORS = {
    "A": "#2ecc71",  # Green
    "T": "#e74c3c",  # Red
    "G": "#f39c12",  # Orange
    "C": "#3498db",  # Blue
    "-": "#ecf0f1",  # Light gray for gaps
}

# TF colors (auto-generated from a palette)
TF_PALETTE = [
    "#FF6B6B",
    "#4ECDC4",
    "#45B7D1",
    "#FFA07A",
    "#98D8C8",
    "#F7DC6F",
    "#BB8FCE",
    "#85C1E2",
    "#F8B88B",
    "#ABEBC6",
    "#F5B7B1",
    "#AED6F1",
    "#F9E79F",
    "#D7BDE2",
    "#F1948A",
]


def load_alignment(fasta_file):
    """Load aligned sequences from FASTA file."""
    recs = list(SeqIO.parse(fasta_file, "fasta"))
    human = str(recs[0].seq)
    chimp = str(recs[1].seq)
    return human, chimp


def get_motif_width_from_meme(meme_file, motif_name):
    """Extract motif width from MEME file."""
    width = 10  # default
    try:
        with open(meme_file, "r") as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if f"NAME {motif_name}" in line or f"MOTIF {motif_name}" in line:
                    # Look for W (width) in next few lines
                    for j in range(i, min(i + 5, len(lines))):
                        if lines[j].startswith("W"):
                            width = int(lines[j].split()[1])
                            break
                    break
    except:
        pass
    return width


def map_position_ungapped_to_aligned(ungapped_pos, posmap):
    """
    Map ungapped position to aligned position using posmap.
    posmap: dict with string keys of ungapped positions -> aligned positions
    """
    pos_str = str(ungapped_pos)
    if pos_str in posmap:
        return posmap[pos_str]

    # Try nearby positions if exact not found
    for offset in range(-2, 3):
        alt_pos = str(ungapped_pos + offset)
        if alt_pos in posmap:
            return posmap[alt_pos]
    return None


def draw_aligned_sequences(ax, human_seq, chimp_seq, max_chars=100):
    """Draw aligned sequences with base coloring."""
    display_len = min(len(human_seq), max_chars)

    # Human sequence
    for i in range(display_len):
        base = human_seq[i]
        color = BASE_COLORS.get(base, "#95a5a6")
        rect = Rectangle(
            (i, 1), 1, 1, facecolor=color, edgecolor="white", linewidth=0.5
        )
        ax.add_patch(rect)

    # Chimp sequence
    for i in range(display_len):
        base = chimp_seq[i]
        color = BASE_COLORS.get(base, "#95a5a6")
        rect = Rectangle(
            (i, 0), 1, 1, facecolor=color, edgecolor="white", linewidth=0.5
        )
        ax.add_patch(rect)

    ax.set_xlim(0, display_len)
    ax.set_ylim(-0.5, 2.5)
    ax.set_yticks([0.5, 1.5])
    ax.set_yticklabels(["Chimp", "Human"], fontsize=9)
    ax.set_xticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    return display_len


def draw_tfbs_on_sequences(ax, hits, human_posmap, chimp_posmap, seq_length, tf_colors):
    """Draw TFBS hits on top of sequences."""

    human_hits = hits[hits["species"] == "human"]
    chimp_hits = hits[hits["species"] == "chimp"]

    y_offset_human = 1.35
    y_offset_chimp = 0.35

    # Draw human TFBS
    for _, hit in human_hits.iterrows():
        aln_pos = map_position_ungapped_to_aligned(hit["position"], human_posmap)
        if aln_pos is None or aln_pos >= seq_length:
            continue

        tf_name = hit["name"]
        color = tf_colors.get(tf_name, "#95a5a6")
        width = 5  # Motif width in aligned coords (approximate)

        rect = Rectangle(
            (aln_pos, y_offset_human),
            width,
            0.25,
            facecolor=color,
            alpha=0.8,
            edgecolor="black",
            linewidth=0.8,
        )
        ax.add_patch(rect)

    # Draw chimp TFBS
    for _, hit in chimp_hits.iterrows():
        aln_pos = map_position_ungapped_to_aligned(hit["position"], chimp_posmap)
        if aln_pos is None or aln_pos >= seq_length:
            continue

        tf_name = hit["name"]
        color = tf_colors.get(tf_name, "#95a5a6")
        width = 5

        rect = Rectangle(
            (aln_pos, y_offset_chimp),
            width,
            0.25,
            facecolor=color,
            alpha=0.8,
            edgecolor="black",
            linewidth=0.8,
        )
        ax.add_patch(rect)


def draw_conservation_track(ax, track, seq_length):
    """Draw conservation/substitution track."""

    colors_map = {
        "conserved": "#2ecc71",  # Green
        "substitution": "#e74c3c",  # Red
        "gap_human": "#f39c12",  # Orange
        "gap_chimp": "#f39c12",  # Orange
        "gap_both": "#ecf0f1",  # Light gray
    }

    for i, state in enumerate(track[:seq_length]):
        color = colors_map.get(state, "#95a5a6")
        rect = Rectangle(
            (i, 0), 1, 1, facecolor=color, edgecolor="white", linewidth=0.3
        )
        ax.add_patch(rect)

    ax.set_xlim(0, seq_length)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_ylabel("Conservation", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)


def plot_conservation_identity(ax, track, seq_length):
    """Plot conservation as a line graph."""

    conservation_vals = [1 if x == "conserved" else 0 for x in track]
    positions = np.arange(seq_length)

    # Moving average for smoother visualization
    window = max(1, seq_length // 100)
    if window > 1:
        smoothed = np.convolve(conservation_vals, np.ones(window) / window, mode="same")
    else:
        smoothed = conservation_vals

    ax.fill_between(positions, smoothed[:seq_length], alpha=0.4, color="#2ecc71")
    ax.plot(positions, smoothed[:seq_length], color="#27ae60", linewidth=1.5)

    ax.set_xlim(0, seq_length)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("Identity\n(proportion)", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", labelsize=8)


def draw_tfbs_summary(ax, hits, title=""):
    """Draw TFBS summary statistics."""

    human_hits = hits[hits["species"] == "human"]
    chimp_hits = hits[hits["species"] == "chimp"]

    human_tfs = human_hits.groupby("name").size().sort_values(ascending=False).head(10)
    chimp_tfs = chimp_hits.groupby("name").size().sort_values(ascending=False).head(10)

    all_tfs = sorted(set(human_tfs.index) | set(chimp_tfs.index))

    x_pos = np.arange(len(all_tfs))
    width = 0.35

    human_counts = [human_tfs.get(tf, 0) for tf in all_tfs]
    chimp_counts = [chimp_tfs.get(tf, 0) for tf in all_tfs]

    ax.bar(
        x_pos - width / 2,
        human_counts,
        width,
        label="Human",
        color="#3498db",
        alpha=0.8,
    )
    ax.bar(
        x_pos + width / 2,
        chimp_counts,
        width,
        label="Chimp",
        color="#e67e22",
        alpha=0.8,
    )

    ax.set_ylabel("Count", fontsize=9)
    ax.set_title(f"{title} - Top TFs", fontsize=10, fontweight="bold")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(all_tfs, rotation=45, ha="right", fontsize=8)
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_rewiring_barplot(ax, parquet_path, title=""):
    """
    Bar plot of binding outcomes per TF for one HAR.
    Reads directly from the HAR's parquet matrix file.
    """
    import pandas as pd
    from pathlib import Path

    if not Path(parquet_path).exists():
        ax.axis("off")
        ax.text(0.5, 0.5, "No matrix file found", transform=ax.transAxes, ha="center")
        return

    df = pd.read_parquet(parquet_path)
    if df.empty:
        ax.axis("off")
        return

    # Count outcomes per TF
    summary = (
        df.groupby(["tf_name", "binding_outcome"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=["shared", "human_only", "chimp_only"], fill_value=0)
    )
    summary = (
        summary[summary.sum(axis=1) > 0]
        .sort_values("human_only", ascending=False)
        .head(15)
    )

    x = np.arange(len(summary))
    w = 0.25
    ax.bar(
        x - w,
        summary["human_only"],
        w,
        label="Human only (gained)",
        color="#3498db",
        alpha=0.85,
    )
    ax.bar(x, summary["shared"], w, label="Shared", color="#2ecc71", alpha=0.85)
    ax.bar(
        x + w,
        summary["chimp_only"],
        w,
        label="Chimp only (lost)",
        color="#e67e22",
        alpha=0.85,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(summary.index, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Positions", fontsize=9)
    ax.set_title(f"{title} — TF Binding Outcomes", fontsize=10, fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def make_figure(
    alignment_fasta,
    posmap_json,
    conservation_json,
    tfbs_tsv,
    output,
    meme_file=None,
    har_name=None,
    max_chars=None,
):
    """Generate comprehensive alignment visualization."""

    human_seq, chimp_seq = load_alignment(alignment_fasta)

    with open(posmap_json) as f:
        maps = json.load(f)

    with open(conservation_json) as f:
        cons = json.load(f)

    hits_df = pd.read_csv(tfbs_tsv, sep="\t", low_memory=False)

    human_region = maps["human_region"]
    hits = hits_df[hits_df["coords_hg38"] == human_region].copy()

    # Use HAR name for titles if provided, otherwise fall back to coordinates
    display_name = har_name if har_name else human_region

    # If max_chars not set, use full sequence length
    seq_display_len = max_chars if max_chars else len(human_seq)

    # Create TF color mapping
    all_tfs = sorted(hits["name"].unique())
    tf_colors = {tf: TF_PALETTE[i % len(TF_PALETTE)] for i, tf in enumerate(all_tfs)}

    # Create figure with multiple subplots
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(4, 2, hspace=0.35, wspace=0.3)

    # Main alignment view
    ax_seqs = fig.add_subplot(gs[0:2, :])
    seq_len = draw_aligned_sequences(
        ax_seqs, human_seq, chimp_seq, max_chars=seq_display_len
    )
    ax_seqs.set_title(
        f"Pairwise Alignment with TFBS Mapping\n{display_name}",
        fontsize=12,
        fontweight="bold",
        pad=10,
    )

    # Draw TFBS overlays
    draw_tfbs_on_sequences(
        ax_seqs,
        hits,
        {str(k): v for k, v in maps["human_pos_map"].items()},
        {str(k): v for k, v in maps["chimp_pos_map"].items()},
        seq_len,
        tf_colors,
    )

    # Conservation track
    ax_cons = fig.add_subplot(gs[2, :])
    draw_conservation_track(ax_cons, cons["track"], seq_len)

    # Identity plot
    ax_id = fig.add_subplot(gs[3, 0])
    plot_conservation_identity(ax_id, cons["track"], seq_len)

    # TF summary
    ax_rewiring = fig.add_subplot(gs[3, 1])
    # derive parquet path from the alignment fasta path
    parquet_path = alignment_fasta.replace(".aligned.fasta", ".parquet").replace(
        "alignments/", "tf_matrices/"
    )
    draw_rewiring_barplot(ax_rewiring, parquet_path, title=display_name)
    # Add legend for conservation states
    legend_elements = [
        mpatches.Patch(facecolor="#2ecc71", label="Conserved (match)"),
        mpatches.Patch(facecolor="#e74c3c", label="Substitution"),
        mpatches.Patch(facecolor="#f39c12", label="Gap"),
        mpatches.Patch(facecolor="#ecf0f1", label="Gap (both)"),
    ]

    # TF legend
    tf_legend_elements = [
        mpatches.Patch(
            facecolor=tf_colors[tf], label=tf, edgecolor="black", linewidth=0.5
        )
        for tf in sorted(tf_colors.keys())[:15]  # Limit to top 15 TFs
    ]

    fig.legend(
        handles=legend_elements,
        loc="upper left",
        fontsize=8,
        bbox_to_anchor=(0.01, 0.98),
        frameon=True,
    )

    fig.suptitle(
        f"HAR Comparative Analysis: Human vs Chimp\n{display_name}",
        fontsize=14,
        fontweight="bold",
        y=0.995,
    )

    plt.savefig(output, dpi=300, bbox_inches="tight", format=output.split(".")[-1])
    print(f"Saved figure: {output}")
    plt.close()


def batch_visualize(alignments_dir, tfbs_tsv, output_dir, meme_file=None):
    """Generate visualizations for all alignments in a directory."""

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Find all alignment FASTA files
    alignment_files = sorted(Path(alignments_dir).glob("*.aligned.fasta"))

    if not alignment_files:
        print(f"No alignment files found in {alignments_dir}")
        return

    print(f"Found {len(alignment_files)} alignment files")

    for fasta_file in alignment_files:
        base_name = fasta_file.stem.replace(".aligned", "")

        # Find corresponding posmap and conservation files
        posmap_file = fasta_file.parent / f"{base_name}.posmaps.json"
        conservation_file = fasta_file.parent / f"{base_name}.conservation.json"

        if not posmap_file.exists() or not conservation_file.exists():
            print(f"Skipping {base_name}: missing posmap or conservation files")
            continue

        output_file = Path(output_dir) / f"{base_name}_visualization.pdf"

        try:
            make_figure(
                str(fasta_file),
                str(posmap_file),
                str(conservation_file),
                tfbs_tsv,
                str(output_file),
                meme_file,
            )
        except Exception as e:
            print(f"Error processing {base_name}: {e}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Visualize alignments with TFBS mapping and conservation analysis"
    )

    # Single figure mode
    parser.add_argument("--alignment", help="Path to aligned FASTA file")
    parser.add_argument("--posmap", help="Path to posmap JSON file")
    parser.add_argument("--conservation", help="Path to conservation JSON file")
    parser.add_argument(
        "--tfbs", required=True, help="Path to TFBS hits TSV from scan.py"
    )
    parser.add_argument("--output", help="Output PDF/PNG file")
    parser.add_argument(
        "--meme", default=None, help="Path to MEME file for motif widths"
    )
    parser.add_argument(
        "--har_name", default=None, help="HAR name for figure titles (e.g. HAR_202)"
    )
    parser.add_argument(
        "--max_chars",
        type=int,
        default=None,
        help="Max alignment positions to display (default: full length)",
    )

    # Batch mode
    parser.add_argument("--alignments_dir", help="Directory with all alignment files")
    parser.add_argument(
        "--output_dir", help="Output directory for batch visualizations"
    )

    args = parser.parse_args()

    if args.alignment and args.posmap and args.conservation and args.output:
        # Single figure mode
        make_figure(
            args.alignment,
            args.posmap,
            args.conservation,
            args.tfbs,
            args.output,
            args.meme,
            args.har_name,
            args.max_chars,
        )
    elif args.alignments_dir and args.output_dir:
        # Batch mode
        batch_visualize(args.alignments_dir, args.tfbs, args.output_dir, args.meme)
    else:
        parser.print_help()
