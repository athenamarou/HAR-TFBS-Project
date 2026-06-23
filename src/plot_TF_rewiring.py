"""
Reads rewiring_report.txt (output of summarize_rewiring.py)
and produces publication-ready bar charts for:
    - Top N TFs gained in human (human_only binding)
    - Top N TFs lost in human  (chimp_only binding)
    - A combined side-by-side figure

Output files (in output_dir/):
    - tf_gained_in_human.pdf
    - tf_lost_in_human.pdf
    - tf_rewiring_combined.pdf

Usage:
    python plot_tf_rewiring.py --report data/rewiring_report/rewiring_report.txt
    python plot_tf_rewiring.py --report data/rewiring_report/rewiring_report.txt --topn 15
    python plot_tf_rewiring.py --report data/rewiring_report/rewiring_report.txt --output_dir figures/
"""

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

_ROOT = Path(__file__).resolve().parent.parent
_RESULTS = _ROOT / "results"

# ---------------------------------------------------------------------------
# TF family classification
# ---------------------------------------------------------------------------

TF_FAMILY = {
    # HOX cluster
    "HOXA1": "HOX",
    "HOXA2": "HOX",
    "HOXA7": "HOX",
    "HOXB2": "HOX",
    "HOXB3": "HOX",
    "HOXB5": "HOX",
    "HOXC8": "HOX",
    # Other homeodomain
    "VAX2": "Homeodomain",
    "BSX": "Homeodomain",
    "UNCX": "Homeodomain",
    "GSX1": "Homeodomain",
    "NKX6-2": "Homeodomain",
    "NKX6-1": "Homeodomain",
    "DLX6": "Homeodomain",
    "MNX1": "Homeodomain",
    "ALX3": "Homeodomain",
    "PDX1": "Homeodomain",
    "RAX2": "Homeodomain",
    "EVX1": "Homeodomain",
    "EVX2": "Homeodomain",
    "HOXA7": "HOX",
    # LIM-homeodomain / Iroquois
    "Lhx4": "LIM-HD",
    "SHOX": "LIM-HD",
    "IRX5": "LIM-HD",
    # Zinc finger
    "ZGLP1": "Zinc finger",
    # TBP / basal
    "TBP": "TBP/Basal",
}

FAMILY_COLORS = {
    "HOX": "#378ADD",
    "Homeodomain": "#7F77DD",
    "LIM-HD": "#1D9E75",
    "Zinc finger": "#BA7517",
    "TBP/Basal": "#D85A30",
    "Other": "#888780",
}


def get_family(tf_name: str) -> str:
    return TF_FAMILY.get(tf_name, "Other")


# ---------------------------------------------------------------------------
# Parse rewiring_report.txt
# ---------------------------------------------------------------------------


def parse_rewiring_report(report_path: str):
    """
    Extract gained and lost TF tables from the text report.
    Returns two lists of dicts: gained, lost
    Each dict: {tf, positions, sub_driven, n_hars, family}
    """
    gained = []
    lost = []

    # Pattern
    pattern = re.compile(
        r"^\s{2}(\S+)\s+positions:\s+(\d+)\s+substitution-driven:\s+(\d+)\s+HARs:\s+(\d+)"
    )

    current_section = None

    with open(report_path) as fh:
        for line in fh:
            if "TFBS GAINED IN HUMAN" in line:
                current_section = "gained"
                continue
            if "TFBS LOST IN HUMAN" in line:
                current_section = "lost"
                continue
            if "TOP HARs BY REWIRING" in line:
                current_section = None
                continue

            m = pattern.match(line)
            if m and current_section:
                tf, pos, sub, hars = (
                    m.group(1),
                    int(m.group(2)),
                    int(m.group(3)),
                    int(m.group(4)),
                )
                entry = {
                    "tf": tf,
                    "positions": pos,
                    "sub_driven": sub,
                    "n_hars": hars,
                    "family": get_family(tf),
                }
                if current_section == "gained":
                    gained.append(entry)
                else:
                    lost.append(entry)

    return gained, lost


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def draw_horizontal_bars(
    ax, data, topn, color_by="family", title="", xlabel="Rewiring positions"
):
    """
    Draw a horizontal bar chart on ax for the top N TFs.
    Each bar is split into substitution-driven (darker) vs other (lighter).
    """
    data = data[:topn]
    max_pos = data[0]["positions"] if data else 1

    y_positions = list(range(len(data)))
    tf_labels = [d["tf"] for d in data]
    colors = [FAMILY_COLORS.get(d["family"], "#888780") for d in data]

    for i, (d, color) in enumerate(zip(data, colors)):
        total = d["positions"]
        sub = d["sub_driven"]
        other = total - sub

        # Full bar (lighter — total)
        ax.barh(i, total, color=color, alpha=0.35, height=0.65)
        # Substitution-driven portion (solid)
        ax.barh(i, sub, color=color, alpha=1.0, height=0.65)

        # Annotation: total positions + n_hars
        ax.text(
            total + max_pos * 0.01,
            i,
            f"{total:,}  ({d['n_hars']} HARs)",
            va="center",
            ha="left",
            fontsize=8,
            color="#444441",
        )

    ax.set_yticks(y_positions)
    ax.set_yticklabels(tf_labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_title(title, fontsize=11, fontweight="normal", pad=8)
    ax.set_xlim(0, max_pos * 1.3)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", linestyle="--", alpha=0.35, linewidth=0.5)

    # Family color legend
    families_present = list(dict.fromkeys(d["family"] for d in data))
    legend_patches = [
        mpatches.Patch(color=FAMILY_COLORS.get(f, "#888780"), label=f)
        for f in families_present
    ]
    # Sub-driven legend
    legend_patches += [
        mpatches.Patch(facecolor="gray", alpha=1.0, label="Substitution-driven"),
        mpatches.Patch(facecolor="gray", alpha=0.35, label="Other (gap/conserved)"),
    ]
    ax.legend(
        handles=legend_patches,
        fontsize=7,
        frameon=False,
        loc="lower right",
        ncol=2,
    )


def plot_single(data, topn, title, xlabel, output_path):
    fig, ax = plt.subplots(figsize=(9, max(4, topn * 0.45 + 1.5)))
    draw_horizontal_bars(ax, data, topn, title=title, xlabel=xlabel)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


def plot_combined(gained, lost, topn, output_path):
    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(17, max(5, topn * 0.45 + 2)),
        gridspec_kw={"wspace": 0.55},
    )
    draw_horizontal_bars(
        ax1,
        gained,
        topn,
        title="TFBS gained in human\n(human-only binding)",
        xlabel="Rewiring positions",
    )
    draw_horizontal_bars(
        ax2,
        lost,
        topn,
        title="TFBS lost in human\n(chimp-only binding)",
        xlabel="Rewiring positions",
    )
    fig.suptitle(
        "Transcription Factor Binding Rewiring — Human vs Chimpanzee HARs",
        fontsize=13,
        fontweight="normal",
        y=1.01,
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Plot top gained/lost TFs from rewiring_report.txt"
    )
    parser.add_argument(
        "--report",
        required=True,
        help="Path to rewiring_report.txt (output of summarize_rewiring.py)",
    )
    parser.add_argument(
        "--topn",
        type=int,
        default=10,
        help="How many top TFs to plot (default: 10)",
    )
    parser.add_argument(
        "--output_dir",
        default=str(_RESULTS / "figures"),
        help="Output directory for PDFs (default: results/figures/)",
    )
    args = parser.parse_args()

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("TF REWIRING BAR CHART GENERATOR")
    print("=" * 60)

    print(f"\nParsing {args.report} ...")
    gained, lost = parse_rewiring_report(args.report)
    print(f"  Found {len(gained)} gained TFs, {len(lost)} lost TFs")

    print(f"\nPlotting top {args.topn} TFs...")

    plot_single(
        gained,
        args.topn,
        title=f"Top {args.topn} TFs gained in human (human-only binding)",
        xlabel="Rewiring positions",
        output_path=f"{args.output_dir}/tf_gained_in_human.pdf",
    )

    plot_single(
        lost,
        args.topn,
        title=f"Top {args.topn} TFs lost in human (chimp-only binding)",
        xlabel="Rewiring positions",
        output_path=f"{args.output_dir}/tf_lost_in_human.pdf",
    )

    plot_combined(
        gained,
        lost,
        args.topn,
        output_path=f"{args.output_dir}/tf_rewiring_combined.pdf",
    )

    print("\n" + "=" * 60)
    print(f"Done. Figures in {args.output_dir}/")
    print("=" * 60)


if __name__ == "__main__":
    main()
