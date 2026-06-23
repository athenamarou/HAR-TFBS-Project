"""

Classifies each HAR into a genomic context category by intersecting
HAR coordinates (hg38 BED) with GENCODE v38 gene annotations.

Categories:
    1. promoter_proximal  : within 2,000 bp of a TSS
    2. exonic             : overlaps a protein-coding exon
    3. intronic           : falls inside a gene body (intron)
    4. intergenic_proximal: intergenic but within 500 kb of the nearest gene
    5. intergenic_distal  : intergenic and > 500 kb from any gene

Output files (in output_dir/):
    - har_categories.tsv            : one row per HAR with category + nearest gene
    - category_summary.tsv          : count and % per category
    - har_category_barchart.pdf     : horizontal bar chart of category counts
    - har_category_piechart.pdf     : pie chart of proportions

Requirements:
    pip install pybedtools pandas matplotlib requests

GENCODE GTF is downloaded automatically on first run and cached locally.
If you already have it, pass --gtf path/to/gencode.v38.annotation.gtf.gz

Commands:
    python categorize_hars.py --bed data/hars_hg38.bed --output_dir data/har_categories/
    python categorize_hars.py --bed data/hars_hg38.bed --gtf data/gencode.v38.annotation.gtf.gz --output_dir data/har_categories/
    python categorize_hars.py --bed data/hars_hg38.bed --promoter_window 1000 --distal_window 250000 #if we have specific difinitions in mind
"""

import argparse
import os
import sys
import subprocess
import urllib.request
from pathlib import Path

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

_ROOT = Path(__file__).resolve().parent.parent
_DATA = _ROOT / "data"
_RESULTS = _ROOT / "results"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

GENCODE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_38/gencode.v38.annotation.gtf.gz"
)
DEFAULT_GTF_CACHE = str(_DATA / "gencode.v38.annotation.gtf.gz")

CATEGORY_ORDER = [
    "promoter_proximal",
    "exonic",
    "intronic",
    "intergenic_proximal",
    "intergenic_distal",
]

CATEGORY_COLORS = {
    "promoter_proximal": "#378ADD",
    "exonic": "#D85A30",
    "intronic": "#7F77DD",
    "intergenic_proximal": "#1D9E75",
    "intergenic_distal": "#888780",
}

CATEGORY_LABELS = {
    "promoter_proximal": "Promoter-proximal (≤2 kb from TSS)",
    "exonic": "Exonic (coding overlap)",
    "intronic": "Intronic (gene body)",
    "intergenic_proximal": "Intergenic proximal (≤500 kb from gene)",
    "intergenic_distal": "Intergenic distal (>500 kb from gene)",
}


# ---------------------------------------------------------------------------
# Obtain and parse GENCODE GTF
# ---------------------------------------------------------------------------


def download_gtf(dest_path: str) -> None:
    """Download GENCODE v38 GTF if not already present."""
    if os.path.exists(dest_path):
        print(f"  Using cached GTF: {dest_path}")
        return
    print(f"  Downloading GENCODE v38 GTF (~1.5 GB uncompressed, ~50 MB gz)...")
    print(f"  URL: {GENCODE_URL}")
    print("  This may take a few minutes on first run.")
    urllib.request.urlretrieve(GENCODE_URL, dest_path)
    print(f"  Saved to {dest_path}")


def gtf_to_beds(gtf_path: str, work_dir: str, promoter_window: int = 2000):
    """
    Parse the GTF and write three BED files:
        - tss.bed         : TSS +- promoter_window  (for promoter_proximal)
        - exons.bed       : all protein-coding exons
        - gene_bodies.bed : full gene spans (for intronic)
        - genes_for_dist.bed : gene bodies extended by 0 (for distance calc)

    """
    import gzip

    tss_path = os.path.join(work_dir, "tss.bed")
    exon_path = os.path.join(work_dir, "exons.bed")
    gene_path = os.path.join(work_dir, "gene_bodies.bed")

    # Skip re-parsing if already done
    if all(os.path.exists(p) for p in [tss_path, exon_path, gene_path]):
        print("  Annotation BED files already built, reusing.")
        return {"tss": tss_path, "exons": exon_path, "genes": gene_path}

    print("  Parsing GTF — extracting TSS, exons, gene bodies...")

    opener = gzip.open if gtf_path.endswith(".gz") else open

    tss_rows = []
    exon_rows = []
    gene_rows = []

    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs = fields
            start = int(start) - 1  # GTF is 1-based → BED is 0-based
            end = int(end)

            # Extract gene_name from attributes string
            gene_name = "."
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_name"):
                    gene_name = attr.split('"')[1] if '"' in attr else attr.split()[1]
                    break

            if feature == "gene":
                # TSS: start of gene on + strand, end on - strand
                if strand == "+":
                    tss_start = max(0, start - promoter_window)
                    tss_end = start + promoter_window
                else:
                    tss_start = max(0, end - promoter_window)
                    tss_end = end + promoter_window

                tss_rows.append(
                    f"{chrom}\t{tss_start}\t{tss_end}\t{gene_name}\t.\t{strand}\n"
                )
                gene_rows.append(f"{chrom}\t{start}\t{end}\t{gene_name}\t.\t{strand}\n")

            elif feature == "exon":
                # Keep only protein-coding exons to reduce false overlaps
                if "protein_coding" in attrs:
                    exon_rows.append(
                        f"{chrom}\t{start}\t{end}\t{gene_name}\t.\t{strand}\n"
                    )

    print(
        f"  Writing {len(tss_rows):,} TSS windows, {len(exon_rows):,} exons, {len(gene_rows):,} gene bodies"
    )

    for path, rows in [
        (tss_path, tss_rows),
        (exon_path, exon_rows),
        (gene_path, gene_rows),
    ]:
        with open(path, "w") as fout:
            fout.writelines(rows)

    return {"tss": tss_path, "exons": exon_path, "genes": gene_path}


# ---------------------------------------------------------------------------
# BED intersections via pybedtools
# ---------------------------------------------------------------------------


def intersect_bed(query_bed: str, subject_bed: str, flags: list = None) -> set:
    """
    Run bedtools intersect and return set of query intervals (chrom:start-end)
    that overlap the subject. Uses -u (report each query once).
    """
    try:
        import pybedtools
    except ImportError:
        raise ImportError("pybedtools is required: pip install pybedtools")

    a = pybedtools.BedTool(query_bed)
    b = pybedtools.BedTool(subject_bed)
    result = a.intersect(b, u=True)
    return {f"{f.chrom}:{f.start}-{f.end}" for f in result}


def nearest_gene_distance(query_bed: str, genes_bed: str) -> dict:
    """
    For each HAR, find the distance to the nearest gene body using bedtools closest.
    Returns dict: {coord_key: (distance, gene_name)}
    """
    try:
        import pybedtools
    except ImportError:
        raise ImportError("pybedtools is required: pip install pybedtools")

    a = pybedtools.BedTool(query_bed)
    b = pybedtools.BedTool(genes_bed).sort()
    a_sorted = a.sort()

    # -d reports distance in last column; -t first reports only 1 hit per query
    closest = a_sorted.closest(b, d=True, t="first")

    result = {}
    for f in closest:
        key = f"{f.chrom}:{f.start}-{f.end}"
        # fields: chrom,start,end,[name],chrom2,start2,end2,gene_name,...,distance
        fields = str(f).rstrip().split("\t")
        try:
            distance = int(fields[-1])
            gene_name = fields[-4] if len(fields) >= 4 else "."
        except (ValueError, IndexError):
            distance = -1
            gene_name = "."
        result[key] = (distance, gene_name)

    return result


# ---------------------------------------------------------------------------
# Classify each HAR
# ---------------------------------------------------------------------------


def classify_hars(
    har_bed: str,
    beds: dict,
    distal_window: int = 500_000,
) -> pd.DataFrame:
    """
    Assign each HAR to exactly one category (priority order).
    Returns a DataFrame with columns:
        chrom, start, end, har_coord, category, nearest_gene, distance_to_gene
    """
    print("  Computing BED intersections...")

    in_promoter = intersect_bed(har_bed, beds["tss"])
    in_exon = intersect_bed(har_bed, beds["exons"])
    in_gene = intersect_bed(har_bed, beds["genes"])

    print("  Computing distances to nearest gene...")
    distances = nearest_gene_distance(har_bed, beds["genes"])

    rows = []
    with open(har_bed) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            key = f"{chrom}:{start}-{end}"

            dist, nearest_gene = distances.get(key, (-1, "."))

            # Priority-based assignment
            if key in in_promoter:
                category = "promoter_proximal"
            elif key in in_exon:
                category = "exonic"
            elif key in in_gene:
                category = "intronic"
            elif 0 <= dist <= distal_window:
                category = "intergenic_proximal"
            else:
                category = "intergenic_distal"

            rows.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "har_coord": key,
                    "category": category,
                    "nearest_gene": nearest_gene,
                    "distance_to_gene": dist,
                }
            )

    df = pd.DataFrame(rows)
    print(f"  Classified {len(df):,} HARs")
    return df


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------


def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    total = len(df)
    summary = (
        df["category"]
        .value_counts()
        .reindex(CATEGORY_ORDER, fill_value=0)
        .reset_index()
    )
    summary.columns = ["category", "count"]
    summary["percent"] = (summary["count"] / total * 100).round(1)
    summary["label"] = summary["category"].map(CATEGORY_LABELS)
    return summary


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------


def plot_barchart(summary: pd.DataFrame, output_path: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 4))

    colors = [CATEGORY_COLORS[c] for c in summary["category"]]
    bars = ax.barh(
        summary["label"],
        summary["count"],
        color=colors,
        edgecolor="white",
        linewidth=0.5,
        height=0.6,
    )

    # Annotate bars with count + %
    for bar, (_, row) in zip(bars, summary.iterrows()):
        w = bar.get_width()
        ax.text(
            w + summary["count"].max() * 0.01,
            bar.get_y() + bar.get_height() / 2,
            f"{int(row['count']):,}  ({row['percent']}%)",
            va="center",
            ha="left",
            fontsize=9,
            color="#444441",
        )

    ax.set_xlabel("Number of HARs", fontsize=10)
    ax.set_title(
        "HAR genomic context classification", fontsize=12, fontweight="normal", pad=10
    )
    ax.invert_yaxis()
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, summary["count"].max() * 1.25)
    ax.grid(axis="x", linestyle="--", alpha=0.4, linewidth=0.5)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Bar chart saved: {output_path}")


def plot_piechart(summary: pd.DataFrame, output_path: str) -> None:
    nonzero = summary[summary["count"] > 0]
    colors = [CATEGORY_COLORS[c] for c in nonzero["category"]]

    fig, ax = plt.subplots(figsize=(7, 5))
    wedges, texts, autotexts = ax.pie(
        nonzero["count"],
        labels=None,
        colors=colors,
        autopct=lambda p: f"{p:.1f}%" if p > 2 else "",
        startangle=90,
        wedgeprops={"edgecolor": "white", "linewidth": 1.2},
        pctdistance=0.75,
    )
    for at in autotexts:
        at.set_fontsize(8)

    legend_patches = [
        mpatches.Patch(
            color=CATEGORY_COLORS[row["category"]],
            label=f"{row['label']}  (n={row['count']:,})",
        )
        for _, row in nonzero.iterrows()
    ]
    ax.legend(
        handles=legend_patches,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.22),
        fontsize=8,
        frameon=False,
        ncol=1,
    )

    ax.set_title(
        "HAR genomic context — proportions", fontsize=12, fontweight="normal", pad=10
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Pie chart saved: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Classify HARs by genomic context (promoter, exon, intron, intergenic)."
    )
    parser.add_argument(
        "--bed",
        required=True,
        help="Path to HAR BED file in hg38 coordinates (e.g. data/hars_hg38.bed)",
    )
    parser.add_argument(
        "--gtf",
        default=None,
        help=(
            "Path to GENCODE v38 GTF (gz or plain). "
            f"Downloaded automatically to {DEFAULT_GTF_CACHE} if not provided."
        ),
    )
    parser.add_argument(
        "--output_dir",
        default=str(_RESULTS / "har_categories"),
        help="Directory for output files (default: results/har_categories/)",
    )
    parser.add_argument(
        "--promoter_window",
        type=int,
        default=2000,
        help="bp upstream/downstream of TSS to define promoter-proximal (default: 2000)",
    )
    parser.add_argument(
        "--distal_window",
        type=int,
        default=500_000,
        help="Distance threshold (bp) separating intergenic-proximal from distal (default: 500000)",
    )
    args = parser.parse_args()

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    work_dir = os.path.join(args.output_dir, "_annotation_beds")
    Path(work_dir).mkdir(exist_ok=True)

    print("=" * 65)
    print("HAR GENOMIC CONTEXT CLASSIFICATION")
    print("=" * 65)

    # 1. GTF
    gtf_path = args.gtf if args.gtf else DEFAULT_GTF_CACHE
    print("\n[1/4] Obtaining GENCODE v38 annotation...")
    download_gtf(gtf_path)

    # 2. Parse GTF to BED files
    print("\n[2/4] Building annotation BED files...")
    beds = gtf_to_beds(gtf_path, work_dir, promoter_window=args.promoter_window)

    # 3. Classify
    print("\n[3/4] Classifying HARs...")
    df = classify_hars(
        har_bed=args.bed,
        beds=beds,
        distal_window=args.distal_window,
    )

    # 4. Save outputs
    print("\n[4/4] Writing outputs...")
    tsv_path = os.path.join(args.output_dir, "har_categories.tsv")
    summary_path = os.path.join(args.output_dir, "category_summary.tsv")
    bar_path = os.path.join(args.output_dir, "har_category_barchart.pdf")
    pie_path = os.path.join(args.output_dir, "har_category_piechart.pdf")

    df.to_csv(tsv_path, sep="\t", index=False)
    print(f"  Per-HAR table saved: {tsv_path}")

    summary = build_summary(df)
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"  Summary table saved: {summary_path}")

    plot_barchart(summary, bar_path)
    plot_piechart(summary, pie_path)

    # Print summary to terminal
    print("\n" + "=" * 65)
    print("CLASSIFICATION SUMMARY")
    print("=" * 65)
    for _, row in summary.iterrows():
        print(
            f"  {row['label']:<45s}  {int(row['count']):>5,}  ({row['percent']:>5.1f}%)"
        )
    print(f"\n  Total HARs classified: {len(df):,}")
    print(f"\nOutput directory: {args.output_dir}/")
    print("=" * 65)


if __name__ == "__main__":
    main()


"""
Output files:
    har_categories/
    ├── har_categories.tsv          # per-HAR: coord, category, nearest_gene, distance
    ├── category_summary.tsv        # counts and percentages per category
    ├── har_category_barchart.pdf   # horizontal bar chart
    ├── har_category_piechart.pdf   # pie chart
    └── _annotation_beds/           # intermediate BED files (cached, safe to delete)
        ├── tss.bed
        ├── exons.bed
        └── gene_bodies.bed
"""
