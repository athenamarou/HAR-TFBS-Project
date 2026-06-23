"""

Result : "Rewiring Summary" that answers:
  - Which TFs gained binding in humans (human_only)?
  - Which TFs lost binding in humans (chimp_only)?
  - Are those binding changes co-occurring with nucleotide substitutions or gaps?
  - Which HARs are completely lost in human (only chimp_only binding, no shared and no
    human_only)? And what are their genomic categories?

Output files:
  - rewiring_summary.tsv       : one row per (HAR, TF, binding_outcome) with event counts
  - rewiring_events.tsv        : one row per individual rewiring position (detailed)
  - rewiring_report.txt        : human-readable summary report
  - completely_lost_hars.tsv   : per-HAR table of HARs lost entirely in human
  - completely_lost_hars_categories.pdf : bar chart of category distribution


"""

import os
import json
import pandas as pd
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ROOT = Path(__file__).resolve().parent.parent
_DATA = _ROOT / "data"
_RESULTS = _ROOT / "results"

MATRICES_DIR = str(_RESULTS / "tf_matrices")
TFBS_TSV = str(_RESULTS / "HAR_tfbs_hits.tsv")
OUTPUT_DIR = str(_RESULTS / "rewiring_report")

HAR_CATEGORIES_TSV = str(_RESULTS / "har_categories" / "har_categories.tsv")
HAR_NAMES_XLSX = str(_DATA / "41586_2025_8622_MOESM4_ESM.xlsx")

# Palette
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


def load_matrix(har, har_dir):
    """Load the consolidated parquet for a given HAR."""
    parquet_path = os.path.join(har_dir, f"{har}.parquet")
    if not os.path.exists(parquet_path):
        return None
    return pd.read_parquet(parquet_path)


def classify_rewiring(df):
    """
    From a HAR's consolidated matrix, extract all rewiring events.
    A rewiring event is any row where binding_outcome is human_only or chimp_only.
    Returns a DataFrame with the rewiring rows enriched with a 'rewiring_type' column.
    """
    rewiring = df[df["binding_outcome"].isin(["human_only", "chimp_only"])].copy()
    # human_only = TF binding gained in human (or lost in chimp)
    # chimp_only = TF binding lost in human (or gained in chimp)
    rewiring["rewiring_type"] = rewiring["binding_outcome"].map(
        {"human_only": "gained_in_human", "chimp_only": "lost_in_human"}
    )
    return rewiring


def build_rewiring_summary(rewiring_df, har_name):
    """
    Aggregate rewiring events per TF for a single HAR.
    For each TF, count how many rewiring positions coincide with substitutions vs gaps.
    """
    rows = []
    for tf_name, group in rewiring_df.groupby("tf_name"):
        for rewiring_type, subgroup in group.groupby("rewiring_type"):
            total = len(subgroup)
            on_substitution = (subgroup["event"] == "substitution").sum()
            on_gap = subgroup["event"].isin(["gap_human", "gap_chimp"]).sum()
            on_conserved = (subgroup["event"] == "conserved").sum()

            rows.append(
                {
                    "har": har_name,
                    "tf_name": tf_name,
                    "rewiring_type": rewiring_type,
                    "total_rewiring_positions": total,
                    "on_substitution": on_substitution,
                    "on_gap": on_gap,
                    "on_conserved": on_conserved,
                    # What fraction of the rewiring is driven by actual nucleotide changes?
                    "substitution_driven_pct": (
                        round(on_substitution / total * 100, 1) if total > 0 else 0
                    ),
                }
            )
    return pd.DataFrame(rows)


def write_report(summary_df, events_df, output_file):
    """Write a human-readable text report."""
    with open(output_file, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("TRANSCRIPTION FACTOR BINDING REWIRING REPORT\n")
        f.write("Human Accelerated Regions — Human vs Chimpanzee\n")
        f.write("=" * 70 + "\n\n")

        total_hars = summary_df["har"].nunique()
        total_tfs = summary_df["tf_name"].nunique()
        total_events = len(events_df)

        f.write(f"HARs with rewiring events: {total_hars}\n")
        f.write(f"Unique TFs involved:        {total_tfs}\n")
        f.write(f"Total rewiring positions:   {total_events}\n\n")

        # Split by direction
        gained = summary_df[summary_df["rewiring_type"] == "gained_in_human"]
        lost = summary_df[summary_df["rewiring_type"] == "lost_in_human"]

        f.write("=" * 70 + "\n")
        f.write("TFBS GAINED IN HUMAN (human_only binding)\n")
        f.write("=" * 70 + "\n\n")
        if not gained.empty:
            # Rank by total rewiring positions across all HARs
            gained_ranked = (
                gained.groupby("tf_name")["total_rewiring_positions"]
                .sum()
                .sort_values(ascending=False)
                .head(20)
            )
            for tf, count in gained_ranked.items():
                tf_rows = gained[gained["tf_name"] == tf]
                sub_driven = tf_rows["on_substitution"].sum()
                n_hars = tf_rows["har"].nunique()
                f.write(
                    f"  {tf:<20s}  positions: {count:4d}  "
                    f"substitution-driven: {sub_driven:4d}  "
                    f"HARs: {n_hars}\n"
                )
        else:
            f.write("  None found.\n")

        f.write("\n")
        f.write("=" * 70 + "\n")
        f.write("TFBS LOST IN HUMAN (chimp_only binding)\n")
        f.write("=" * 70 + "\n\n")
        if not lost.empty:
            lost_ranked = (
                lost.groupby("tf_name")["total_rewiring_positions"]
                .sum()
                .sort_values(ascending=False)
                .head(20)
            )
            for tf, count in lost_ranked.items():
                tf_rows = lost[lost["tf_name"] == tf]
                sub_driven = tf_rows["on_substitution"].sum()
                n_hars = tf_rows["har"].nunique()
                f.write(
                    f"  {tf:<20s}  positions: {count:4d}  "
                    f"substitution-driven: {sub_driven:4d}  "
                    f"HARs: {n_hars}\n"
                )
        else:
            f.write("  None found.\n")

        f.write("\n")
        f.write("=" * 70 + "\n")
        f.write("TOP HARs BY REWIRING ACTIVITY\n")
        f.write("=" * 70 + "\n\n")
        har_activity = (
            summary_df.groupby("har")["total_rewiring_positions"]
            .sum()
            .sort_values(ascending=False)
            .head(15)
        )
        for har, count in har_activity.items():
            har_rows = summary_df[summary_df["har"] == har]
            n_tfs = har_rows["tf_name"].nunique()
            f.write(
                f"  {har:<25s}  rewiring positions: {count:4d}  TFs affected: {n_tfs}\n"
            )

        f.write("\n")
        f.write("=" * 70 + "\n")
        f.write("SUBSTITUTION-DRIVEN REWIRING (most likely causal)\n")
        f.write(
            "Positions where binding change coincides with a nucleotide substitution\n"
        )
        f.write("=" * 70 + "\n\n")
        sub_driven = (
            summary_df[summary_df["on_substitution"] > 0]
            .sort_values("on_substitution", ascending=False)
            .head(20)
        )
        if not sub_driven.empty:
            for _, row in sub_driven.iterrows():
                f.write(
                    f"  {row['har']:<25s}  {row['tf_name']:<20s}  "
                    f"{row['rewiring_type']:<20s}  "
                    f"substitution positions: {row['on_substitution']}\n"
                )
        else:
            f.write("  None found.\n")

    print(f"Report written to {output_file}")


# ---------------------------------------------------------------------------
# Completely-lost-in-human detection
# ---------------------------------------------------------------------------


def build_name_to_coord(xlsx_path):
    """
    Map HAR_N -> 'chrN:start-end' from the Cui et al. MOESM4 xlsx.
    Returns {} if the file is missing or unreadable.
    """
    if not os.path.exists(xlsx_path):
        return {}
    try:
        df = pd.read_excel(xlsx_path, sheet_name=0, header=2)
    except Exception as e:
        print(f"  [WARN] Could not read {xlsx_path}: {e}")
        return {}

    mapping = {}
    for _, row in df.iterrows():
        try:
            name = str(row["Names"]).strip()
            chrom = str(row["chr_hg38"]).strip()
            start = int(row["start_hg38"])
            end = int(row["end_hg38"])
            mapping[name] = f"{chrom}:{start}-{end}"
        except (ValueError, KeyError):
            continue
    return mapping


def enrich_lost_with_categories(lost_df, categories_tsv, name_to_coord):

    lost_df = lost_df.copy()
    lost_df["har_coord"] = lost_df["har"].map(name_to_coord).fillna("")
    lost_df["category"] = ""
    lost_df["nearest_gene"] = ""
    lost_df["distance_to_gene"] = -1  # sentinel for "unmapped"

    if not os.path.exists(categories_tsv) or not name_to_coord:
        return lost_df

    try:
        cat_df = pd.read_csv(categories_tsv, sep="\t")
    except Exception as e:
        print(f"  [WARN] Could not read {categories_tsv}: {e}")
        return lost_df

    cat_lookup = cat_df.set_index("har_coord").to_dict("index")

    for idx, row in lost_df.iterrows():
        coord = row["har_coord"]
        if coord and coord in cat_lookup:
            rec = cat_lookup[coord]
            lost_df.at[idx, "category"] = rec.get("category", "")
            lost_df.at[idx, "nearest_gene"] = rec.get("nearest_gene", "")
            lost_df.at[idx, "distance_to_gene"] = rec.get("distance_to_gene", "")

    return lost_df


def append_lost_hars_report(output_file, lost_df):
    """
    'COMPLETELY LOST IN HUMAN' section
    """
    with open(output_file, "a") as f:
        f.write("\n")
        f.write("=" * 70 + "\n")
        f.write("HARs COMPLETELY LOST IN HUMAN\n")
        f.write("(only chimp_only binding — no shared, no human_only hits)\n")
        f.write("=" * 70 + "\n\n")

        if lost_df.empty:
            f.write("  None found. Every HAR retains at least one shared or\n")
            f.write("  human-specific TFBS hit.\n")
            return

        f.write(f"  Total HARs completely lost in human: {len(lost_df)}\n")
        f.write(
            f"  Total chimp-only TFBS positions across these HARs: "
            f"{int(lost_df['n_chimp_only_positions'].sum())}\n\n"
        )

        # Category breakdown
        has_categories = (lost_df["category"] != "").any()
        if has_categories:
            f.write("  Genomic category breakdown:\n")
            cat_counts = lost_df["category"].value_counts()
            total = len(lost_df)
            for cat, n in cat_counts.items():
                if not cat:
                    cat_label = "(unmapped)"
                else:
                    cat_label = cat
                pct = round(n / total * 100, 1)
                f.write(f"    {cat_label:<25s}  {n:>4d}  ({pct:>5.1f}%)\n")
            f.write("\n")

        f.write("  Per-HAR details (sorted by chimp-only positions, descending):\n\n")
        display_cols = [
            "har",
            "n_chimp_only_positions",
            "n_unique_tfs_lost",
            "category",
            "nearest_gene",
        ]
        f.write(
            "  " + lost_df[display_cols].to_string(index=False).replace("\n", "\n  ")
        )
        f.write("\n")


def plot_lost_hars_category_chart(lost_df, output_path):
    """
    Horizontal bar chart of genomic-category distribution for completely-lost
    HARs. Skipped silently if no category enrichment is available.
    """
    if lost_df.empty or not (lost_df["category"] != "").any():
        return

    cat_counts = lost_df["category"].value_counts()
    total = len(lost_df)

    cats = [c for c in CATEGORY_COLORS if c in cat_counts.index]
    counts = [int(cat_counts[c]) for c in cats]
    labels = [CATEGORY_LABELS[c] for c in cats]
    colors = [CATEGORY_COLORS[c] for c in cats]

    fig, ax = plt.subplots(figsize=(9, 4))
    bars = ax.barh(
        labels, counts, color=colors, edgecolor="white", linewidth=0.5, height=0.6
    )
    for bar, n in zip(bars, counts):
        pct = round(n / total * 100, 1)
        ax.text(
            bar.get_width() + max(counts) * 0.01,
            bar.get_y() + bar.get_height() / 2,
            f"{n:,}  ({pct}%)",
            va="center",
            ha="left",
            fontsize=9,
            color="#444441",
        )

    ax.set_xlabel("Number of HARs completely lost in human", fontsize=10)
    ax.set_title(
        "Genomic context of HARs with total loss of TFBS in human",
        fontsize=12,
        fontweight="normal",
        pad=10,
    )
    ax.invert_yaxis()
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, max(counts) * 1.25)
    ax.grid(axis="x", linestyle="--", alpha=0.4, linewidth=0.5)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Lost-HARs category chart written to {output_path}")


def main():
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    # Load the original TFBS hits to enrich with relative_score
    print(f"Loading TFBS hits from {TFBS_TSV} ...")
    hits_df = pd.read_csv(TFBS_TSV, sep="\t", low_memory=False)

    # Build a lookup: (region, tf_name, species, position) -> relative_score
    # so we can filter for high-confidence hits (score > 0.8) as Cui et al. did
    hits_df["_key"] = (
        hits_df["region"].astype(str)
        + "|"
        + hits_df["name"].astype(str)
        + "|"
        + hits_df["species"].astype(str)
        + "|"
        + hits_df["position"].astype(str)
    )
    score_lookup = hits_df.set_index("_key")["relative_score"].to_dict()

    # Discover all HAR matrix files
    har_dirs = [
        d
        for d in os.listdir(MATRICES_DIR)
        if os.path.isdir(os.path.join(MATRICES_DIR, d))
    ]

    # Also support flat structure (one TSV per HAR directly in MATRICES_DIR)
    flat_files = list(Path(MATRICES_DIR).glob("*.parquet"))
    flat_hars = [f.stem for f in flat_files]

    all_summary_rows = []
    all_event_rows = []
    lost_records = []  # HARs where binding is exclusively chimp_only

    print(f"Processing {len(flat_hars)} HAR matrix files...\n")

    for har in sorted(flat_hars):
        matrix = load_matrix(har, MATRICES_DIR)
        if matrix is None or matrix.empty:
            continue

        # Detect "completely lost in human": only chimp_only rows in this HAR
        outcomes = set(matrix["binding_outcome"].unique()) - {"absent"}
        if outcomes == {"chimp_only"}:
            top_tfs = matrix["tf_name"].value_counts().head(5).index.tolist()
            lost_records.append(
                {
                    "har": har,
                    "n_chimp_only_positions": int(len(matrix)),
                    "n_unique_tfs_lost": int(matrix["tf_name"].nunique()),
                    "n_positions_on_substitution": int(
                        (matrix["event"] == "substitution").sum()
                    ),
                    "n_positions_on_gap": int(
                        matrix["event"].isin(["gap_human", "gap_chimp"]).sum()
                    ),
                    "top5_lost_tfs": ", ".join(top_tfs),
                }
            )

        rewiring = classify_rewiring(matrix)
        if rewiring.empty:
            continue

        # Enrich each rewiring row with relative_score from the original TSV
        def get_score(row):
            species = "human" if row["tf_human"] == 1 else "chimp"
            key = f"{har}|{row['tf_name']}|{species}|{row['aln_pos']}"
            return score_lookup.get(key, None)

        rewiring = rewiring.copy()
        rewiring["relative_score"] = rewiring.apply(get_score, axis=1)
        rewiring.insert(0, "har", har)

        # Filter for high-confidence hits (relative_score > 0.8)
        # Keep all but flag low-confidence ones
        rewiring["high_confidence"] = rewiring["relative_score"].apply(
            lambda s: True if s is None or s >= 0.8 else False
        )

        all_event_rows.append(rewiring)

        # Build per-HAR summary
        summary = build_rewiring_summary(rewiring, har)
        all_summary_rows.append(summary)

    if not all_summary_rows:
        print(
            "No rewiring events found. Check that tf_matrices/ contains *_tf_matrix.tsv files."
        )
        return

    summary_df = pd.concat(all_summary_rows, ignore_index=True)
    events_df = pd.concat(all_event_rows, ignore_index=True)

    # Save outputs
    summary_path = os.path.join(OUTPUT_DIR, "rewiring_summary.tsv")
    events_path = os.path.join(OUTPUT_DIR, "rewiring_events.tsv")
    report_path = os.path.join(OUTPUT_DIR, "rewiring_report.txt")

    summary_df.to_csv(summary_path, sep="\t", index=False)
    events_df.to_csv(events_path, sep="\t", index=False)
    write_report(summary_df, events_df, report_path)

    # Completely-lost-in-human: enrich with categories, write TSV, append to report
    lost_df = pd.DataFrame(lost_records)
    if not lost_df.empty:
        lost_df = lost_df.sort_values(
            "n_chimp_only_positions", ascending=False
        ).reset_index(drop=True)
        name_to_coord = build_name_to_coord(HAR_NAMES_XLSX)
        lost_df = enrich_lost_with_categories(
            lost_df, HAR_CATEGORIES_TSV, name_to_coord
        )
        lost_path = os.path.join(OUTPUT_DIR, "completely_lost_hars.tsv")
        lost_df.to_csv(lost_path, sep="\t", index=False)
        print(f"Completely-lost HARs written to {lost_path}")
        chart_path = os.path.join(OUTPUT_DIR, "completely_lost_hars_categories.pdf")
        plot_lost_hars_category_chart(lost_df, chart_path)
    append_lost_hars_report(report_path, lost_df)

    print(f"\n{'='*60}")
    print("REWIRING SUMMARY COMPLETE")
    print(f"{'='*60}")
    print(f"HARs processed:           {events_df['har'].nunique()}")
    print(f"Total rewiring events:    {len(events_df)}")
    print(f"Unique TFs:               {events_df['tf_name'].nunique()}")
    print(f"HARs completely lost in human: {len(lost_df)}")
    print(f"\nOutput files:")
    print(f"  {summary_path}")
    print(f"  {events_path}")
    print(f"  {report_path}")
    if not lost_df.empty:
        print(f"  {os.path.join(OUTPUT_DIR, 'completely_lost_hars.tsv')}")
        print(f"  {os.path.join(OUTPUT_DIR, 'completely_lost_hars_categories.pdf')}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
