"""
generate_tf_comparison_matrix.py

Builds a TF binding comparison matrix for each HAR.

This script outputs ONE file per HAR containing all TFs as rows, with columns:

    aln_pos | event | tf_name | tf_human | tf_chimp | binding_outcome

The `binding_outcome` column classifies each hit position:
    - shared        : both species have the TF bound here
    - human_only    : only human has the TF bound here
    - chimp_only    : only chimp has the TF bound here
    - absent        : neither species has the TF bound (rows are only written for hit positions)

Output options (set via --format):
    - parquet  (default) : one .parquet per HAR — fast I/O, minimal disk space
    - tsv                : one .tsv per HAR — human-readable
    - single_tsv         : one global HAR_tf_matrix.tsv with a `har` column added

Usage:
    python generate_tf_comparison_matrix.py
    python generate_tf_comparison_matrix.py --format tsv
    python generate_tf_comparison_matrix.py --format single_tsv
    python generate_tf_comparison_matrix.py --path ./data --format parquet
"""

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parent.parent
_RESULTS = _ROOT / "results"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def unique_hars(alignments_dir: str) -> set:
    """Return the set of HAR base names found in the alignments directory."""
    hars = set()
    for filename in os.listdir(alignments_dir):
        full = os.path.join(alignments_dir, filename)
        if os.path.isfile(full):
            hars.add(filename.split(".")[0])
    return hars


def load_alignment_support(posmap_path: str, cons_path: str):
    """Load and return (h_map, c_map, aln_len, track) from JSON files."""
    with open(posmap_path) as f:
        maps = json.load(f)
    with open(cons_path) as f:
        cons = json.load(f)

    h_map = {int(k): v for k, v in maps["human_pos_map"].items()}
    c_map = {int(k): v for k, v in maps["chimp_pos_map"].items()}
    return h_map, c_map, cons["aln_length"], cons["track"]


def binding_outcome(h: int, c: int) -> str:
    if h and c:
        return "shared"
    if h:
        return "human_only"
    if c:
        return "chimp_only"
    return "absent"


# ---------------------------------------------------------------------------
# Core: build the matrix for one HAR
# ---------------------------------------------------------------------------


def build_har_matrix(
    har_name: str, har_hits: pd.DataFrame, posmap_path: str, cons_path: str
) -> pd.DataFrame:
    """
    Returns a long-format DataFrame with one row per (aln_pos, tf_name) where
    at least one species has a hit — plus absent rows for positions that appear
    in the alignment track but have no binding at all.

    Columns: aln_pos, event, tf_name, tf_human, tf_chimp, binding_outcome
    """
    h_map, c_map, aln_len, track = load_alignment_support(posmap_path, cons_path)

    rows = []

    for tf_name in har_hits["name"].unique():
        tf_hits = har_hits[har_hits["name"] == tf_name]

        # Build binary vectors for this TF
        h_tf = [0] * aln_len
        c_tf = [0] * aln_len

        for _, hit in tf_hits[tf_hits["species"] == "human"].iterrows():
            pos = h_map.get(int(hit["position"]))
            if pos is not None and pos < aln_len:
                h_tf[pos] = 1

        for _, hit in tf_hits[tf_hits["species"] == "chimp"].iterrows():
            pos = c_map.get(int(hit["position"]))
            if pos is not None and pos < aln_len:
                c_tf[pos] = 1

        # Only emit rows where at least one species has a hit (keeps files small)
        for aln_pos in range(aln_len):
            if h_tf[aln_pos] or c_tf[aln_pos]:
                rows.append(
                    {
                        "aln_pos": aln_pos,
                        "event": track[aln_pos],
                        "tf_name": tf_name,
                        "tf_human": h_tf[aln_pos],
                        "tf_chimp": c_tf[aln_pos],
                        "binding_outcome": binding_outcome(
                            h_tf[aln_pos], c_tf[aln_pos]
                        ),
                    }
                )

    if not rows:
        return pd.DataFrame(
            columns=[
                "aln_pos",
                "event",
                "tf_name",
                "tf_human",
                "tf_chimp",
                "binding_outcome",
            ]
        )

    df = pd.DataFrame(rows)
    df.sort_values(["tf_name", "aln_pos"], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Build consolidated TF comparison matrices (one file per HAR)."
    )
    parser.add_argument(
        "--path", default=str(_RESULTS), help="Root results directory (default: results/)"
    )
    parser.add_argument(
        "--format",
        choices=["parquet", "tsv", "single_tsv"],
        default="parquet",
        help=(
            "Output format. "
            "'parquet' = one .parquet per HAR (recommended). "
            "'tsv' = one .tsv per HAR. "
            "'single_tsv' = one global TSV with a 'har' column."
        ),
    )
    args = parser.parse_args()

    alignments_dir = os.path.join(args.path, "alignments")
    tfbs_tsv = os.path.join(args.path, "HAR_tfbs_hits.tsv")
    output_dir = os.path.join(args.path, "tf_matrices")
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Load the full TFBS hit table once
    print(f"Loading TFBS hits from {tfbs_tsv} ...")
    hits_df = pd.read_csv(tfbs_tsv, sep="\t", low_memory=False)
    print(f"  {len(hits_df):,} total hits, {hits_df['name'].nunique()} unique TFs")

    har_set = unique_hars(alignments_dir)
    print(f"Found {len(har_set)} HARs in {alignments_dir}\n")

    all_frames = []  # used only for single_tsv mode
    total_rows = 0
    skipped = 0

    for har in sorted(har_set):
        posmap_path = os.path.join(alignments_dir, f"{har}.posmaps.json")
        cons_path = os.path.join(alignments_dir, f"{har}.conservation.json")

        if not (os.path.exists(posmap_path) and os.path.exists(cons_path)):
            print(f"  [SKIP] {har}: missing .posmaps.json or .conservation.json")
            skipped += 1
            continue

        # Filter hits to this HAR only
        har_hits = hits_df[hits_df["region"] == har]
        if har_hits.empty:
            print(f"  [SKIP] {har}: no TFBS hits in TSV")
            skipped += 1
            continue

        matrix = build_har_matrix(har, har_hits, posmap_path, cons_path)
        total_rows += len(matrix)

        if args.format == "parquet":
            out_path = os.path.join(output_dir, f"{har}.parquet")
            matrix.to_parquet(out_path, index=False)
        elif args.format == "tsv":
            out_path = os.path.join(output_dir, f"{har}_tf_matrix.tsv")
            matrix.to_csv(out_path, sep="\t", index=False)
        elif args.format == "single_tsv":
            matrix.insert(0, "har", har)
            all_frames.append(matrix)

        print(
            f"  {har}: {len(matrix):>6,} rows  "
            f"({har_hits['name'].nunique()} TFs, "
            f"{(matrix['binding_outcome']=='shared').sum()} shared, "
            f"{(matrix['binding_outcome']=='human_only').sum()} human-only, "
            f"{(matrix['binding_outcome']=='chimp_only').sum()} chimp-only)"
        )

    # Write single global TSV if requested
    if args.format == "single_tsv" and all_frames:
        global_path = os.path.join(output_dir, "HAR_tf_matrix.tsv")
        global_df = pd.concat(all_frames, ignore_index=True)
        global_df.to_csv(global_path, sep="\t", index=False)
        print(f"\nGlobal matrix saved to {global_path}  ({len(global_df):,} rows)")

    print(f"\n{'='*60}")
    print(f"Done.  {len(har_set) - skipped} HARs processed, {skipped} skipped.")
    print(f"Total rows written: {total_rows:,}")
    print(f"Output directory:   {output_dir}")

    if args.format == "parquet":
        print("\nTo read a HAR matrix in Python:")
        print("  df = pd.read_parquet('data/tf_matrices/HAR_1.parquet')")
    print("=" * 60)


if __name__ == "__main__":
    main()
