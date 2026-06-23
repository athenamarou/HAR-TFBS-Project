"""
align.py

This script forms the core of the HAR (Human Accelerated Region) comparative analysis pipeline.
It takes orthologous genomic sequences from human and chimpanzee, performs pairwise sequence
alignments using MAFFT, and calculates nucleotide-level evolutionary changes.

It handles sequence pairing safely (avoiding liftOver dropouts) and generates
coordinate mapping dictionaries so downstream tools can correctly find Transcription Factor
Binding Sites (TFBS) onto the gapped alignment.
"""

import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path
import pandas as pd
import openpyxl  # Required by pandas for reading .xlsx files
from Bio import SeqIO

_ROOT = Path(__file__).resolve().parent.parent
_DATA = _ROOT / "data"
_RESULTS = _ROOT / "results"


def run_mafft(human_seq, chimp_seq, region_id, mafft_bin="mafft"):
    """
    Executes the external MAFFT aligner on a pair of sequences.

    Since MAFFT requires file inputs, this function securely creates a temporary
    FASTA file, writes the sequences to it, calls the MAFFT binary via the shell,
    parses the output, and then deletes the temporary file.
    """
    # Create a temporary file to hold the unaligned sequences.
    # delete=False is used so we can safely close it and let MAFFT read it.
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        tmp_path = tmp.name
        # Write both sequences in standard FASTA format
        tmp.write(f">human|{region_id}\n{human_seq}\n")
        tmp.write(f">chimp|{region_id}\n{chimp_seq}\n")

    try:
        # Run MAFFT.
        # --auto lets MAFFT choose the best alignment strategy based on sequence length.
        # --quiet suppresses the standard error progress output.
        result = subprocess.run(
            [mafft_bin, "--auto", "--quiet", tmp_path], capture_output=True, text=True
        )
    except FileNotFoundError:
        # Catch the error if MAFFT isn't installed or not in the system PATH
        raise RuntimeError("MAFFT binary not found.")
    finally:
        # Ensure the temporary file is deleted even if MAFFT crashes
        os.unlink(tmp_path)

    # If MAFFT failed (e.g., out of memory), return None
    if result.returncode != 0:
        return None, None

    # Parse MAFFT's stdout (which is a multi-FASTA string of the aligned sequences)
    seqs = {}
    current_id = None
    current_parts = []

    for line in result.stdout.strip().splitlines():
        if line.startswith(">"):
            # If we hit a new header, save the previous sequence chunks
            if current_id is not None:
                seqs[current_id] = "".join(current_parts).upper()
            # Extract species name (e.g., 'human' or 'chimp')
            current_id = line[1:].split("|")[0]
            current_parts = []
        else:
            # Accumulate sequence lines
            current_parts.append(line.strip())

    # Save the final sequence after the loop ends
    if current_id is not None:
        seqs[current_id] = "".join(current_parts).upper()

    return seqs.get("human", ""), seqs.get("chimp", "")


def build_pos_map(aligned_seq):
    """
    Creates a mapping from ungapped (original) sequence coordinates to
    gapped (alignment) coordinates.

    If the TFBS scanner found a motif at position 10, but
    MAFFT inserted 2 gaps before that position, the visualization needs to draw
    the box at position 12. This dictionary acts as the translation layer.
    """
    pos_map = {}
    ungapped = 0
    # Enumerate through every column in the alignment
    for aln_pos, char in enumerate(aligned_seq):
        if char != "-":
            # If it's a real nucleotide, map its original index to this alignment column
            pos_map[ungapped] = aln_pos
            ungapped += 1
    return pos_map


def conservation_track(aligned_human, aligned_chimp):
    """
    Compares two aligned sequences column-by-column to determine the evolutionary
    event at every single position. Returns a list of classification strings.
    """
    track = []
    # zip() iterates through both sequences simultaneously
    for h, c in zip(aligned_human, aligned_chimp):
        if h == "-" and c == "-":
            # Very rare artifact of alignment, usually ignored
            track.append("gap_both")
        elif h == "-":
            # Human has a gap, meaning chimp has an insertion (or human deletion)
            track.append("gap_human")
        elif c == "-":
            # Chimp has a gap, meaning human has an insertion (or chimp deletion)
            track.append("gap_chimp")
        elif h == c:
            # Nucleotides match perfectly
            track.append("conserved")
        else:
            # Nucleotides differ (e.g., A in human, G in chimp)
            track.append("substitution")
    return track


def build_har_name_lookup(xlsx_path, sheet=0):
    """
    Reads the Cui et al. supplementary Excel table to map raw genomic coordinates
    (e.g., chr1:2920237-2920259) to standardized names (e.g., HAR_1).
    """
    # header=2 skips the title rows in the specific Cui et al. supplementary table format
    df = pd.read_excel(xlsx_path, sheet_name=sheet, header=2)
    lookup = {}

    for _, row in df.iterrows():
        try:
            # Clean and format coordinates to create the dictionary key
            chrom = str(row["chr_hg38"]).strip()
            start = int(row["start_hg38"])
            end = int(row["end_hg38"])
            name = str(row["Names"]).strip()
            key = f"{chrom}:{start}-{end}"
            lookup[key] = name
        except (ValueError, KeyError):
            # Skip rows that are malformed or missing data
            continue

    print(f"Loaded {len(lookup)} HAR name mappings from {xlsx_path}")
    return lookup


def align_all_hars(
    tsv_path,
    output_dir,
    mafft_bin="mafft",
    fasta_dir=".",
    bed_dir=".",
    har_names_xlsx=None,
):
    """
    It loads data, resolves LiftOver dropouts,
    runs alignments, builds mapping tracks, and saves everything to disk.
    """
    # Ensure output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # 1. Load the HAR name lookup table if provided
    har_name_lookup = {}
    if har_names_xlsx and os.path.exists(har_names_xlsx):
        har_name_lookup = build_har_name_lookup(har_names_xlsx)

    def region_label(coord_str):
        # Returns standard name (HAR_1) if found, otherwise returns raw coordinates
        return har_name_lookup.get(coord_str, coord_str)

    # 2. Parse the local FASTA files containing the raw sequences
    print("Loading local FASTA files (Bypassing NCBI to prevent HTTP 400 errors)...")
    human_fasta_path = os.path.join(fasta_dir, "human_hars.fasta")
    chimp_fasta_path = os.path.join(fasta_dir, "chimp_hars.fasta")

    # Store sequences in memory as dictionaries mapped by their sequence ID
    human_seqs = {
        rec.id: str(rec.seq).upper() for rec in SeqIO.parse(human_fasta_path, "fasta")
    }
    chimp_seqs = {
        rec.id: str(rec.seq).upper() for rec in SeqIO.parse(chimp_fasta_path, "fasta")
    }

    # 3. Robust Sequence Pairing
    # We must explicitly pair human and chimp regions using the BED files.
    # Naive pairing by index breaks if some regions failed to map (liftOver dropouts).
    human_bed = os.path.join(bed_dir, "hars_hg38.bed")
    chimp_bed = os.path.join(bed_dir, "hars_panTro6.bed")
    unmapped_bed = os.path.join(bed_dir, "unmapped.bed")

    # Collect all human coordinates that failed liftOver mapping
    unmapped_coords = set()
    if os.path.exists(unmapped_bed):
        with open(unmapped_bed) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    parts = line.strip().split()
                    unmapped_coords.add(f"{parts[0]}:{parts[1]}-{parts[2]}")

    paired_regions = []
    if os.path.exists(human_bed) and os.path.exists(chimp_bed):
        with open(human_bed) as fh, open(chimp_bed) as fc:
            human_lines = [l.strip().split() for l in fh if l.strip()]
            chimp_lines = [l.strip().split() for l in fc if l.strip()]

            c_idx = 0
            for h in human_lines:
                h_coord = f"{h[0]}:{h[1]}-{h[2]}"
                # If this human region failed mapping, skip it to keep lists synchronized
                if h_coord in unmapped_coords:
                    continue

                # Create the tuple matching the human coordinate to its specific chimp ortholog
                if c_idx < len(chimp_lines):
                    c = chimp_lines[c_idx]
                    c_coord = f"{c[0]}:{c[1]}-{c[2]}"
                    paired_regions.append((h_coord, c_coord))
                    c_idx += 1
        print(
            f"Reconstructed {len(paired_regions)} mathematically perfect pairs directly from BED files."
        )
    else:
        print("CRITICAL ERROR: Original BED files not found. Cannot reconstruct pairs.")
        return

    results = {}

    # 4. Iterate over pairs and process alignments
    for human_region, chimp_region in paired_regions:
        har_name = region_label(human_region)
        print(f"Aligning: {har_name} ({human_region}) <-> {chimp_region}")

        # Create a filesystem-safe identifier (removes colons and dashes)
        safe_id = har_name.replace(":", "_").replace("-", "_")

        human_seq = human_seqs.get(human_region)
        chimp_seq = chimp_seqs.get(chimp_region)

        # Skip if either sequence is missing from the parsed FASTA data
        if not human_seq or not chimp_seq:
            print(f"  Skipping {human_region}: Sequence missing from FASTA files.")
            continue

        # Execute alignment
        aligned_human, aligned_chimp = run_mafft(
            human_seq, chimp_seq, safe_id, mafft_bin
        )

        if aligned_human is None:
            print(f"  Skipping {human_region}: MAFFT failed.")
            continue

        # Export 1: Write raw alignment to FASTA
        fasta_path = os.path.join(output_dir, f"{safe_id}.aligned.fasta")
        with open(fasta_path, "w") as f:
            f.write(f">human|{human_region}\n{aligned_human}\n")
            f.write(f">chimp|{chimp_region}\n{aligned_chimp}\n")

        # Export 2: Generate and save position maps (ungapped -> gapped coordinates)
        human_pos_map = build_pos_map(aligned_human)
        chimp_pos_map = build_pos_map(aligned_chimp)

        posmap_path = os.path.join(output_dir, f"{safe_id}.posmaps.json")
        with open(posmap_path, "w") as f:
            json.dump(
                {
                    "human_region": human_region,
                    "chimp_region": chimp_region,
                    "human_pos_map": {str(k): v for k, v in human_pos_map.items()},
                    "chimp_pos_map": {str(k): v for k, v in chimp_pos_map.items()},
                },
                f,
            )

        # Export 3: Generate and save conservation track
        track = conservation_track(aligned_human, aligned_chimp)
        conservation_path = os.path.join(output_dir, f"{safe_id}.conservation.json")
        with open(conservation_path, "w") as f:
            json.dump(
                {
                    "human_region": human_region,
                    "chimp_region": chimp_region,
                    "aln_length": len(track),
                    "track": track,
                },
                f,
            )

        # Accumulate metrics for the summary table
        results[safe_id] = {
            "har_name": har_name,
            "human_region": human_region,
            "chimp_region": chimp_region,
            "aln_length": len(aligned_human),
            "n_conserved": track.count("conserved"),
            "n_substitutions": track.count("substitution"),
            "n_gaps": track.count("gap_human") + track.count("gap_chimp"),
        }

    # 5. Export Summary TSV
    # Convert dictionary of results to a DataFrame (one row per HAR) and save
    summary_df = pd.DataFrame(results).T
    summary_path = os.path.join(output_dir, "alignment_summary.tsv")
    summary_df.to_csv(summary_path, sep="\t")
    print(f"\nAlignment summary saved to {summary_path}")

    return results


if __name__ == "__main__":
    # Command Line Interface (CLI) configuration
    parser = argparse.ArgumentParser(
        description="Pairwise MAFFT alignment of human/chimp HAR sequences."
    )
    # Define required and optional arguments mapping to function parameters
    parser.add_argument("--tsv", required=True, help="TFBS hits TSV from scan.py")
    parser.add_argument(
        "--output_dir", default=str(_RESULTS / "alignments"), help="Directory for alignment outputs"
    )
    parser.add_argument("--mafft_bin", default="mafft", help="Path to MAFFT binary")
    parser.add_argument("--fasta_dir", default=str(_DATA), help="Directory with FASTA files")
    parser.add_argument(
        "--bed_dir", default=str(_DATA), help="Directory with liftOver BED files"
    )
    parser.add_argument(
        "--har_names",
        default=None,
        help="Path to Cui et al. MOESM4 xlsx for HAR name lookup",
    )
    args = parser.parse_args()

    # Execute the primary function using the parsed arguments
    align_all_hars(
        tsv_path=args.tsv,
        output_dir=args.output_dir,
        mafft_bin=args.mafft_bin,
        fasta_dir=args.fasta_dir,
        bed_dir=args.bed_dir,
        har_names_xlsx=args.har_names,
    )

"""
Example run from command line:

python align.py \
  --tsv data/HAR_tfbs_hits.tsv \
  --fasta_dir data/ --bed_dir data/ \
  --har_names data/41586_2025_8622_MOESM4_ESM.xlsx \
  --output_dir data/alignments/

"""
