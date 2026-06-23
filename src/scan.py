from retrieve_motifs import parse_meme_file, calculate_background_frequencies

from Bio import motifs, SeqIO
from Bio.Seq import Seq
from Bio import Entrez, SeqIO
import argparse
import subprocess
import os
from pathlib import Path
import pandas as pd
import openpyxl  # required for xlsx reading

_ROOT = Path(__file__).resolve().parent.parent
_DATA = _ROOT / "data"
_RESULTS = _ROOT / "results"

# ----------------------------------------------------------------
#               Argument Parsing
# ----------------------------------------------------------------

parser = argparse.ArgumentParser(
    description="""

                                 Map TFBS from JASPAR on any human genomic sequence and its orthologous regions for other species
                                 Sequence can be provided in 2 ways:

                                 1. by coordinate (fetched by NCBI):
                                 python scan.py --email youremail
                                                --chrom chr14
                                                --start 33576362
                                                --end   33576566

                                2. by FASTA file (if we already have the sequence):
                                python scan.py --email youremail
                                               --human_fasta human_seq.fasta
                                               --chimp_fasta chimp_seq.fasta

                                For chimp:
                                - provide chimp_fasta directly
                                - provide --chimp_start and --chimp_end from liftOver output (manually)
                                - provide --liftover_bin
                                 """,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

# Entrez
parser.add_argument(
    "--email", required=True, type=str, help="Type your email for Entrez"
)


# Input
# Option 1: Coordinates
parser.add_argument("--chrom", default=None, help="Chromosome e.g. chr14 (hg38)")

parser.add_argument("--start", type=int, default=None, help="Start coordinate")

parser.add_argument("--end", type=int, default=None, help="End coordinate")


# Option 2: Fasta Files
parser.add_argument("--human_fasta", default=None, help="Fata file with human seq")

parser.add_argument("--chimp_fasta", default=None, help="Fata file with chimp seq")


# Chimp coordinates if liftOver has already been run
parser.add_argument(
    "--chimp_chrom", default=None, help="Chimp chromosome e.g chr14 (panTro6)"
)

parser.add_argument(
    "--chimp_start", type=int, default=None, help="Chimp start coordinate"
)

parser.add_argument("--chimp_end", type=int, default=None, help="Chimp end coordinate")

# LiftOver execution
parser.add_argument(
    "--liftover_bin", default=None, help="path to liftover binary e.g ./liftOver"
)

parser.add_argument(
    "--liftover_chain", default=None, help="path to hg38topanTro6.over.chain.gz"
)

# JASPAR motif file
parser.add_argument(
    "--meme_file",
    default=str(_DATA / "meme_file.txt"),
    help="Path to JASPAR MEME format file",
)

# Scanner settings
parser.add_argument(
    "--threshold", type=float, default=0.8, help="Threshold as fraction of max score"
)

# Output
parser.add_argument(
    "--output",
    default=str(_RESULTS / "HAR_tfbs_hits.tsv"),
    help="Output file name. Name it as prefered",
)
parser.add_argument(
    "--har_names",
    default=None,
    help="Path to Cui et al. MOESM4 xlsx for HAR name lookup (e.g. data/41586_2025_8622_MOESM4_ESM.xlsx)",
)
parser.add_argument(
    "--bed_dir",
    default=None,
    help="Directory with hars_hg38.bed, hars_panTro6.bed, unmapped.bed — used to build "
    "correct human/chimp FASTA pairs when counts differ due to liftOver dropouts.",
)
args = parser.parse_args()


# ----------------------------------------------------------------
#               HAR NAME LOOKUP (Cui et al. supplementary table)
# ----------------------------------------------------------------


def build_har_name_lookup(xlsx_path, sheet=0):
    """
    Build a dict: "chr1:2920237-2920259" -> "HAR_1"
    from the Cui et al. supplementary table (MOESM4).
    The xlsx has 2 title rows then column headers on row 3.
    """
    df = pd.read_excel(xlsx_path, sheet_name=sheet, header=2)
    lookup = {}
    for _, row in df.iterrows():
        try:
            chrom = str(row["chr_hg38"]).strip()
            start = int(row["start_hg38"])
            end = int(row["end_hg38"])
            name = str(row["Names"]).strip()
            key = f"{chrom}:{start}-{end}"
            lookup[key] = name
        except (ValueError, KeyError):
            continue
    print(f"  Loaded {len(lookup)} HAR name mappings from {xlsx_path}")
    return lookup


def build_fasta_pairs_from_bed(human_fasta, chimp_fasta, bed_dir):
    """
    Pair human and chimp SeqRecords correctly using the original BED files,
    skipping HARs that failed liftOver. Returns list of (human_record, chimp_record).
    Without this, a 3-record count mismatch shifts all downstream pairs after
    each dropout.
    """
    import os

    human_bed = os.path.join(bed_dir, "hars_hg38.bed")
    chimp_bed = os.path.join(bed_dir, "hars_panTro6.bed")
    unmapped_bed = os.path.join(bed_dir, "unmapped.bed")

    if not (os.path.exists(human_bed) and os.path.exists(chimp_bed)):
        print(
            "  Warning: BED files not found in bed_dir — falling back to index pairing."
        )
        return None

    # Collect unmapped human coords
    unmapped = set()
    if os.path.exists(unmapped_bed):
        with open(unmapped_bed) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    p = line.strip().split()
                    unmapped.add(f"{p[0]}:{p[1]}-{p[2]}")

    # Build ordered list of (human_coord, chimp_coord) skipping unmapped
    pairs = []
    with open(human_bed) as fh, open(chimp_bed) as fc:
        human_lines = [l.strip().split() for l in fh if l.strip()]
        chimp_lines = [l.strip().split() for l in fc if l.strip()]
        c_idx = 0
        for h in human_lines:
            h_coord = f"{h[0]}:{h[1]}-{h[2]}"
            if h_coord in unmapped:
                continue
            if c_idx < len(chimp_lines):
                c = chimp_lines[c_idx]
                c_coord = f"{c[0]}:{c[1]}-{c[2]}"
                pairs.append((h_coord, c_coord))
                c_idx += 1

    # Index the FASTA records by their ID
    from Bio import SeqIO

    human_idx = {rec.id: rec for rec in SeqIO.parse(human_fasta, "fasta")}
    chimp_idx = {rec.id: rec for rec in SeqIO.parse(chimp_fasta, "fasta")}

    paired = []
    skipped = 0
    for h_coord, c_coord in pairs:
        h_rec = human_idx.get(h_coord)
        c_rec = chimp_idx.get(c_coord)
        if h_rec and c_rec:
            paired.append((h_rec, c_rec))
        else:
            skipped += 1
    print(
        f"  BED-based pairing: {len(paired)} valid pairs, {skipped} skipped (missing FASTA records)."
    )
    return paired


# ----------------------------------------------------------------
#               NCBI ACCESSION MAPS
# ----------------------------------------------------------------

# retrieved from : https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt
accession_hg38 = {
    "chr1": "NC_000001.11",
    "chr2": "NC_000002.12",
    "chr3": "NC_000003.12",
    "chr4": "NC_000004.12",
    "chr5": "NC_000005.10",
    "chr6": "NC_000006.12",
    "chr7": "NC_000007.14",
    "chr8": "NC_000008.11",
    "chr9": "NC_000009.12",
    "chr10": "NC_000010.11",
    "chr11": "NC_000011.10",
    "chr12": "NC_000012.12",
    "chr13": "NC_000013.11",
    "chr14": "NC_000014.9",
    "chr15": "NC_000015.10",
    "chr16": "NC_000016.10",
    "chr17": "NC_000017.11",
    "chr18": "NC_000018.10",
    "chr19": "NC_000019.10",
    "chr20": "NC_000020.11",
    "chr21": "NC_000021.9",
    "chr22": "NC_000022.11",
}

# retrieved from: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/515/GCF_000001515.7_Pan_tro_3.0/GCF_000001515.7_Pan_tro_3.0_assembly_report.txt
accession_PanTro = {
    "chr1": "NC_006468.4",
    "chr2A": "NC_006469.4",
    "chr2B": "NC_006470.4",
    "chr3": "NC_006471.4",
    "chr4": "NC_006472.4",
    "chr5": "NC_006473.4",
    "chr6": "NC_006474.4",
    "chr7": "NC_006475.4",
    "chr8": "NC_006476.4",
    "chr9": "NC_006477.4",
    "chr10": "NC_006478.4",
    "chr11": "NC_006480.4",
    "chr12": "NC_006479.4",
    "chr13": "NC_006482.4",
    "chr14": "NC_006481.4",
    "chr15": "NC_006483.4",
    "chr16": "NC_006484.4",
    "chr17": "NC_006485.4",
    "chr18": "NC_006486.4",
    "chr19": "NC_006487.4",
    "chr20": "NC_006488.4",
    "chr21": "NC_006489.4",
    "chr22": "NC_006490.4",
}

# ----------------------------------------------------------------
#               FUNCTIONS
# ----------------------------------------------------------------


def fetch_sequence(chrom, start, end, accession_map):
    """Get any genomic region from NCBI by coordinates"""

    accession = accession_map.get(chrom)
    print(f"Fetching {chrom}:{start}-{end} [{accession}] from NCBI")

    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="fasta",
        retmode="text",
        seq_start=start,
        seq_stop=end,
    )
    return SeqIO.read(handle, "fasta").seq


def run_liftover(chrom, start, end, liftover_bin, chain_file):
    """
    Run UCSC liftOver to convert hg38 coords to panTro6.
    binary and chain file are required and should be downloaded from UCSC.
    Returns (chrom, start, end) for chimp, or None if it fails.

    ftps:
    chain_file-> https://hgdownload.gi.ucsc.edu/goldenPath/hg38/liftOver/hg38ToPanTro6.over.chain.gz
    binary_human --> https://hgdownload.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
    binary_chimp -->https://hgdownload.gi.ucsc.edu/goldenPath/panTro6/bigZips/panTro6.2bit
    """

    # BED files for liftOver
    _RESULTS.mkdir(exist_ok=True)
    bed_in = str(_RESULTS / "liftover_input.bed")
    bed_out = str(_RESULTS / "liftover_output.bed")
    bed_unmapped = str(_RESULTS / "liftover_unmapped.bed")  # not everything will be aligned between the two species, so we need to know which regions are unmapped

    with open(bed_in, "w") as f:
        f.write(f"{chrom}\t{start-1}\t{end}\tregion\n")  # BED is 0-based

    result = subprocess.run(
        [liftover_bin, bed_in, chain_file, bed_out, bed_unmapped],
        capture_output=True,
        text=True,
    )
    if not os.path.exists(bed_out) or os.path.getsize(bed_out) == 0:
        print(f"LiftOver failed - region unmapped")
        print(f" stderr:{result.stderr}")
        return None

    with open(bed_out) as f:
        line = f.readline().strip().split("\t")

    chimp_chrom = line[0]
    chimp_start = int(line[1]) + 1
    chimp_end = int(line[2])

    print(
        f"LiftOver completed. Result: {chimp_chrom}:{chimp_start} - {chimp_end} (panTro6)"
    )
    return chimp_chrom, chimp_start, chimp_end


def map_TFBS(sequence, motif_list, background, threshold_pct=0.8):
    """
    Map TFBS from JASPAR on any genomic sequence.
    Both strands are scanned.
    """
    hits = []
    for record in motif_list:
        m = motifs.Motif(alphabet="ACGT", counts=record["probs"])
        pwm = m.counts.normalize(pseudocounts=0.5)
        pssm = pwm.log_odds(background)
        threshold = pssm.max * threshold_pct

        # Skip if sequence is shorter than the motif width
        if len(sequence) < record["width"]:
            continue

        for position, score in pssm.search(sequence, threshold=threshold, both=True):
            strand = "+" if position >= 0 else "-"
            pos = position if position >= 0 else -position
            hits.append(
                {
                    "motif_id": record["matrix_id"],
                    "name": record["name"],
                    "position": pos,
                    "strand": strand,
                    "score": round(score, 4),
                    "relative_score": round(score / pssm.max, 4),
                }
            )
    return hits


# ----------------------------------------------------------------
#               MAIN
# ----------------------------------------------------------------

Entrez.email = args.email

# Load HAR name lookup if provided
har_name_lookup = {}
if args.har_names:
    print("Loading HAR name lookup...")
    har_name_lookup = build_har_name_lookup(args.har_names)


def region_label(coord_str):
    """Return HAR name if known, else the coordinate string."""
    return har_name_lookup.get(coord_str, coord_str)


# Load JASPAR motifs

print(f"Loading JASPAR motifs from {args.meme_file}")

motif_list, _ = parse_meme_file(args.meme_file)

background = calculate_background_frequencies(motif_list)
print(f"Loaded {len(motif_list)} motifs with background frequencies: {background}")


all_hits = []

# Get human seq

print("Human Sequence from hg38")

# If fasta file with the sequence is provided
if args.human_fasta:
    human_records = list(SeqIO.parse(args.human_fasta, "fasta"))
    print(
        f"Loaded human sequence from {args.human_fasta} with length {len(human_records)}"
    )
    for rec in human_records:
        hits = map_TFBS(rec.seq, motif_list, background, args.threshold)
        har_id = region_label(rec.id)  # HAR_1 if lookup available, else coords
        for h in hits:
            h["species"] = "human"
            h["region"] = har_id
            h["coords_hg38"] = rec.id
        all_hits.extend(hits)

    print(f"  Total human hits: {len([h for h in all_hits if h['species']=='human'])}")

elif args.chrom and args.start and args.end:
    # Using provided coordinates from NCBI
    human_seq = fetch_sequence(args.chrom, args.start, args.end, accession_hg38)
    coord_str = f"{args.chrom}:{args.start}-{args.end}"
    hits = map_TFBS(human_seq, motif_list, background, args.threshold)
    for h in hits:
        h["species"] = "human"
        h["region"] = region_label(coord_str)
        h["coords_hg38"] = coord_str
    all_hits.extend(hits)
    print(f"  Human hits: {len(hits)}")
else:
    raise ValueError("Provide either fasta or coordinates for human sequence")

# Get chimp seq
print("\nChimp sequences:")
chimp_seq = None

if args.chimp_fasta:
    chimp_records = list(SeqIO.parse(args.chimp_fasta, "fasta"))
    print(f"  Loaded {len(chimp_records)} sequences from {args.chimp_fasta}")

    # Use BED-based pairing if bed_dir provided (handles liftOver dropouts correctly)
    paired_records = None
    if args.bed_dir and args.human_fasta:
        print("  Using BED-based pairing to handle liftOver dropouts...")
        paired_records = build_fasta_pairs_from_bed(
            args.human_fasta, args.chimp_fasta, args.bed_dir
        )

    if paired_records is not None:
        # BED-based: each pair is (human_rec, chimp_rec) correctly matched
        for human_rec, chimp_rec in paired_records:
            hits = map_TFBS(chimp_rec.seq, motif_list, background, args.threshold)
            human_region = region_label(human_rec.id)
            for h in hits:
                h["species"] = "chimp"
                h["region"] = human_region
                h["coords_hg38"] = human_rec.id
                h["chimp_region"] = chimp_rec.id
            all_hits.extend(hits)
    else:
        # Fallback: index-based pairing (only safe if counts match)
        if args.human_fasta:
            human_records_for_pairing = list(SeqIO.parse(args.human_fasta, "fasta"))
        else:
            human_records_for_pairing = []
        if len(human_records_for_pairing) != len(chimp_records):
            print(
                f"  WARNING: human ({len(human_records_for_pairing)}) and chimp "
                f"({len(chimp_records)}) FASTA counts differ. "
                f"Provide --bed_dir for correct pairing."
            )
        for i, rec in enumerate(chimp_records):
            hits = map_TFBS(rec.seq, motif_list, background, args.threshold)
            if i < len(human_records_for_pairing):
                human_coord = human_records_for_pairing[i].id
                human_region = region_label(human_coord)
                human_coords = human_coord
            else:
                human_region = region_label(rec.id)
                human_coords = rec.id
            for h in hits:
                h["species"] = "chimp"
                h["region"] = human_region
                h["coords_hg38"] = human_coords
                h["chimp_region"] = rec.id
            all_hits.extend(hits)

    print(f"  Total chimp hits: {len([h for h in all_hits if h['species']=='chimp'])}")
    chimp_seq = True  # flag for comparison block below

elif args.chimp_chrom and args.chimp_start and args.chimp_end:
    chimp_seq = fetch_sequence(
        args.chimp_chrom, args.chimp_start, args.chimp_end, accession_PanTro
    )
    hits = map_TFBS(chimp_seq, motif_list, background, args.threshold)
    human_coord = (
        f"{args.chrom}:{args.start}-{args.end}"
        if args.chrom
        else f"{args.chimp_chrom}:{args.chimp_start}-{args.chimp_end}"
    )
    for h in hits:
        h["species"] = "chimp"
        h["region"] = region_label(human_coord)
        h["coords_hg38"] = human_coord
        h["chimp_region"] = f"{args.chimp_chrom}:{args.chimp_start}-{args.chimp_end}"
    all_hits.extend(hits)
    print(f"  Chimp hits: {len(hits)}")

elif args.liftover_bin and args.liftover_chain and args.chrom and args.start:
    print("  Running liftOver to get panTro6 coordinates...")
    result = run_liftover(
        args.chrom, args.start, args.end, args.liftover_bin, args.liftover_chain
    )
    if result:
        chimp_chrom, chimp_start, chimp_end = result
        chimp_seq = fetch_sequence(
            chimp_chrom, chimp_start, chimp_end, accession_PanTro
        )
        hits = map_TFBS(chimp_seq, motif_list, background, args.threshold)
        human_coord = f"{args.chrom}:{args.start}-{args.end}"
        for h in hits:
            h["species"] = "chimp"
            h["region"] = region_label(human_coord)
            h["coords_hg38"] = human_coord
            h["chimp_region"] = f"{chimp_chrom}:{chimp_start}-{chimp_end}"
        all_hits.extend(hits)
        print(f"  Chimp hits: {len(hits)}")
else:
    print("  No chimp sequence provided — scanning human only.")

# --- Export TSV ---
df = pd.DataFrame(all_hits)
df.sort_values(
    ["species", "region", "score"], ascending=[True, True, False], inplace=True
)
df.to_csv(args.output, sep="\t", index=False)
print(f"\nSaved {len(df)} total hits to {args.output}")

# --- Species comparison ---
if chimp_seq is not None:
    human_tfs = set(df[df["species"] == "human"]["name"])
    chimp_tfs = set(df[df["species"] == "chimp"]["name"])
    print(f"\nHuman-only TFs: {sorted(human_tfs - chimp_tfs)}")
    print(f"Chimp-only TFs: {sorted(chimp_tfs - human_tfs)}")
    print(f"Shared TFs:     {len(human_tfs & chimp_tfs)}")


"""

EXAMPLE RUNS FRO COMMAND LINE:
# Human + chimp via coordinates + automatic liftOver
python scan.py \
  --email you@email.com \
  --chrom chr14 --start 33576362 --end 33576566 \
  --liftover_bin ./liftOver \
  --liftover_chain hg38ToPanTro6.over.chain.gz \
  --output HAR202_hits.tsv

# Human + chimp with manual liftOver coords already known
python scan.py \
  --email you@email.com \
  --chrom chr14 --start 33576362 --end 33576566 \
  --chimp_chrom chr14 --chimp_start 33116156 --chimp_end 33116360 \
  --output HAR202_hits.tsv

# Human only from FASTA
python scan.py \
  --email you@email.com \
  --human_fasta my_region.fasta \
  --output hits.tsv
"""
