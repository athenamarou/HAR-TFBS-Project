# HAR TFBS Pipeline

Genome-wide analysis of **Transcription Factor Binding Site (TFBS) rewiring** in Human Accelerated Regions (HARs), comparing regulatory landscapes between human (hg38) and chimpanzee (panTro6).


---

## Purpose

Human Accelerated Regions are genomic sequences conserved across mammals that show a burst of substitutions specifically in the human lineage. This pipeline asks: **do those substitutions disrupt or create transcription factor binding sites?**

For each of 3,254 HARs it:
1. Aligns the human and chimp sequences (MAFFT)
2. Scans both for TFBS (JASPAR 2024 motifs)
3. Identifies gained, lost, and shared binding events
4. Annotates each HAR by genomic context (promoter / exon / intron / intergenic)
5. Produces per-HAR alignment figures and genome-wide rewiring statistics

**Key results:**
- 3,235 of 3,254 HARs show TFBS rewiring (288,029 positions)
- 891 unique TFs involved; homeodomain factors dominate (HOXC8, HOXB2, BSX, HOXA1)
- 7 HARs completely lost all chimp binding sites in human — 4 intronic, 2 intergenic, 1 promoter-proximal
- 91.7% average sequence conservation between human and chimp across HARs

---

## Repository structure

```
HAR-TFBS-Pipeline/
├── src/                          # Pipeline scripts
│   ├── align.py                  # MAFFT pairwise alignment (human ↔ chimp)
│   ├── scan.py                   # JASPAR TFBS scanning
│   ├── retrieve_motifs.py        # MEME/JASPAR motif parser
│   ├── HARs_categories.py        # Genomic annotation via GENCODE
│   ├── generate_tf_comparison_matrix.py
│   ├── summarize_rewiring.py     # Gained / lost / shared classification
│   ├── plot_TF_rewiring.py       # Rewiring heatmaps and bar charts
│   ├── visualize_alignment.py    # Per-HAR alignment + TFBS overlay figures
│   └── generate_alignment_report.py
│
├── data/                         # Input data (reference genomes excluded — see data/README.md)
│   ├── hars_hg38.bed             # HAR coordinates (hg38)
│   ├── hars_panTro6.bed          # HAR coordinates (panTro6, after liftOver)
│   ├── human_hars.fasta          # Human sequences
│   ├── chimp_hars.fasta          # Chimp sequences
│   ├── meme_file.txt             # JASPAR 2024 motifs
│   └── supplementary/            # Cui et al. 2025 supplementary tables
│
└── results/                      # Key outputs (bulk alignments/figures excluded)
    ├── analysis_report/          # Alignment and TFBS statistics
    ├── har_categories/           # Genomic annotation counts and charts
    ├── rewiring_report/          # Rewiring summaries and per-HAR tables
    └── example_figures/          # Representative output figures
```

---

## Setup

```bash
git clone https://github.com/YOUR_USERNAME/HAR-TFBS-Pipeline.git
cd HAR-TFBS-Pipeline
pip install -r requirements.txt
```

External tools needed: **MAFFT** (`brew install mafft`) and **UCSC liftOver** (download from UCSC Genome Browser). Download reference genomes as described in [data/README.md](data/README.md).

---

## Running the pipeline

```bash
# 1. Lift over HAR coordinates from hg38 to panTro6
./liftOver data/hars_hg38.bed data/reference/hg38ToPanTro6.over.chain.gz \
           data/hars_panTro6.bed data/unmapped.bed

# 2. Pairwise alignment for each HAR
python src/align.py \
  --human_fasta data/human_hars.fasta \
  --chimp_fasta data/chimp_hars.fasta \
  --output_dir results/alignments/

# 3. TFBS scanning with JASPAR motifs
python src/scan.py \
  --human_fasta data/human_hars.fasta \
  --chimp_fasta data/chimp_hars.fasta \
  --motifs data/meme_file.txt \
  --output_dir results/tf_matrices/

# 4. Build comparison matrices and summarise rewiring
python src/generate_tf_comparison_matrix.py
python src/summarize_rewiring.py --output_dir results/rewiring_report/

# 5. Annotate HARs by genomic context
python src/HARs_categories.py \
  --bed data/hars_hg38.bed \
  --annotation data/reference/gencode.v38.annotation.gtf.gz \
  --output results/har_categories/

# 6. Generate figures
python src/plot_TF_rewiring.py \
  --rewiring results/rewiring_report/rewiring_events.tsv \
  --output_dir results/figures/

python src/generate_alignment_report.py \
  --alignments_dir results/alignments/ \
  --output_dir results/figures/
```

---

## Example outputs

**TF rewiring overview** (`results/example_figures/tf_rewiring_combined.pdf`)

Summary of TFs with the most gained and lost binding positions across all HARs. Top gained: ZGLP1 (3,500 positions), VAX2 (1,568), HOXC8 (1,425). Top lost: ZGLP1 (3,523), VAX2 (1,810), HOXC8 (1,671).

**HAR genomic categories** (`results/example_figures/` → `har_category_barchart.pdf`)

Distribution of 3,254 HARs: 51.2% intronic, 35.0% intergenic-proximal, 9.1% promoter-proximal, 4.7% exonic.

**Per-HAR alignment figure** (`results/example_figures/HAR_1304_visualization.pdf`)

Example of a HAR that completely lost all TFBS in human — intronic in the *DCC* gene, losing 8 unique TF binding sites.

**Per-HAR alignment figure** (`results/example_figures/HAR_202_visualization.pdf`)

Example of a HAR with mixed rewiring — shows the colour-coded alignment, TFBS overlay, and conservation track.

---

## Data sources

| File | Source |
|------|--------|
| HAR coordinates | Girskis et al. 2021, *Cell Stem Cell* |
| hg38, panTro6 | UCSC Genome Browser |
| LiftOver chain | UCSC Genome Browser |
| JASPAR motifs | JASPAR 2024 vertebrate CORE |
| GENCODE annotation | GENCODE v38 |
| Supplementary tables | Cui et al. 2025, *Nature* |

---

## Citation

```
Marougka, A. (2026). HAR TFBS Pipeline: Transcription Factor Binding Site
Rewiring in Human Accelerated Regions. MSc Bioinformatics, Algorithms course.
https://github.com/athenamarou/HAR-TFBS-Pipeline
```
