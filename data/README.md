# Data

Small input files are tracked here. Large reference genomes must be downloaded separately.

## Tracked files

| File | Description |
|------|-------------|
| `hars_hg38.bed` | HAR coordinates in hg38 |
| `hars_hg38.tsv` | HAR metadata table |
| `hars_panTro6.bed` | HAR coordinates in panTro6 (after liftOver) |
| `unmapped.bed` | HARs that could not be lifted over |
| `human_hars.fasta` | Human HAR sequences |
| `chimp_hars.fasta` | Chimpanzee HAR sequences |
| `meme_file.txt` | JASPAR 2024 vertebrate CORE motifs (MEME format) |
| `supplementary/` | Cui et al. 2025 (Nature) supplementary tables |

## Reference genomes (not tracked — download before running)

```bash
# Human genome (hg38)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz && samtools faidx hg38.fa
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit

# Chimpanzee genome (panTro6)
wget https://hgdownload.soe.ucsc.edu/goldenPath/panTro6/bigZips/panTro6.fa.gz
gunzip panTro6.fa.gz && samtools faidx panTro6.fa
wget https://hgdownload.soe.ucsc.edu/goldenPath/panTro6/bigZips/panTro6.2bit

# LiftOver chain
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToPanTro6.over.chain.gz

# GENCODE annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
```

Place all downloaded files in `data/reference/`.
