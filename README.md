
# Sequence Diversity Analysis - README

This folder contains a single R script that performs comparative sequence diversity
analyses between **Tetrasphaera sp. Soli756 (GCA_001428065)** and **Escherichia coli K-12 MG1655 (GCF_000005845.2)**.

## Files in this package
- `sequence_diversity_analysis.R` : The main R script that downloads CDS & protein FASTA files from NCBI,
  computes summary statistics, makes plots, and writes CSV outputs. Outputs are saved into the `seq_diversity_outputs/` directory created by the script.
- `seq_diversity_outputs/` : (created when the script runs) contains CSV tables and PNG plots:
  - `cds_summary_table.csv` : number of CDS and total CDS bp per organism
  - `cds_length_summary.csv`   : mean and median CDS length per organism
  - `cds_length_boxplot.png`   : boxplot of CDS lengths
  - `nucleotide_counts.csv` and `nucleotide_counts.png` : nucleotide counts (A,C,G,T) across CDS
  - `amino_acid_counts.csv` and `amino_acid_counts.png` : amino acid counts across proteins
  - `codon_usage_table.csv` and `codon_counts.png` : codon counts and comparison
  - `kmer_counts_k3.csv` ... `kmer_counts_k5.csv` : full k-mer counts for k=3,4,5
  - `top10_over_k3.png`, `top10_under_k3.png`, etc. : plots for top over/under-represented k-mers
  - `sequence_diversity_results.RData` : saved R workspace with objects for inspection

## How the script works (short)
1. Uses the `rentrez` package to query NCBI Assembly for the assembly FTP path.
2. Downloads `_cds_from_genomic.fna.gz` and `_protein.faa.gz` files from the assembly FTP directory.
3. Reads CDS (DNA) and proteins (AA) using `Biostrings` into R objects.
4. Computes counts, lengths, frequencies, codon usage and k-mer comparisons and saves results and plots to disk.

## Requirements
- R (>=4.0)
- Internet connection to download files from NCBI
- R packages: `Biostrings`, `rentrez`, `dplyr`, `ggplot2`, `tidyr`, `stringr`, `seqinr` (the script will attempt to install missing packages automatically)

## Running the script
1. Open RStudio in this folder.
2. Run the script (Source or `Rscript sequence_diversity_analysis.R`).
3. Check the `seq_diversity_outputs/` directory for CSVs and plots.

## Notes & troubleshooting
- If NCBI changes assembly FTP layout, the script may fail to find files; check the output FTP path printed in the console.
- If assemblies lack `_cds_from_genomic.fna.gz` or `_protein.faa.gz` files, use alternate files or download manually and place into `seq_diversity_outputs/` with the expected names.

