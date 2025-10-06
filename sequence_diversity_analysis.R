
# sequence_diversity_analysis.R
# Comparative sequence diversity analysis: Tetrasphaera sp. Soli756 (GCA_001428065) vs Escherichia coli K-12 MG1655 (GCF_000005845.2)
# 
# This script will:
# - Download CDS and protein FASTA files for each assembly (using rentrez to find FTP paths)
# - Count CDS, compute lengths and total coding DNA
# - Produce CDS length boxplot and report mean/median
# - Compute nucleotide and amino-acid frequencies and plot them
# - Compute codon usage tables and plot codon usage bias
# - Identify top 10 over- and under-represented protein k-mers (k = 3..5) in Tetrasphaera vs E. coli
#
# NOTE: This script requires an internet connection and the following R packages:
#   Biostrings, rentrez, dplyr, ggplot2, tidyr, stringr, seqinr
#
# Run: in R/RStudio, set working directory to the folder containing this script and run:
#   source("sequence_diversity_analysis.R")
#
options(stringsAsFactors = FALSE)

# --------------------
# 0. Load libraries (install if necessary)
# --------------------
pkgs <- c("Biostrings", "rentrez", "dplyr", "ggplot2", "tidyr", "stringr", "seqinr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

# --------------------
# Helper: get FTP path for an assembly accession using rentrez summary
# --------------------
get_assembly_ftp <- function(assembly_accession) {
  message("Searching NCBI Assembly for ", assembly_accession)
  s <- rentrez::entrez_search(db = "assembly", term = assembly_accession)
  if (length(s$ids) == 0) stop("Assembly not found: ", assembly_accession)
  sumry <- rentrez::entrez_summary(db = "assembly", id = s$ids[1])
  # The summary object may include FtpPath_RefSeq or FtpPath_GenBank
  if (!is.null(sumry$`FtpPath_RefSeq`)) {
    return(sumry$`FtpPath_RefSeq`)
  } else if (!is.null(sumry$`FtpPath_GenBank`)) {
    return(sumry$`FtpPath_GenBank`)
  } else if (!is.null(sumry$`FtpPath`)) {
    return(sumry$`FtpPath`)
  } else {
    # Try to inspect the whole summary
    warning("No FtpPath field found in assembly summary; returning NULL")
    return(NULL)
  }
}

# --------------------
# Helper: download CDS and protein files from assembly FTP
# --------------------
download_assembly_cds_protein <- function(ftp_path, prefix = NULL, destdir = ".") {
  # ftp_path should be like: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/...
  if (is.null(ftp_path) || ftp_path == "") stop("Invalid ftp_path")
  if (is.null(prefix)) {
    prefix <- basename(ftp_path)  # e.g. GCA_001428065.1_ASM...
  }
  # filenames commonly available
  cds_file <- paste0(prefix, "_cds_from_genomic.fna.gz")
  prot_file <- paste0(prefix, "_protein.faa.gz")
  urls <- c(
    file.path(ftp_path, cds_file),
    file.path(ftp_path, prot_file)
  )
  dests <- file.path(destdir, basename(urls))
  for (i in seq_along(urls)) {
    url <- urls[i]
    dest <- dests[i]
    if (!file.exists(dest)) {
      message("Downloading ", url)
      tryCatch({
        utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
      }, error = function(e) {
        warning("Failed to download ", url, " (", e$message, ")")
      })
    } else {
      message("File already exists: ", dest)
    }
  }
  return(dests)
}

# --------------------
# Data sources: assembly accessions
# --------------------
tetra_acc <- "GCA_001428065"     # Tetrasphaera sp. Soli756 (user-specified)
ecoli_acc <- "GCF_000005845.2"   # E. coli K-12 MG1655 (RefSeq)

outdir <- "seq_diversity_outputs"
dir.create(outdir, showWarnings = FALSE)

# --------------------
# 1) Download (or locate) CDS and protein FASTA files for each organism
# --------------------
get_and_read <- function(acc, label) {
  ftp <- get_assembly_ftp(acc)
  if (is.null(ftp)) stop("FTP path not found for ", acc)
  prefix <- basename(ftp)
  files <- download_assembly_cds_protein(ftp, prefix = prefix, destdir = outdir)
  # read files (if present)
  cds_path <- files[1]
  prot_path <- files[2]
  cds <- NULL; prot <- NULL
  if (file.exists(cds_path)) {
    message("Reading CDS fasta: ", cds_path)
    cds <- Biostrings::readDNAStringSet(cds_path, format = "fasta")
  } else {
    warning("CDS fasta not found for ", label)
  }
  if (file.exists(prot_path)) {
    message("Reading protein fasta: ", prot_path)
    prot <- Biostrings::readAAStringSet(prot_path, format = "fasta")
  } else {
    warning("Protein fasta not found for ", label)
  }
  return(list(cds = cds, prot = prot, ftp = ftp, files = files))
}

message("Accessing Tetrasphaera assembly...")
tetra <- get_and_read(tetra_acc, "Tetrasphaera")
message("Accessing E. coli assembly...")
ecoli <- get_and_read(ecoli_acc, "E.coli")

# --------------------
# 2) Number of CDS and total coding DNA
# --------------------
compute_cds_stats <- function(cds) {
  if (is.null(cds)) return(NULL)
  n_cds <- length(cds)
  cds_lengths <- Biostrings::width(cds)
  total_bp <- sum(cds_lengths, na.rm = TRUE)
  return(list(n_cds = n_cds, lengths = cds_lengths, total_bp = total_bp))
}

tetra_stats <- compute_cds_stats(tetra$cds)
ecoli_stats <- compute_cds_stats(ecoli$cds)

# Create summary table
summary_table <- data.frame(
  Organism = c("Tetrasphaera_sp_Soli756", "Escherichia_coli_K12_MG1655"),
  Assembly = c(tetra_acc, ecoli_acc),
  CDS_count = c(ifelse(is.null(tetra_stats), NA, tetra_stats$n_cds),
                ifelse(is.null(ecoli_stats), NA, ecoli_stats$n_cds)),
  Total_CDS_bp = c(ifelse(is.null(tetra_stats), NA, tetra_stats$total_bp),
                   ifelse(is.null(ecoli_stats), NA, ecoli_stats$total_bp))
)
print(summary_table)
write.csv(summary_table, file = file.path(outdir, "cds_summary_table.csv"), row.names = FALSE)

# --------------------
# 3) CDS length distribution, boxplot, mean & median
# --------------------
# Combine lengths into a dataframe for plotting
lengths_df <- data.frame(
  length = c(ifelse(is.null(tetra_stats), numeric(0), tetra_stats$lengths),
             ifelse(is.null(ecoli_stats), numeric(0), ecoli_stats$lengths)),
  organism = c(rep("Tetrasphaera", length(ifelse(is.null(tetra_stats), numeric(0), tetra_stats$lengths))),
               rep("E_coli", length(ifelse(is.null(ecoli_stats), numeric(0), ecoli_stats$lengths))))
)

# Compute mean & median
length_summary <- lengths_df %>% group_by(organism) %>% summarise(
  mean_length = mean(length, na.rm = TRUE),
  median_length = median(length, na.rm = TRUE),
  n = n()
)
print(length_summary)
write.csv(length_summary, file = file.path(outdir, "cds_length_summary.csv"), row.names = FALSE)

# Boxplot
p_box <- ggplot(lengths_df, aes(x = organism, y = length, fill = organism)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "CDS length distribution", x = "", y = "Length (bp)")
print(p_box)
ggsave(filename = file.path(outdir, "cds_length_boxplot.png"), plot = p_box, width = 7, height = 5)

# --------------------
# 4) Base frequency (DNA) and amino-acid frequency (protein)
# --------------------
# Nucleotide frequency over concatenated CDS
count_nucleotides <- function(cds) {
  if (is.null(cds)) return(NULL)
  # alphabetFrequency can sum across a set
  af <- Biostrings::alphabetFrequency(cds, baseOnly = TRUE, as.prob = FALSE)
  sum_af <- colSums(af[, c("A","C","G","T"), drop = FALSE], na.rm = TRUE)
  return(sum_af)
}

tetra_nt <- count_nucleotides(tetra$cds)
ecoli_nt <- count_nucleotides(ecoli$cds)
nt_df <- data.frame(
  Organism = rep(c("Tetrasphaera","E_coli"), each = 4),
  Base = rep(c("A","C","G","T"), times = 2),
  Count = c(tetra_nt, ecoli_nt)
)
print(nt_df)
write.csv(nt_df, file = file.path(outdir, "nucleotide_counts.csv"), row.names = FALSE)

# Protein freq
count_aas <- function(prot) {
  if (is.null(prot)) return(NULL)
  af <- Biostrings::alphabetFrequency(prot, baseOnly = FALSE, as.prob = FALSE)
  aa_cols <- intersect(colnames(af), Biostrings::AMINO_ACID_CODE)
  sum_af <- colSums(af[, aa_cols, drop = FALSE], na.rm = TRUE)
  return(sum_af)
}

tetra_aa <- count_aas(tetra$prot)
ecoli_aa <- count_aas(ecoli$prot)
# Ensure same amino acid order
aa_names <- sort(unique(c(names(tetra_aa), names(ecoli_aa))))
tetra_aa_all <- sapply(aa_names, function(x) ifelse(!is.null(tetra_aa[x]), tetra_aa[x], 0))
ecoli_aa_all <- sapply(aa_names, function(x) ifelse(!is.null(ecoli_aa[x]), ecoli_aa[x], 0))
aa_df <- data.frame(
  Organism = rep(c("Tetrasphaera","E_coli"), each = length(aa_names)),
  AA = rep(aa_names, times = 2),
  Count = c(as.numeric(tetra_aa_all), as.numeric(ecoli_aa_all))
)
print(head(aa_df))
write.csv(aa_df, file = file.path(outdir, "amino_acid_counts.csv"), row.names = FALSE)

# Barplots for nucleotides and amino acids
p_nt <- ggplot(nt_df, aes(x = Base, y = Count, fill = Organism)) +
  geom_col(position = "dodge") + theme_minimal() + labs(title = "Nucleotide counts in CDS", y = "Count")
print(p_nt)
ggsave(filename = file.path(outdir, "nucleotide_counts.png"), plot = p_nt, width = 6, height = 4)

p_aa <- ggplot(aa_df, aes(x = AA, y = Count, fill = Organism)) +
  geom_col(position = "dodge") + theme_minimal() + labs(title = "Amino acid counts in proteins", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_aa)
ggsave(filename = file.path(outdir, "amino_acid_counts.png"), plot = p_aa, width = 10, height = 5)

# --------------------
# 5) Codon usage table (count codons in-frame)
# --------------------
count_codons_from_cds <- function(cds_set) {
  if (is.null(cds_set)) return(NULL)
  codon_counts <- integer(64)
  names(codon_counts) <- as.vector(outer(c("T","C","A","G"), c("T","C","A","G"), paste0))
  # Iterate sequences and extract codons in frame 1
  for (i in seq_along(cds_set)) {
    seq_chr <- as.character(cds_set[[i]])
    L <- nchar(seq_chr)
    max_full <- floor(L / 3) * 3
    if (max_full < 3) next
    cods <- substring(seq_chr, first = seq(1, max_full, by = 3), last = seq(3, max_full, by = 3))
    cods <- toupper(cods)
    tab <- table(cods)
    codon_counts[names(tab)] <- codon_counts[names(tab)] + as.integer(tab)
  }
  return(codon_counts)
}

tetra_codon <- count_codons_from_cds(tetra$cds)
ecoli_codon <- count_codons_from_cds(ecoli$cds)
codon_df <- data.frame(
  Codon = names(tetra_codon),
  Tetrasphaera = as.integer(tetra_codon),
  E_coli = as.integer(ecoli_codon)
)
codon_df <- codon_df %>% mutate(Tet_per_1000 = 1000 * Tetrasphaera / sum(Tetrasphaera, na.rm=TRUE),
                                Eco_per_1000 = 1000 * E_coli / sum(E_coli, na.rm=TRUE))
write.csv(codon_df, file = file.path(outdir, "codon_usage_table.csv"), row.names = FALSE)
codon_long <- codon_df %>% pivot_longer(cols = c("Tetrasphaera","E_coli"), names_to = "Organism", values_to = "Count")
p_codon <- ggplot(codon_long, aes(x = Codon, y = Count, fill = Organism)) +
  geom_col(position = "dodge") + theme_minimal() + labs(title = "Codon counts (per codon)", y = "Count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p_codon)
ggsave(filename = file.path(outdir, "codon_counts.png"), plot = p_codon, width = 12, height = 6)

# --------------------
# 6) Protein k-mers (k = 3..5): over- and under-representation in Tetrasphaera vs E. coli
# --------------------
count_kmers_protein <- function(prot_set, k) {
  if (is.null(prot_set)) return(NULL)
  freq <- integer(0)
  for (i in seq_along(prot_set)) {
    seq_chr <- as.character(prot_set[[i]])
    L <- nchar(seq_chr)
    if (L < k) next
    starts <- seq(1, L - k + 1)
    kmers <- substring(seq_chr, starts, starts + k - 1)
    kmers <- toupper(kmers)
    tab <- table(kmers)
    for (nm in names(tab)) {
      if (!(nm %in% names(freq))) freq[nm] <- 0L
      freq[nm] <- freq[nm] + as.integer(tab[nm])
    }
  }
  return(freq)
}

`%||%` <- function(a,b) if (!is.null(a)) a else b

tetra_kmers <- list(); ecoli_kmers <- list()
for (k in 3:5) {
  message("Counting protein k-mers, k=", k)
  tetra_kmers[[as.character(k)]] <- count_kmers_protein(tetra$prot, k)
  ecoli_kmers[[as.character(k)]] <- count_kmers_protein(ecoli$prot, k)
}

pseudocount <- 1
kmer_results <- list()
for (k in names(tetra_kmers)) {
  tvec <- tetra_kmers[[k]]
  evec <- ecoli_kmers[[k]]
  all_k <- union(names(tvec), names(evec))
  tcounts <- as.numeric(sapply(all_k, function(x) ifelse(!is.null(tvec[x]), tvec[x], 0)))
  ecounts <- as.numeric(sapply(all_k, function(x) ifelse(!is.null(evec[x]), evec[x], 0)))
  fold <- (tcounts + pseudocount) / (ecounts + pseudocount)
  res_df <- data.frame(kmer = all_k, tetra = tcounts, ecoli = ecounts, fold = fold, stringsAsFactors = FALSE)
  over <- res_df %>% arrange(desc(fold)) %>% head(10)
  under <- res_df %>% arrange(fold) %>% head(10)
  kmer_results[[k]] <- list(over = over, under = under, all = res_df)
  write.csv(res_df, file = file.path(outdir, paste0("kmer_counts_k", k, ".csv")), row.names = FALSE)
  write.csv(over, file = file.path(outdir, paste0("top10_over_k", k, ".csv")), row.names = FALSE)
  write.csv(under, file = file.path(outdir, paste0("top10_under_k", k, ".csv")), row.names = FALSE)
  p_over <- ggplot(over, aes(x = reorder(kmer, fold), y = fold)) + geom_col() + coord_flip() +
    labs(title = paste0("Top 10 over-represented k=", k, " in Tetrasphaera vs E.coli"), x = "", y = "Fold-change (Tetra/Ecoli)") +
    theme_minimal()
  ggsave(filename = file.path(outdir, paste0("top10_over_k", k, ".png")), plot = p_over, width = 7, height = 5)
  p_under <- ggplot(under, aes(x = reorder(kmer, fold), y = fold)) + geom_col() + coord_flip() +
    labs(title = paste0("Top 10 under-represented k=", k, " in Tetrasphaera vs E.coli"), x = "", y = "Fold-change (Tetra/Ecoli)") +
    theme_minimal()
  ggsave(filename = file.path(outdir, paste0("top10_under_k", k, ".png")), plot = p_under, width = 7, height = 5)
}

save.image(file = file.path(outdir, "sequence_diversity_results.RData"))
message("All done. Outputs written to: ", normalizePath(outdir))
