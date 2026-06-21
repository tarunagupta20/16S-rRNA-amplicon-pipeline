# =============================================================================
# dada2_denoise.R — ASV Generation via DADA2 Denoising
# =============================================================================
# Description : Performs quality filtering, denoising, merging, and chimera
#               removal on paired-end 16S rRNA amplicon reads using DADA2.
#               Produces an ASV count table and representative sequences.
#
# Input       : Directory containing trimmed paired-end FASTQ files
#               (output from qc_reads.sh)
#
# Output      : - asv_table.rds       — ASV count matrix (samples x ASVs)
#               - asv_sequences.fasta — Representative ASV sequences
#               - dada2_summary.tsv   — Read tracking table per sample
#               - quality_plots/      — Per-sample quality profiles
#
# Usage       : Rscript scripts/dada2_denoise.R \
#                   --input data/processed/ \
#                   --output results/dada2/ \
#                   --truncF 230 \
#                   --truncR 200
#
# Author      : Taruna Gupta
# GitHub      : github.com/tarunagupta20
# Date        : 2026-06-21
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dada2)
  library(optparse)
  library(ggplot2)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--input",  type="character", help="Path to trimmed FASTQ directory"),
  make_option("--output", type="character", default="results/dada2", help="Output directory"),
  make_option("--truncF", type="integer",   default=230, help="Truncation length for forward reads (default: 230)"),
  make_option("--truncR", type="integer",   default=200, help="Truncation length for reverse reads (default: 200)"),
  make_option("--maxEE",  type="numeric",   default=2,   help="Max expected errors per read (default: 2)"),
  make_option("--threads",type="integer",   default=4,   help="Number of threads (default: 4)")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ── Validate input ─────────────────────────────────────────────────────────────
if (is.null(opt$input)) stop("ERROR: --input directory is required")
if (!dir.exists(opt$input)) stop(paste("ERROR: Input directory not found:", opt$input))

dir.create(opt$output, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(opt$output, "quality_plots"), showWarnings=FALSE)

cat("=============================================\n")
cat("  DADA2 Denoising Pipeline\n")
cat(paste("  Input:   ", opt$input, "\n"))
cat(paste("  Output:  ", opt$output, "\n"))
cat(paste("  truncF:  ", opt$truncF, "\n"))
cat(paste("  truncR:  ", opt$truncR, "\n"))
cat("=============================================\n\n")

# ── STEP 1: Get file paths ─────────────────────────────────────────────────────
cat("[STEP 1] Loading FASTQ files...\n")

fnFs <- sort(list.files(opt$input, pattern="_R1_trimmed.fastq.gz", full.names=TRUE))
fnRs <- sort(list.files(opt$input, pattern="_R2_trimmed.fastq.gz", full.names=TRUE))

if (length(fnFs) == 0) stop("ERROR: No R1 trimmed FASTQ files found in input directory")
if (length(fnFs) != length(fnRs)) stop("ERROR: Mismatched number of R1 and R2 files")

# Extract sample names from filenames
sample_names <- sub("_R1_trimmed.fastq.gz", "", basename(fnFs))
cat(paste("  Found", length(sample_names), "samples\n\n"))

# ── STEP 2: Quality profiles ───────────────────────────────────────────────────
cat("[STEP 2] Plotting quality profiles...\n")

pdf(file.path(opt$output, "quality_plots", "forward_quality.pdf"))
print(plotQualityProfile(fnFs[1:min(4, length(fnFs))]))
dev.off()

pdf(file.path(opt$output, "quality_plots", "reverse_quality.pdf"))
print(plotQualityProfile(fnRs[1:min(4, length(fnRs))]))
dev.off()

cat("  Quality plots saved\n\n")

# ── STEP 3: Filter and trim ────────────────────────────────────────────────────
cat("[STEP 3] Quality filtering and trimming...\n")

filt_dir <- file.path(opt$output, "filtered")
filtFs <- file.path(filt_dir, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(sample_names, "_R_filt.fastq.gz"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen  = c(opt$truncF, opt$truncR),
  maxN      = 0,
  maxEE     = c(opt$maxEE, opt$maxEE),
  truncQ    = 2,
  rm.phix   = TRUE,
  compress  = TRUE,
  multithread = opt$threads
)

cat(paste("  Reads passing filter:", sum(out[,2]), "/", sum(out[,1]), "\n\n"))

# ── STEP 4: Learn error rates ──────────────────────────────────────────────────
cat("[STEP 4] Learning error rates...\n")

errF <- learnErrors(filtFs, multithread=opt$threads)
errR <- learnErrors(filtRs, multithread=opt$threads)

pdf(file.path(opt$output, "quality_plots", "error_rates.pdf"))
print(plotErrors(errF, nominalQ=TRUE))
dev.off()

cat("  Error rates learned\n\n")

# ── STEP 5: Denoise (core DADA2 step) ─────────────────────────────────────────
cat("[STEP 5] Denoising with DADA2...\n")

dadaFs <- dada(filtFs, err=errF, multithread=opt$threads)
dadaRs <- dada(filtRs, err=errR, multithread=opt$threads)

# ── STEP 6: Merge paired reads ─────────────────────────────────────────────────
cat("[STEP 6] Merging paired reads...\n")

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=FALSE)

# ── STEP 7: Make ASV table ─────────────────────────────────────────────────────
cat("[STEP 7] Building ASV count table...\n")

seqtab <- makeSequenceTable(mergers)
cat(paste("  ASVs before chimera removal:", ncol(seqtab), "\n"))

# ── STEP 8: Remove chimeras ────────────────────────────────────────────────────
cat("[STEP 8] Removing chimeric sequences...\n")

seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=opt$threads)
cat(paste("  ASVs after chimera removal: ", ncol(seqtab_nochim), "\n"))
cat(paste("  Reads retained:             ",
          round(sum(seqtab_nochim)/sum(seqtab)*100, 1), "%\n\n"))

# ── STEP 9: Track reads through pipeline ──────────────────────────────────────
cat("[STEP 9] Generating read tracking summary...\n")

get_n <- function(x) sum(getUniques(x))

track <- cbind(
  out,
  sapply(dadaFs,  get_n),
  sapply(dadaRs,  get_n),
  sapply(mergers, get_n),
  rowSums(seqtab_nochim)
)
colnames(track) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
rownames(track) <- sample_names

write.table(track,
            file      = file.path(opt$output, "dada2_summary.tsv"),
            sep       = "\t",
            quote     = FALSE,
            col.names = NA)

# ── STEP 10: Save outputs ──────────────────────────────────────────────────────
cat("[STEP 10] Saving outputs...\n")

saveRDS(seqtab_nochim, file.path(opt$output, "asv_table.rds"))

# Save ASV sequences as FASTA
asv_seqs    <- colnames(seqtab_nochim)
asv_headers <- paste0(">ASV", seq_along(asv_seqs))
fasta_lines <- c(rbind(asv_headers, asv_seqs))
writeLines(fasta_lines, file.path(opt$output, "asv_sequences.fasta"))

# ── Summary ────────────────────────────────────────────────────────────────────
cat("=============================================\n")
cat("  DADA2 COMPLETE\n")
cat(paste("  ASV table:    ", file.path(opt$output, "asv_table.rds"),       "\n"))
cat(paste("  FASTA:        ", file.path(opt$output, "asv_sequences.fasta"), "\n"))
cat(paste("  Read summary: ", file.path(opt$output, "dada2_summary.tsv"),   "\n"))
cat(paste("  Quality plots:", file.path(opt$output, "quality_plots/"),      "\n"))
cat("=============================================\n")
