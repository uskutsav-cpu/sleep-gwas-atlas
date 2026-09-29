#!/usr/bin/env Rscript
# Extract only candidate-to-lead signed LD edges for the range-filtered PLACO clumps.
# Raw values are preserved; this script does not clip correlations.
suppressPackageStartupMessages(library(LAVA))
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3L) stop("usage: extract_placo_candidate_lead_ld.R <reference-prefix> <candidate-tsv> <out-dir>")
prefix <- args[[1L]]
candidate_file <- args[[2L]]
out_dir <- args[[3L]]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
candidates <- utils::read.delim(candidate_file, stringsAsFactors=FALSE, check.names=FALSE)
required <- c("pair_id", "SNP", "CHR", "BP", "lead_SNP")
if (!all(required %in% names(candidates))) stop("candidate table requires pair_id, SNP, CHR, BP, lead_SNP")
candidates$SNP <- tolower(candidates$SNP)
candidates$lead_SNP <- tolower(candidates$lead_SNP)
candidates$CHR <- as.integer(candidates$CHR)
candidates$BP <- as.integer(candidates$BP)
if (anyNA(candidates[, required]) || any(!candidates$CHR %in% 1:22)) stop("invalid candidate identifiers/coordinates")
if (as.character(utils::packageVersion("LAVA")) != "0.1.5") stop("Expected pinned LAVA 0.1.5")
edge_rows <- list()
ref_rows <- list()
sample_sizes <- numeric()
for (chr in sort(unique(candidates$CHR))) {
  chr_rows <- candidates[candidates$CHR == chr, , drop=FALSE]
  requested <- unique(c(chr_rows$SNP, chr_rows$lead_SNP))
  ref <- LAVA:::load.reference(list(prefix=prefix, mode="ld", chromosomes=chr))
  sample_sizes <- c(sample_sizes, ref$sample.size)
  info <- ref$snp.info[ref$snp.info$SNP %in% requested, , drop=FALSE]
  if (!nrow(info)) next
  info$CHR <- as.integer(info$CHR)
  info$POS <- as.integer(info$POS)
  info$SNP <- tolower(info$SNP)
  ref_rows[[length(ref_rows)+1L]] <- info[, c("SNP", "CHR", "POS", "A1", "A2")]
  common_ids <- info$SNP
  ld <- LAVA:::read.ld(ref, common_ids)
  mat <- as.matrix(ld$ld)
  if (nrow(mat) != ncol(mat) || !all(is.finite(mat))) stop(paste("invalid LD matrix on chromosome", chr))
  if (max(abs(mat - t(mat))) > 1e-8) stop(paste("asymmetric LD matrix on chromosome", chr))
  if (any(abs(diag(mat) - 1) > 1e-6)) stop(paste("LD diagonal differs from one on chromosome", chr))
  row_ids <- tolower(rownames(mat)); col_ids <- tolower(colnames(mat))
  if (!setequal(row_ids, common_ids) || !setequal(col_ids, common_ids)) stop(paste("LD IDs differ from reference info on chromosome", chr))
  mat <- mat[match(common_ids, row_ids), match(common_ids, col_ids), drop=FALSE]
  for (i in seq_len(nrow(chr_rows))) {
    snp <- chr_rows$SNP[i]; lead <- chr_rows$lead_SNP[i]
    row_idx <- match(snp, common_ids); col_idx <- match(lead, common_ids)
    if (is.na(row_idx) || is.na(col_idx)) next
    edge_rows[[length(edge_rows)+1L]] <- data.frame(
      pair_id=chr_rows$pair_id[i], SNP=snp, lead_SNP=lead, CHR=chr,
      R=mat[row_idx, col_idx], stringsAsFactors=FALSE)
  }
}
edge_out <- if (length(edge_rows)) do.call(rbind, edge_rows) else data.frame(pair_id=character(), SNP=character(), lead_SNP=character(), CHR=integer(), R=numeric())
ref_out <- if (length(ref_rows)) unique(do.call(rbind, ref_rows)) else data.frame(SNP=character(), CHR=integer(), POS=integer(), A1=character(), A2=character())
utils::write.table(edge_out, file=file.path(out_dir, "candidate_lead_ld.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
utils::write.table(ref_out, file=file.path(out_dir, "reference_alleles.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
writeLines(format(as.numeric(max(sample_sizes)), scientific=FALSE, trim=TRUE), file.path(out_dir, "reference_sample_size.txt"))
