#!/usr/bin/env Rscript
# Extract raw LAVA-reference LD only for candidate pairs in the frozen clumping window.
# Values outside [-1, 1] are retained for QC and are never clamped here.
suppressPackageStartupMessages(library(LAVA))
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4L) stop("usage: extract_placo_windowed_ld.R <reference-prefix> <candidate-tsv> <out-dir> <window-bp>")
prefix <- args[[1L]]
candidate_file <- args[[2L]]
out_dir <- args[[3L]]
window_bp <- as.integer(args[[4L]])
if (is.na(window_bp) || window_bp < 1L) stop("window-bp must be a positive integer")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
candidates <- utils::read.delim(candidate_file, stringsAsFactors=FALSE, check.names=FALSE)
required <- c("SNP", "CHR", "BP")
if (!all(required %in% names(candidates))) stop("candidate table requires SNP, CHR, BP")
candidates$SNP <- tolower(candidates$SNP)
candidates$CHR <- as.integer(candidates$CHR)
candidates$BP <- as.integer(candidates$BP)
if (anyNA(candidates[, required]) || any(!candidates$CHR %in% 1:22)) stop("invalid candidate identifiers/coordinates")
if (anyDuplicated(paste(candidates$CHR, candidates$SNP, sep=":"))) stop("candidate SNPs must be unique by chromosome")
if (as.character(utils::packageVersion("LAVA")) != "0.1.5") stop("Expected pinned LAVA 0.1.5")
ref_rows <- list()
ld_rows <- list()
ref_n <- numeric()
for (chr in sort(unique(candidates$CHR))) {
  requested <- candidates$SNP[candidates$CHR == chr]
  ref <- LAVA:::load.reference(list(prefix=prefix, mode="ld", chromosomes=chr))
  ref_n <- c(ref_n, ref$sample.size)
  info <- ref$snp.info[ref$snp.info$SNP %in% requested, , drop=FALSE]
  if (!nrow(info)) next
  info$CHR <- as.integer(info$CHR)
  info$POS <- as.integer(info$POS)
  info$SNP <- tolower(info$SNP)
  pos <- setNames(info$POS, info$SNP)
  ref_rows[[length(ref_rows)+1L]] <- info[, c("SNP", "CHR", "POS", "A1", "A2")]
  ld <- LAVA:::read.ld(ref, info$SNP)
  mat <- as.matrix(ld$ld)
  if (nrow(mat) != ncol(mat) || !all(is.finite(mat))) stop(paste("invalid LD matrix on chromosome", chr))
  if (max(abs(mat - t(mat))) > 1e-8) stop(paste("asymmetric LD matrix on chromosome", chr))
  if (any(abs(diag(mat) - 1) > 1e-6)) stop(paste("LD diagonal differs from one on chromosome", chr))
  names_l <- tolower(rownames(mat)); names_r <- tolower(colnames(mat))
  if (!setequal(names_l, names(pos)) || !setequal(names_r, names(pos))) stop(paste("LD IDs differ from reference info on chromosome", chr))
  mat <- mat[match(names(pos), names_l), match(names(pos), names_r), drop=FALSE]
  pairs <- vector("list", nrow(mat))
  for (i in seq_len(nrow(mat))) {
    # Candidate order in reference metadata need not be genomic.
    js <- which(seq_len(nrow(mat)) >= i & abs(pos - pos[i]) <= window_bp)
    if (length(js)) pairs[[i]] <- data.frame(SNP1=names(pos)[i], SNP2=names(pos)[js], R=mat[i, js], stringsAsFactors=FALSE)
  }
  pairs <- Filter(Negate(is.null), pairs)
  if (length(pairs)) ld_rows[[length(ld_rows)+1L]] <- do.call(rbind, pairs)
}
info_out <- if (length(ref_rows)) do.call(rbind, ref_rows) else data.frame(SNP=character(), CHR=integer(), POS=integer(), A1=character(), A2=character())
ld_out <- if (length(ld_rows)) do.call(rbind, ld_rows) else data.frame(SNP1=character(), SNP2=character(), R=numeric())
utils::write.table(info_out, file=file.path(out_dir, "reference_matches.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
utils::write.table(ld_out, file=file.path(out_dir, "ld_pairs_within_window.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
writeLines(format(as.numeric(max(ref_n)), scientific=FALSE, trim=TRUE), file.path(out_dir, "reference_sample_size.txt"))
