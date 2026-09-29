#!/usr/bin/env Rscript
# Extract sparse UKB-reference LD only for frozen PLACO candidate SNPs.
# LAVA 0.1.5 internals are used because the reference is its native .bcor format.
suppressPackageStartupMessages(library(LAVA))
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3L) stop("usage: extract_placo_ld.R <reference-prefix> <candidate-tsv> <out-dir>")
prefix <- args[[1L]]
candidate_file <- args[[2L]]
out_dir <- args[[3L]]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
v <- as.character(utils::packageVersion("LAVA"))
if (v != "0.1.5") stop(paste0("Expected pinned LAVA 0.1.5; got ", v))
candidates <- utils::read.delim(candidate_file, stringsAsFactors=FALSE, check.names=FALSE)
required <- c("SNP", "CHR", "BP")
if (!all(required %in% names(candidates))) stop("candidate table requires SNP, CHR, BP")
candidates$SNP <- tolower(candidates$SNP)
candidates$CHR <- as.integer(candidates$CHR)
candidates$BP <- as.integer(candidates$BP)
if (anyNA(candidates[, required])) stop("candidate identifiers/coordinates must be nonmissing")
if (any(!candidates$CHR %in% 1:22)) stop("only autosomal candidates are supported")
ref_rows <- list()
ld_rows <- list()
ref_n <- numeric()
for (chr in sort(unique(candidates$CHR))) {
  requested <- unique(candidates$SNP[candidates$CHR == chr])
  ref <- LAVA:::load.reference(list(prefix=prefix, mode="ld", chromosomes=chr))
  ref_n <- c(ref_n, ref$sample.size)
  info <- ref$snp.info[ref$snp.info$SNP %in% requested, , drop=FALSE]
  if (nrow(info) > 0L) {
    info$CHR <- as.integer(info$CHR)
    info$POS <- as.integer(info$POS)
    ref_rows[[length(ref_rows)+1L]] <- info[, c("SNP", "CHR", "POS", "A1", "A2")]
    ld <- LAVA:::read.ld(ref, info$SNP)
    mat <- as.matrix(ld$ld)
    if (nrow(mat) != ncol(mat) || !all(is.finite(mat))) stop(paste("invalid LD matrix on chromosome", chr))
    if (max(abs(mat - t(mat))) > 1e-8) stop(paste("asymmetric LD matrix on chromosome", chr))
    if (any(abs(diag(mat) - 1) > 1e-6)) stop(paste("LD diagonal differs from one on chromosome", chr))
    ij <- which(upper.tri(mat, diag=TRUE), arr.ind=TRUE)
    ld_rows[[length(ld_rows)+1L]] <- data.frame(
      SNP1=rownames(mat)[ij[,1L]], SNP2=colnames(mat)[ij[,2L]], R=mat[ij],
      stringsAsFactors=FALSE)
  }
}
info_out <- if (length(ref_rows)) do.call(rbind, ref_rows) else data.frame(SNP=character(), CHR=integer(), POS=integer(), A1=character(), A2=character())
ld_out <- if (length(ld_rows)) do.call(rbind, ld_rows) else data.frame(SNP1=character(), SNP2=character(), R=numeric())
utils::write.table(info_out, file=file.path(out_dir, "reference_matches.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
utils::write.table(ld_out, file=file.path(out_dir, "ld_pairs.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
writeLines(format(as.numeric(max(ref_n)), scientific=FALSE, trim=TRUE),
           file.path(out_dir, "reference_sample_size.txt"))
