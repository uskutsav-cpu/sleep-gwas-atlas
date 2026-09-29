#!/usr/bin/env Rscript
# Independently reload every preserved PLACO LD range exception via pinned LAVA.
# This is a reader replay only; it does not modify candidate/clumping outputs.
suppressPackageStartupMessages(library(LAVA))
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3L) stop("usage: replay_placo_ld_range_exceptions.R <reference-prefix> <exceptions.tsv> <output.tsv>")
prefix <- args[[1L]]
exceptions_file <- args[[2L]]
output_file <- args[[3L]]
if (as.character(utils::packageVersion("LAVA")) != "0.1.5") stop("Expected pinned LAVA 0.1.5")
exceptions <- utils::read.delim(exceptions_file, stringsAsFactors=FALSE, check.names=FALSE)
required <- c("CHR", "SNP1", "SNP2", "lava_r_raw")
if (!all(required %in% names(exceptions))) stop("Exception table schema is incomplete")
exceptions$CHR <- as.integer(exceptions$CHR)
exceptions$SNP1 <- tolower(exceptions$SNP1)
exceptions$SNP2 <- tolower(exceptions$SNP2)
exceptions$lava_r_raw <- as.numeric(exceptions$lava_r_raw)
if (anyNA(exceptions[, required]) || any(!exceptions$CHR %in% 1:22) ||
    any(!is.finite(exceptions$lava_r_raw))) stop("Invalid exception rows")
edges <- unique(exceptions[, c("CHR", "SNP1", "SNP2", "lava_r_raw")])
if (anyDuplicated(paste(edges$CHR, edges$SNP1, edges$SNP2, sep=":"))) {
  stop("Duplicate raw LD values disagree across pair-specific exception rows")
}
summary <- list()
for (chr in sort(unique(edges$CHR))) {
  chr_edges <- edges[edges$CHR == chr, , drop=FALSE]
  requested <- unique(c(chr_edges$SNP1, chr_edges$SNP2))
  ref <- LAVA:::load.reference(list(prefix=prefix, mode="ld", chromosomes=chr))
  ids <- ref$snp.info$SNP[ref$snp.info$SNP %in% requested]
  if (length(ids) != length(requested) || !setequal(ids, requested)) {
    stop(paste("Exception SNPs do not map exactly to reference rows on chromosome", chr))
  }
  ld <- LAVA:::read.ld(ref, ids)
  mat <- as.matrix(ld$ld)
  if (!identical(rownames(mat), ids) || !identical(colnames(mat), ids)) {
    stop(paste("LAVA changed requested SNP order on chromosome", chr))
  }
  if (any(!is.finite(mat)) || max(abs(mat - t(mat))) > 1e-8 ||
      max(abs(diag(mat) - 1)) > 1e-6) stop(paste("Invalid replay matrix on chromosome", chr))
  i <- match(chr_edges$SNP1, ids)
  j <- match(chr_edges$SNP2, ids)
  replay <- mat[cbind(i, j)]
  delta <- abs(replay - chr_edges$lava_r_raw)
  if (any(delta > 1e-12)) stop(paste("Exception replay differs from archived extraction on chromosome", chr))
  summary[[length(summary)+1L]] <- data.frame(
    chromosome=chr, exception_edges=nrow(chr_edges), unique_exception_snps=length(ids),
    replay_matches=nrow(chr_edges), max_abs_replay_delta=max(delta),
    max_abs_r_in_exception_submatrix=max(abs(mat)),
    count_abs_r_gt_1=sum(abs(mat) > 1), minimum_eigenvalue=min(eigen(mat, symmetric=TRUE, only.values=TRUE)$values),
    diagonal_max_abs_deviation=max(abs(diag(mat)-1)), symmetry_max_abs=max(abs(mat-t(mat))),
    stringsAsFactors=FALSE)
}
result <- do.call(rbind, summary)
if (sum(result$exception_edges) != nrow(edges) ||
    sum(result$replay_matches) != nrow(edges)) stop("Replay did not cover all unique exceptions")
dir.create(dirname(output_file), recursive=TRUE, showWarnings=FALSE)
utils::write.table(result, output_file, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
cat(sprintf("LAVA_EXCEPTION_REPLAY\tstatus=PASS\tunique_edges=%d\tchromosomes=%d\tmax_abs_delta=%.17g\n",
            nrow(edges), nrow(result), max(result$max_abs_replay_delta)))
