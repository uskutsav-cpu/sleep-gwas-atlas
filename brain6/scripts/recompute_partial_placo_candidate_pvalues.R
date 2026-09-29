#!/usr/bin/env Rscript
# Recompute candidate PLACO+ p-values with the frozen upstream implementation.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("usage: recompute_partial_placo_candidate_pvalues.R input.tsv parameters.tsv placo_source.R output.tsv")
if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required")
data.table::fread(args[[1]]) -> d
parameters <- data.table::fread(args[[2]])
source(args[[3]])
if (!all(c("pair_id", "Z1", "Z2") %in% names(d))) stop("candidate input schema mismatch")
if (!all(c("pair_id", "VarZ1", "VarZ2", "CorZ", "AbsTol") %in% names(parameters))) stop("parameter schema mismatch")
if (anyDuplicated(parameters$pair_id)) stop("duplicate pair parameters")
parameter_index <- match(d$pair_id, parameters$pair_id)
if (anyNA(parameter_index)) stop("missing pair parameters")
d[, P_PLACO_RECOMPUTED := NA_real_]
for (i in seq_len(nrow(d))) {
  p <- parameters[parameter_index[[i]]]
  z <- c(d$Z1[[i]], d$Z2[[i]])
  d$P_PLACO_RECOMPUTED[[i]] <- if (any(z == 0)) 1.0 else
    placo.plus(z, VarZ = c(p$VarZ1, p$VarZ2), CorZ = p$CorZ,
               AbsTol = p$AbsTol)$p.placo.plus
}
data.table::fwrite(d, args[[4]], sep = "\t", quote = FALSE, na = "NA")
