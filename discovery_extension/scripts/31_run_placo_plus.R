#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("usage: 31_run_placo_plus.R MANIFEST PAIR_ID VARIANT_OUT SUMMARY_OUT ABS_TOL")
}

manifest_path <- args[[1]]
pair_id <- args[[2]]
variant_out <- args[[3]]
summary_out <- args[[4]]
abs_tol <- as.numeric(args[[5]])
if (!is.finite(abs_tol) || abs_tol <= 0) stop("ABS_TOL must be positive and finite")

manifest <- read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
selected <- manifest[manifest$pair_id == pair_id, , drop = FALSE]
if (nrow(selected) != 1) stop("PAIR_ID must identify exactly one locked manifest row")
if (selected$results_accessed_before_lock != "NO") stop("manifest is not result-free")

source(selected$placo_source_path)
required_functions <- c("var.placo", "cor.pearson", "placo.plus")
if (!all(vapply(required_functions, exists, logical(1), inherits = TRUE))) {
  stop("pinned PLACO+ source lacks required functions")
}

input_path <- selected$merged_genomewide_path
connection <- if (grepl("\\.(gz|bgz)$", input_path)) gzfile(input_path, "rt") else file(input_path, "rt")
on.exit(close(connection), add = TRUE)
dat <- read.delim(connection, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("SNP", "CHR", "BP", "A1", "A2", "Z1", "Z2", "P1", "P2")
if (!all(required %in% names(dat))) stop("merged genome-wide input lacks required columns")
if (nrow(dat) != as.integer(selected$eligible_variant_count)) stop("eligible variant count differs from lock")

z_matrix <- as.matrix(dat[, c("Z1", "Z2")])
p_matrix <- as.matrix(dat[, c("P1", "P2")])
storage.mode(z_matrix) <- "double"
storage.mode(p_matrix) <- "double"
if (any(!is.finite(z_matrix)) || any(!is.finite(p_matrix)) || any(p_matrix < 0 | p_matrix > 1)) {
  stop("non-finite or invalid Z/P values")
}
if (any(z_matrix^2 > 80)) stop("primary PLACO+ input retained Z squared above 80")

marginal_threshold <- 1e-4
var_z <- var.placo(z_matrix, p_matrix, p.threshold = marginal_threshold)
cor_z <- cor.pearson(z_matrix, p_matrix, p.threshold = marginal_threshold, returnMatrix = FALSE)
if (length(var_z) != 2 || any(!is.finite(var_z)) || any(var_z <= 0)) stop("invalid PLACO+ variance estimates")
if (length(cor_z) != 1 || !is.finite(cor_z) || cor_z <= -1 || cor_z >= 1) stop("invalid PLACO+ CorZ")

evaluate_variant <- function(index) {
  tryCatch({
    result <- placo.plus(z_matrix[index, ], VarZ = var_z, CorZ = cor_z, AbsTol = abs_tol)
    c(T_PLACO_PLUS = as.numeric(result$T.placo.plus), P_PLACO_PLUS = as.numeric(result$p.placo.plus), status = 1)
  }, error = function(error) {
    c(T_PLACO_PLUS = NA_real_, P_PLACO_PLUS = NA_real_, status = 0)
  })
}
evaluated <- t(vapply(seq_len(nrow(dat)), evaluate_variant, numeric(3)))
variant_results <- data.frame(
  pair_id = pair_id, SNP = dat$SNP, CHR = dat$CHR, BP = dat$BP,
  A1 = dat$A1, A2 = dat$A2, Z1 = dat$Z1, Z2 = dat$Z2, P1 = dat$P1, P2 = dat$P2,
  T_PLACO_PLUS = evaluated[, "T_PLACO_PLUS"], P_PLACO_PLUS = evaluated[, "P_PLACO_PLUS"],
  analysis_status = ifelse(evaluated[, "status"] == 1, "PLACO_PLUS_COMPLETE", "FAILED_NUMERICAL_INTEGRATION"),
  stringsAsFactors = FALSE
)
dir.create(dirname(variant_out), recursive = TRUE, showWarnings = FALSE)
write.table(variant_results, variant_out, quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

summary <- data.frame(
  pair_id = pair_id, input_variant_count = as.integer(selected$input_variant_count),
  eligible_variant_count = nrow(dat), z2_excluded_count = as.integer(selected$z2_excluded_count),
  VarZ1 = var_z[[1]], VarZ2 = var_z[[2]], CorZ = cor_z,
  marginal_p_threshold = marginal_threshold, analysis_status = "PLACO_PLUS_COMPLETE",
  stringsAsFactors = FALSE
)
dir.create(dirname(summary_out), recursive = TRUE, showWarnings = FALSE)
write.table(summary, summary_out, quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")
cat(sprintf("PLACO_PLUS_PAIR_OK pair=%s variants=%d failures=%d CorZ=%.8g\n", pair_id, nrow(dat), sum(evaluated[, "status"] == 0), cor_z))
