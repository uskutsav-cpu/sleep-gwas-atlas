#!/usr/bin/env Rscript

# One immutable Track B pair x official block per fresh R process.  Model and
# prior failures are serialized as retained QC rows; this engine never invokes
# coloc.abf and never changes parameters or retries a nonconverged model.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: 150_run_track_b_susie_coloc.R TASK_TSV EMPTY_OUTPUT_DIRECTORY")
}
task_path <- args[[1L]]
output_directory <- args[[2L]]

task_fields <- c(
  "analysis_id", "locus_entry_id", "pair_id", "family_role", "trait1", "trait2",
  "CHR", "START", "STOP", "ld_block_id", "variant_count", "summary1_path",
  "summary1_sha256", "summary2_path", "summary2_sha256", "variant_order_path",
  "variant_order_sha256", "signed_ld_path", "signed_ld_sha256", "ld_qc_path",
  "ld_qc_sha256", "ld_validation_path", "ld_validation_sha256",
  "ld_runtime_package_lock_path", "ld_runtime_package_lock_sha256",
  "ld_runtime_attestation_path", "ld_runtime_attestation_sha256",
  "runtime_package_lock_path", "runtime_package_lock_sha256",
  "execution_amendment_sha256", "trait1_type", "trait1_ncase", "trait1_ncontrol",
  "trait1_case_fraction", "trait1_sdY", "trait1_scalar_N", "trait1_scalar_N_rule",
  "trait1_per_snp_N_sha256", "trait1_per_snp_N_count", "trait1_per_snp_N_min",
  "trait1_per_snp_N_q1", "trait1_per_snp_N_median", "trait1_per_snp_N_q3",
  "trait1_per_snp_N_max", "trait1_per_snp_N_max_to_min_ratio",
  "trait1_per_snp_N_fraction_below_90pct_max",
  "trait1_per_snp_N_fraction_below_50pct_max", "trait1_N_dispersion_status",
  "trait1_INFO_status", "trait2_type", "trait2_ncase",
  "trait2_ncontrol", "trait2_case_fraction", "trait2_sdY", "trait2_scalar_N",
  "trait2_scalar_N_rule", "trait2_per_snp_N_sha256", "trait2_per_snp_N_count",
  "trait2_per_snp_N_min", "trait2_per_snp_N_q1", "trait2_per_snp_N_median",
  "trait2_per_snp_N_q3", "trait2_per_snp_N_max",
  "trait2_per_snp_N_max_to_min_ratio",
  "trait2_per_snp_N_fraction_below_90pct_max",
  "trait2_per_snp_N_fraction_below_50pct_max", "trait2_N_dispersion_status",
  "trait2_INFO_status",
  "reference_id", "reference_state", "reference_sample_size", "prior_method",
  "maximum_causal_signals", "credible_set_coverage", "minimum_absolute_correlation",
  "maximum_iterations", "estimate_residual_variance", "p1", "p2", "p12_grid",
  "single_signal_fallback", "maximum_locus_variants", "locus_splitting",
  "variant_thinning", "lead_centered_truncation", "claim_limit"
)
summary_fields <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "EAF", "INFO", "N")
order_fields <- c("SNP", "CHR", "BP", "A1", "A2")
raw_variant_fields <- c(
  "analysis_id", "pair_id", "family_role", "locus_entry_id", "locus_id", "CHR",
  "START", "STOP", "trait_id", "trait_role", "SNP", "BP", "A1", "A2", "BETA",
  "SE", "P", "EAF", "INFO", "N", "prior_role", "prior_weight", "PIP",
  "credible_set_ids", "credible_set_sizes", "credible_set_requested_coverage",
  "credible_set_achieved_coverage", "credible_set_min_abs_corr",
  "credible_set_mean_abs_corr", "max_alpha_component", "model_converged", "niter",
  "rss_ld_s", "kriging_allele_switch_outlier", "diagnostic_status", "error",
  "claim_limit"
)
credible_fields <- c(
  "analysis_id", "pair_id", "family_role", "locus_entry_id", "trait_id",
  "trait_role", "signal_id", "component_index", "lead_snp", "lead_pip",
  "credible_set_size", "credible_set_snps", "requested_coverage",
  "achieved_coverage", "min_abs_corr", "mean_abs_corr", "median_abs_corr",
  "cs_log10bf", "model_converged", "analysis_status", "claim_limit"
)
diagnostic_fields <- c(
  "analysis_id", "pair_id", "family_role", "locus_entry_id", "dataset_role",
  "trait_id", "variant_count", "ld_variant_count", "allele_match_count",
  "model_converged", "niter", "credible_set_count", "max_pip", "rss_ld_s",
  "kriging_allele_switch_outlier_count", "kriging_allele_switch_outliers",
  "ld_symmetry_max_abs", "ld_diagonal_max_abs", "ld_minimum_eigenvalue",
  "reference_sample_size", "diagnostic_status", "error"
)
finemap_qc_fields <- c(
  "analysis_id", "pair_id", "family_role", "locus_entry_id", "trait_id",
  "trait_role", "variant_count", "model_converged", "niter", "credible_set_count",
  "max_pip", "rss_ld_s", "kriging_allele_switch_outlier_count",
  "ld_symmetry_max_abs", "ld_diagonal_max_abs", "ld_minimum_eigenvalue",
  "reference_sample_size", "analysis_status", "diagnostic_status", "error"
)
raw_coloc_fields <- c(
  "analysis_id", "pair_id", "family_role", "locus_entry_id", "locus_id", "CHR",
  "START", "STOP", "trait1", "trait2", "signal1", "signal2", "coloc_method", "p1",
  "p2", "p12", "prior_role", "nsnps", "PP_H0", "PP_H1", "PP_H2", "PP_H3",
  "PP_H4", "PP_H4_over_PP_H3", "top_shared_variant",
  "top_shared_variant_PP_H4", "fine_mapping_qc", "ld_qc", "engine_status", "error",
  "claim_limit"
)
engine_status_fields <- c(
  "analysis_id", "pair_id", "locus_entry_id", "engine_status", "trait1_status",
  "trait2_status", "coloc_status", "error"
)
sample_size_fields <- c(
  "analysis_id", "pair_id", "family_role", "locus_entry_id", "dataset_role",
  "trait_id", "per_snp_N_sha256", "per_snp_N_count", "per_snp_N_min",
  "per_snp_N_q1", "per_snp_N_median", "per_snp_N_q3", "per_snp_N_max",
  "per_snp_N_max_to_min_ratio", "per_snp_N_fraction_below_90pct_max",
  "per_snp_N_fraction_below_50pct_max", "scalar_N", "scalar_N_rule",
  "dispersion_status", "claim_cap"
)
runtime_lock_fields <- c(
  "record_type", "name", "version", "path", "file_count", "bytes", "sha256"
)
runtime_attestation_fields <- c(
  "analysis_id", "locus_entry_id", "execution_amendment_sha256",
  "runtime_package_lock_sha256", "precheck_status", "postcheck_status",
  "loaded_namespace_count"
)

empty_frame <- function(fields) {
  as.data.frame(
    setNames(replicate(length(fields), character(0), simplify = FALSE), fields),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

sha256_file <- function(path) {
  command <- Sys.which("sha256sum")
  command_args <- path
  if (!nzchar(command)) {
    command <- Sys.which("shasum")
    command_args <- c("-a", "256", path)
  }
  if (!nzchar(command)) stop("no SHA-256 utility is available")
  output <- system2(command, command_args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("SHA-256 utility failed")
  token <- strsplit(trimws(output[[1L]]), "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", token)) stop("invalid SHA-256 output")
  token
}

verify_file <- function(path, expected, label) {
  if (!file.exists(path) || !grepl("^[0-9a-f]{64}$", expected) || sha256_file(path) != expected) {
    stop(sprintf("locked %s is missing or hash-mismatched", label))
  }
}

runtime_namespace_closure <- c(
  "base", "cli", "coloc", "compiler", "crayon", "data.table", "datasets",
  "dplyr", "farver", "generics", "ggplot2", "glue", "graphics", "grDevices",
  "grid", "gridExtra", "gtable", "irlba", "lattice", "lifecycle", "magrittr",
  "Matrix", "matrixStats", "methods", "mixsqp", "pillar", "pkgconfig", "plyr",
  "R6", "RColorBrewer", "Rcpp", "reshape", "rlang", "S7", "scales", "stats",
  "susieR", "tibble", "tidyselect", "tools", "utils", "vctrs", "viridis",
  "viridisLite"
)

package_tree_identity <- function(package_path) {
  relative <- sort(list.files(
    package_path, all.files = TRUE, recursive = TRUE, full.names = FALSE,
    include.dirs = FALSE, no.. = TRUE
  ), method = "radix")
  full <- file.path(package_path, relative)
  info <- file.info(full)
  keep <- !is.na(info$isdir) & !info$isdir
  relative <- relative[keep]
  full <- full[keep]
  info <- info[keep, , drop = FALSE]
  if (!length(relative) || any(grepl("[\t\r\n]", relative))) {
    stop("runtime package tree is empty or has an unsafe path")
  }
  links <- Sys.readlink(full)
  if (any(!is.na(links) & nzchar(links))) stop("runtime package tree contains a symbolic link")
  hashes <- vapply(full, sha256_file, character(1L))
  sizes <- format(info$size, scientific = FALSE, trim = TRUE)
  manifest <- paste0(paste(relative, sizes, hashes, sep = "\t", collapse = "\n"), "\n")
  temporary <- tempfile(pattern = ".runtime-tree.", tmpdir = output_directory)
  on.exit(unlink(temporary), add = TRUE)
  connection <- file(temporary, open = "wb")
  writeBin(charToRaw(enc2utf8(manifest)), connection)
  close(connection)
  list(
    file_count = length(relative), bytes = sum(info$size), sha256 = sha256_file(temporary)
  )
}

verify_runtime_lock <- function(phase, require_loaded_closure = FALSE) {
  lock <- read.delim(
    task$runtime_package_lock_path, stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = "character"
  )
  if (!identical(names(lock), runtime_lock_fields) || nrow(lock) != length(runtime_namespace_closure) + 1L) {
    stop("runtime package lock schema/closure count drifted")
  }
  executable <- lock[lock$record_type == "R_EXECUTABLE", , drop = FALSE]
  packages <- lock[lock$record_type == "PACKAGE_TREE", , drop = FALSE]
  if (
    nrow(executable) != 1L || executable$name != "Rscript" || executable$version != "R-4.3.3" ||
    executable$file_count != "1" || normalizePath(executable$path, winslash = "/", mustWork = TRUE) !=
      normalizePath(".r-env/bin/Rscript", winslash = "/", mustWork = TRUE) ||
    as.character(file.info(executable$path)$size) != executable$bytes ||
    sha256_file(executable$path) != executable$sha256 ||
    !identical(packages$name, runtime_namespace_closure)
  ) {
    stop(sprintf("%s runtime executable or namespace family differs from lock", phase))
  }
  for (index in seq_len(nrow(packages))) {
    row <- packages[index, , drop = FALSE]
    located <- find.package(row$name, quiet = TRUE)
    if (
      !nzchar(located) || normalizePath(located, winslash = "/", mustWork = TRUE) !=
        normalizePath(row$path, winslash = "/", mustWork = TRUE) ||
      as.character(utils::packageDescription(
        row$name, lib.loc = dirname(row$path), fields = "Version"
      )) != row$version
    ) stop(sprintf("%s loaded-package path/version mismatch for %s", phase, row$name))
    # At PRE, validate the actual namespace path for every package already
    # loaded by R itself. At POST, every locked namespace must be loaded and is
    # checked here, followed by the exact-closure assertion below.
    if (row$name %in% loadedNamespaces()) {
      loaded_path <- if (row$name == "base") {
        file.path(R.home("library"), "base")
      } else {
        getNamespaceInfo(asNamespace(row$name), "path")
      }
      if (normalizePath(loaded_path, winslash = "/", mustWork = TRUE) !=
          normalizePath(row$path, winslash = "/", mustWork = TRUE)) {
        stop(sprintf("%s actual namespace path mismatch for %s", phase, row$name))
      }
    } else if (require_loaded_closure) {
      stop(sprintf("%s locked namespace is not actually loaded: %s", phase, row$name))
    }
    tree <- package_tree_identity(row$path)
    if (
      as.character(tree$file_count) != row$file_count ||
      as.character(tree$bytes) != row$bytes || tree$sha256 != row$sha256
    ) stop(sprintf("%s installed-package byte mismatch for %s", phase, row$name))
  }
  if (require_loaded_closure && !identical(sort(loadedNamespaces()), sort(runtime_namespace_closure))) {
    stop(sprintf("%s actually loaded namespace closure differs from lock", phase))
  }
  invisible(lock)
}

write_tsv_no_replace <- function(data, name, fields) {
  path <- file.path(output_directory, name)
  if (file.exists(path)) stop(sprintf("refusing to replace engine output %s", name))
  if (!identical(names(data), fields)) stop(sprintf("internal schema mismatch for %s", name))
  temporary <- tempfile(pattern = paste0(".", name, "."), tmpdir = output_directory)
  on.exit(unlink(temporary), add = TRUE)
  write.table(
    data, temporary, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t",
    na = "NA", digits = 17
  )
  if (!file.link(temporary, path)) {
    stop(sprintf("exclusive no-replace publication failed for engine output %s", name))
  }
  unlink(temporary)
}

if (!file.exists(task_path)) stop("locked task is absent")
if (!dir.exists(output_directory)) stop("output directory must be pre-created")
if (length(list.files(output_directory, all.files = TRUE, no.. = TRUE)) != 0L) {
  stop("output directory is not empty")
}
task_table <- read.delim(task_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(names(task_table), task_fields) || nrow(task_table) != 1L) {
  stop("locked Track B task schema/row count drifted")
}
task <- task_table[1L, , drop = FALSE]
p12_grid <- as.numeric(strsplit(task$p12_grid, ",", fixed = TRUE)[[1L]])
if (
  task$analysis_id != "track-b-v1.0-finemapping-trait-coloc" ||
  task$prior_method != "FLAT" || task$maximum_causal_signals != 10L ||
  task$credible_set_coverage != 0.95 || task$minimum_absolute_correlation != 0.5 ||
  task$maximum_iterations != 1000L || task$estimate_residual_variance != "FALSE" ||
  task$p1 != 1e-4 || task$p2 != 1e-4 ||
  !identical(p12_grid, c(1e-6, 5e-6, 1e-5, 5e-5)) ||
  task$single_signal_fallback != "FORBIDDEN" || task$maximum_locus_variants != "NONE" ||
  task$locus_splitting != "FORBIDDEN" || task$variant_thinning != "FORBIDDEN" ||
  task$lead_centered_truncation != "FORBIDDEN" ||
  !grepl("^[0-9a-f]{64}$", task$execution_amendment_sha256)
) {
  stop("Track B primary model/prior/no-subsetting contract drifted")
}

for (specification in list(
  c("summary1_path", "summary1_sha256"), c("summary2_path", "summary2_sha256"),
  c("variant_order_path", "variant_order_sha256"),
  c("signed_ld_path", "signed_ld_sha256"), c("ld_qc_path", "ld_qc_sha256"),
  c("ld_validation_path", "ld_validation_sha256"),
  c("ld_runtime_package_lock_path", "ld_runtime_package_lock_sha256"),
  c("ld_runtime_attestation_path", "ld_runtime_attestation_sha256"),
  c("runtime_package_lock_path", "runtime_package_lock_sha256")
)) {
  verify_file(task[[specification[[1L]]]], task[[specification[[2L]]]], specification[[1L]])
}

verify_runtime_lock("PRE", require_loaded_closure = FALSE)

summary1 <- read.delim(gzfile(task$summary1_path), stringsAsFactors = FALSE, check.names = FALSE)
summary2 <- read.delim(gzfile(task$summary2_path), stringsAsFactors = FALSE, check.names = FALSE)
order <- read.delim(task$variant_order_path, stringsAsFactors = FALSE, check.names = FALSE)
if (
  !identical(names(summary1), summary_fields) || !identical(names(summary2), summary_fields) ||
  !identical(names(order), order_fields) || nrow(order) != task$variant_count ||
  anyDuplicated(order$SNP) || nrow(summary1) != nrow(order) || nrow(summary2) != nrow(order)
) {
  stop("summary/order family is incomplete or schema-drifted")
}
identity_fields <- c("SNP", "CHR", "BP", "A1", "A2")
if (
  !identical(summary1[, identity_fields], order) ||
  !identical(summary2[, identity_fields], order) ||
  any(!is.finite(summary1$BETA)) || any(!is.finite(summary2$BETA)) ||
  any(!is.finite(summary1$SE)) || any(!is.finite(summary2$SE)) ||
  any(summary1$SE <= 0) || any(summary2$SE <= 0) ||
  any(!is.finite(summary1$P)) || any(!is.finite(summary2$P)) ||
  any(summary1$P < 0 | summary1$P > 1) || any(summary2$P < 0 | summary2$P > 1) ||
  any(!is.finite(summary1$EAF)) || any(!is.finite(summary2$EAF)) ||
  any(pmin(summary1$EAF, 1 - summary1$EAF) < 0.01) ||
  any(pmin(summary2$EAF, 1 - summary2$EAF) < 0.01) ||
  any(!is.finite(summary1$N)) || any(!is.finite(summary2$N)) ||
  any(summary1$N <= 0) || any(summary2$N <= 0)
) {
  stop("summary identity, allele, or finite-value validation failed")
}

make_sample_size_row <- function(summary, prefix, dataset_role, trait_id) {
  values <- as.numeric(summary$N)
  maximum <- max(values)
  observed <- c(
    count = length(values), min = min(values), q1 = quantile(values, 0.25, names = FALSE, type = 7),
    median = median(values), q3 = quantile(values, 0.75, names = FALSE, type = 7),
    max = maximum, max_to_min_ratio = maximum / min(values),
    fraction_below_90pct_max = mean(values < 0.90 * maximum),
    fraction_below_50pct_max = mean(values < 0.50 * maximum)
  )
  locked <- as.numeric(vapply(names(observed), function(suffix) {
    as.numeric(task[[paste0(prefix, "_per_snp_N_", suffix)]])
  }, numeric(1L)))
  tolerance <- pmax(1e-12, abs(locked) * 1e-12)
  dispersion <- if (all(values == values[[1L]])) {
    "EXACT_CONSTANT_PER_SNP_N"
  } else {
    "VARIABLE_PER_SNP_N_MEDIAN_SCALAR_APPROXIMATION_CLAIM_CAPPED"
  }
  scalar_rule <- if (dispersion == "EXACT_CONSTANT_PER_SNP_N") {
    "EXACT_CONSTANT_ELIGIBLE_PER_SNP_N"
  } else {
    "PRE_RESULT_MEDIAN_ELIGIBLE_PER_SNP_N_CLAIM_CAPPED"
  }
  if (
    any(!is.finite(locked)) || any(abs(observed - locked) > tolerance) ||
    abs(as.numeric(task[[paste0(prefix, "_scalar_N")]]) - observed[["median"]]) >
      max(1e-12, abs(observed[["median"]]) * 1e-12) ||
    as.character(task[[paste0(prefix, "_scalar_N_rule")]]) != scalar_rule ||
    as.character(task[[paste0(prefix, "_N_dispersion_status")]]) != dispersion ||
    !grepl("^[0-9a-f]{64}$", task[[paste0(prefix, "_per_snp_N_sha256")]])
  ) stop("full-vector sample-size diagnostics/scalar convention drifted")
  data.frame(
    analysis_id = task$analysis_id, pair_id = task$pair_id,
    family_role = task$family_role, locus_entry_id = task$locus_entry_id,
    dataset_role = dataset_role, trait_id = trait_id,
    per_snp_N_sha256 = task[[paste0(prefix, "_per_snp_N_sha256")]],
    per_snp_N_count = task[[paste0(prefix, "_per_snp_N_count")]],
    per_snp_N_min = task[[paste0(prefix, "_per_snp_N_min")]],
    per_snp_N_q1 = task[[paste0(prefix, "_per_snp_N_q1")]],
    per_snp_N_median = task[[paste0(prefix, "_per_snp_N_median")]],
    per_snp_N_q3 = task[[paste0(prefix, "_per_snp_N_q3")]],
    per_snp_N_max = task[[paste0(prefix, "_per_snp_N_max")]],
    per_snp_N_max_to_min_ratio = task[[paste0(prefix, "_per_snp_N_max_to_min_ratio")]],
    per_snp_N_fraction_below_90pct_max = task[[paste0(prefix, "_per_snp_N_fraction_below_90pct_max")]],
    per_snp_N_fraction_below_50pct_max = task[[paste0(prefix, "_per_snp_N_fraction_below_50pct_max")]],
    scalar_N = task[[paste0(prefix, "_scalar_N")]],
    scalar_N_rule = task[[paste0(prefix, "_scalar_N_rule")]],
    dispersion_status = dispersion,
    claim_cap = dispersion != "EXACT_CONSTANT_PER_SNP_N",
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

sample_size_diagnostics <- rbind(
  make_sample_size_row(summary1, "trait1", "TRAIT1_SLEEP", task$trait1),
  make_sample_size_row(summary2, "trait2", "TRAIT2_EXTERNAL", task$trait2)
)

ld_raw <- read.delim(gzfile(task$signed_ld_path), stringsAsFactors = FALSE,
                     check.names = FALSE, row.names = 1L)
R <- as.matrix(ld_raw)
storage.mode(R) <- "double"
if (
  nrow(R) != nrow(order) || ncol(R) != nrow(order) ||
  !identical(rownames(R), as.character(order$SNP)) ||
  !identical(colnames(R), as.character(order$SNP)) || any(!is.finite(R))
) {
  stop("signed-LD dimensions/names/order/finite QC failed")
}
ld_symmetry <- max(abs(R - t(R)))
ld_diagonal <- max(abs(diag(R) - 1))
ld_minimum_eigenvalue <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
if (
  ld_symmetry > 1e-8 || ld_diagonal > 1e-6 ||
  !is.finite(ld_minimum_eigenvalue) || ld_minimum_eigenvalue < -1e-6
) {
  stop("signed-LD symmetry/diagonal/PSD QC failed")
}

suppressPackageStartupMessages(library(susieR, lib.loc = ".r-env/lib/R/library"))
suppressPackageStartupMessages(library(coloc, lib.loc = ".r-env/lib/R/library"))
if (
  as.character(packageVersion("susieR")) != "0.14.2" ||
  as.character(packageVersion("coloc")) != "5.2.3"
) {
  stop("runtime package versions differ from susieR 0.14.2 / coloc 5.2.3")
}
if (!identical(sort(loadedNamespaces()), sort(runtime_namespace_closure))) {
  stop("actually loaded namespace closure differs from the runtime lock")
}

make_dataset <- function(summary, prefix) {
  dataset <- list(
    beta = as.numeric(summary$BETA), varbeta = as.numeric(summary$SE)^2,
    snp = as.character(summary$SNP), position = as.integer(summary$BP),
    MAF = pmin(as.numeric(summary$EAF), 1 - as.numeric(summary$EAF)),
    type = as.character(task[[paste0(prefix, "_type")]]),
    N = as.numeric(task[[paste0(prefix, "_scalar_N")]]), LD = R
  )
  # Exact case fractions are retained from frozen source counts for metadata;
  # runsusie uses beta/varbeta/z and scalar N. Quantitative sdY is deliberately
  # absent rather than estimated or invented.
  if (dataset$type == "cc") {
    dataset$s <- as.numeric(task[[paste0(prefix, "_case_fraction")]])
  } else if (dataset$type != "quant") {
    stop("dataset type is neither cc nor quant")
  }
  dataset
}

diagnose_dataset <- function(dataset) {
  tryCatch({
    z <- dataset$beta / sqrt(dataset$varbeta)
    s <- susieR::estimate_s_rss(z, dataset$LD, n = dataset$N)
    if (!is.finite(s) || s < 0 || s > 1) stop("estimate_s_rss returned invalid s")
    kriging <- suppressWarnings(
      susieR::kriging_rss(z, dataset$LD, n = dataset$N, s = s)
    )
    conditional <- kriging$conditional_dist
    if (
      is.null(conditional) || !all(c("logLR", "z") %in% names(conditional)) ||
      nrow(conditional) != length(z) || any(!is.finite(conditional$logLR)) ||
      any(!is.finite(conditional$z))
    ) stop("kriging_rss returned incomplete diagnostics")
    outliers <- which(conditional$logLR > 2 & abs(conditional$z) > 2)
    list(ok = TRUE, s = as.numeric(s), outliers = outliers, error = "")
  }, error = function(error) {
    list(ok = FALSE, s = NA_real_, outliers = integer(0), error = conditionMessage(error))
  })
}

fit_dataset <- function(dataset, suffix) {
  tryCatch({
    prior <- rep(1 / length(dataset$snp), length(dataset$snp))
    fit <- coloc::runsusie(
      dataset, suffix = suffix, maxit = 1000L, repeat_until_convergence = FALSE,
      L = 10L, coverage = 0.95, min_abs_corr = 0.5,
      estimate_residual_variance = FALSE, prior_weights = prior
    )
    if (!isTRUE(fit$converged) || is.null(fit$niter) || fit$niter > 1000L) {
      stop("SuSiE-RSS did not converge within exactly 1000 iterations")
    }
    if (length(fit$pip) != length(dataset$snp) || any(!is.finite(fit$pip)) ||
        any(fit$pip < 0 | fit$pip > 1)) stop("SuSiE-RSS returned invalid PIPs")
    list(ok = TRUE, fit = fit, error = "")
  }, error = function(error) {
    list(ok = FALSE, fit = NULL, error = conditionMessage(error))
  })
}

dataset1 <- make_dataset(summary1, "trait1")
dataset2 <- make_dataset(summary2, "trait2")
diagnostic1 <- diagnose_dataset(dataset1)
diagnostic2 <- diagnose_dataset(dataset2)
fit1 <- fit_dataset(dataset1, "trait1")
fit2 <- fit_dataset(dataset2, "trait2")

diagnostic_status <- function(diagnostic, attempt) {
  if (!diagnostic$ok || is.na(diagnostic$s)) return("FAILED_REQUIRED_DIAGNOSTIC")
  if (length(diagnostic$outliers)) return("FAILED_KRIGING_QC")
  if (!attempt$ok || is.null(attempt$fit) || !isTRUE(attempt$fit$converged)) {
    return("FAILED_MODEL_CONVERGENCE")
  }
  "PASS"
}
status1 <- diagnostic_status(diagnostic1, fit1)
status2 <- diagnostic_status(diagnostic2, fit2)

fine_claim <- paste(
  "PIPs and credible sets are model-based candidate-variant evidence and do not prove",
  "a causal variant or biological mechanism"
)
coloc_claim <- paste(
  "Trait colocalization supports a shared-signal model under stated data, LD, and priors;",
  "it does not prove causality, mediation, a causal gene, or mechanism"
)

fit_cs <- function(attempt) {
  if (!attempt$ok || is.null(attempt$fit) || is.null(attempt$fit$sets$cs)) list() else attempt$fit$sets$cs
}

make_diagnostic_row <- function(summary, prefix, trait_id, trait_role, diagnostic, attempt, status) {
  fit <- attempt$fit
  cs <- fit_cs(attempt)
  error <- paste(c(diagnostic$error, attempt$error)[nzchar(c(diagnostic$error, attempt$error))], collapse = ";")
  # susieR documents larger s as stronger z/LD inconsistency but does not
  # define 0.10 (or any value) as a validity cutoff. Keep this predeclared
  # warning conspicuous without converting it to INVALID_INPUT.
  if (diagnostic$ok && diagnostic$s > 0.10) {
    error <- paste(c(error, "WARNING_RSS_LD_S_GT_0.10_NONFATAL")[nzchar(c(error, "WARNING_RSS_LD_S_GT_0.10_NONFATAL"))], collapse = ";")
  }
  if (task[[paste0(prefix, "_N_dispersion_status")]] != "EXACT_CONSTANT_PER_SNP_N") {
    warning <- "WARNING_VARIABLE_PER_SNP_N_MEDIAN_SCALAR_APPROXIMATION_CLAIM_CAPPED"
    error <- paste(c(error, warning)[nzchar(c(error, warning))], collapse = ";")
  }
  if (length(diagnostic$outliers)) {
    error <- paste(c(error, "KRIGING_ALLELE_SWITCH_OUTLIERS_PRESENT")[nzchar(c(error, "KRIGING_ALLELE_SWITCH_OUTLIERS_PRESENT"))], collapse = ";")
  }
  data.frame(
    analysis_id = task$analysis_id, pair_id = task$pair_id,
    family_role = task$family_role, locus_entry_id = task$locus_entry_id,
    dataset_role = trait_role, trait_id = trait_id, variant_count = nrow(summary),
    ld_variant_count = nrow(R), allele_match_count = nrow(summary),
    model_converged = !is.null(fit) && isTRUE(fit$converged),
    niter = if (is.null(fit)) NA else fit$niter, credible_set_count = length(cs),
    max_pip = if (is.null(fit)) NA else max(fit$pip), rss_ld_s = diagnostic$s,
    kriging_allele_switch_outlier_count = length(diagnostic$outliers),
    kriging_allele_switch_outliers = if (length(diagnostic$outliers)) {
      paste(summary$SNP[diagnostic$outliers], collapse = ";")
    } else "NONE",
    ld_symmetry_max_abs = ld_symmetry, ld_diagonal_max_abs = ld_diagonal,
    ld_minimum_eigenvalue = ld_minimum_eigenvalue,
    reference_sample_size = task$reference_sample_size,
    diagnostic_status = status, error = if (nzchar(error)) error else "NA",
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

diagnostics <- rbind(
  make_diagnostic_row(summary1, "trait1", task$trait1, "TRAIT1_SLEEP", diagnostic1, fit1, status1),
  make_diagnostic_row(summary2, "trait2", task$trait2, "TRAIT2_EXTERNAL", diagnostic2, fit2, status2)
)

make_qc_row <- function(row, trait_role) {
  data.frame(
    analysis_id = row$analysis_id, pair_id = row$pair_id, family_role = row$family_role,
    locus_entry_id = row$locus_entry_id, trait_id = row$trait_id,
    trait_role = trait_role, variant_count = row$variant_count,
    model_converged = row$model_converged, niter = row$niter,
    credible_set_count = row$credible_set_count, max_pip = row$max_pip,
    rss_ld_s = row$rss_ld_s,
    kriging_allele_switch_outlier_count = row$kriging_allele_switch_outlier_count,
    ld_symmetry_max_abs = row$ld_symmetry_max_abs,
    ld_diagonal_max_abs = row$ld_diagonal_max_abs,
    ld_minimum_eigenvalue = row$ld_minimum_eigenvalue,
    reference_sample_size = row$reference_sample_size,
    analysis_status = if (row$diagnostic_status == "PASS") "COMPLETE" else "FAILED_QC",
    diagnostic_status = row$diagnostic_status, error = row$error,
    stringsAsFactors = FALSE, check.names = FALSE
  )
}
finemap_qc <- rbind(
  make_qc_row(diagnostics[1L, , drop = FALSE], "TRAIT1_SLEEP"),
  make_qc_row(diagnostics[2L, , drop = FALSE], "TRAIT2_EXTERNAL")
)

make_credible_rows <- function(summary, trait_id, trait_role, attempt, status) {
  fit <- attempt$fit
  cs <- fit_cs(attempt)
  if (!length(cs)) return(empty_frame(credible_fields))
  output <- vector("list", length(cs))
  for (index in seq_along(cs)) {
    component <- fit$sets$cs_index[[index]]
    members <- as.integer(cs[[index]])
    purity <- fit$sets$purity[index, ]
    lead_index <- which.max(fit$alpha[component, ])
    output[[index]] <- data.frame(
      analysis_id = task$analysis_id, pair_id = task$pair_id,
      family_role = task$family_role, locus_entry_id = task$locus_entry_id,
      trait_id = trait_id, trait_role = trait_role, signal_id = paste0("L", component),
      component_index = component, lead_snp = summary$SNP[[lead_index]],
      lead_pip = fit$pip[[lead_index]], credible_set_size = length(members),
      credible_set_snps = paste(summary$SNP[members], collapse = ";"),
      requested_coverage = fit$sets$requested_coverage,
      achieved_coverage = fit$sets$coverage[[index]], min_abs_corr = purity$min.abs.corr,
      mean_abs_corr = purity$mean.abs.corr, median_abs_corr = purity$median.abs.corr,
      cs_log10bf = if (!is.null(fit$lbf)) fit$lbf[[component]] / log(10) else NA_real_,
      model_converged = fit$converged,
      analysis_status = if (status == "PASS") "COMPLETE" else "FAILED_QC",
      claim_limit = fine_claim, stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  do.call(rbind, output)
}

credible1 <- make_credible_rows(summary1, task$trait1, "TRAIT1_SLEEP", fit1, status1)
credible2 <- make_credible_rows(summary2, task$trait2, "TRAIT2_EXTERNAL", fit2, status2)
credible_sets <- rbind(credible1, credible2)

make_variant_rows <- function(summary, trait_id, trait_role, diagnostic, attempt, status, retained_error) {
  fit <- attempt$fit
  count <- nrow(summary)
  memberships <- vector("list", count)
  sizes <- vector("list", count)
  requested <- vector("list", count)
  achieved <- vector("list", count)
  minimum <- vector("list", count)
  means <- vector("list", count)
  if (!is.null(fit) && length(fit_cs(attempt))) {
    for (index in seq_along(fit$sets$cs)) {
      component <- fit$sets$cs_index[[index]]
      members <- as.integer(fit$sets$cs[[index]])
      purity <- fit$sets$purity[index, ]
      for (variant in members) {
        memberships[[variant]] <- c(memberships[[variant]], paste0("L", component))
        sizes[[variant]] <- c(sizes[[variant]], length(members))
        requested[[variant]] <- c(requested[[variant]], fit$sets$requested_coverage)
        achieved[[variant]] <- c(achieved[[variant]], fit$sets$coverage[[index]])
        minimum[[variant]] <- c(minimum[[variant]], purity$min.abs.corr)
        means[[variant]] <- c(means[[variant]], purity$mean.abs.corr)
      }
    }
  }
  collapse_or <- function(values, missing = "NA") {
    vapply(values, function(value) if (length(value)) paste(value, collapse = ";") else missing, character(1L))
  }
  max_component <- if (is.null(fit)) rep("NA", count) else paste0("L", apply(fit$alpha, 2L, which.max))
  data.frame(
    analysis_id = task$analysis_id, pair_id = task$pair_id,
    family_role = task$family_role, locus_entry_id = task$locus_entry_id,
    locus_id = task$ld_block_id, CHR = task$CHR, START = task$START, STOP = task$STOP,
    trait_id = trait_id, trait_role = trait_role, SNP = summary$SNP, BP = summary$BP,
    A1 = summary$A1, A2 = summary$A2, BETA = summary$BETA, SE = summary$SE,
    P = summary$P, EAF = summary$EAF, INFO = summary$INFO, N = summary$N,
    prior_role = "PRIMARY_FLAT", prior_weight = rep(1 / count, count),
    PIP = if (is.null(fit)) rep(NA_real_, count) else fit$pip,
    credible_set_ids = collapse_or(memberships, "NONE"),
    credible_set_sizes = collapse_or(sizes),
    credible_set_requested_coverage = collapse_or(requested),
    credible_set_achieved_coverage = collapse_or(achieved),
    credible_set_min_abs_corr = collapse_or(minimum),
    credible_set_mean_abs_corr = collapse_or(means),
    max_alpha_component = max_component,
    model_converged = !is.null(fit) && isTRUE(fit$converged),
    niter = if (is.null(fit)) NA else fit$niter, rss_ld_s = diagnostic$s,
    kriging_allele_switch_outlier = seq_len(count) %in% diagnostic$outliers,
    diagnostic_status = status, error = retained_error,
    claim_limit = fine_claim, stringsAsFactors = FALSE, check.names = FALSE
  )
}

variants1 <- make_variant_rows(
  summary1, task$trait1, "TRAIT1_SLEEP", diagnostic1, fit1, status1,
  as.character(diagnostics$error[[1L]])
)
variants2 <- make_variant_rows(
  summary2, task$trait2, "TRAIT2_EXTERNAL", diagnostic2, fit2, status2,
  as.character(diagnostics$error[[2L]])
)

add_coloc_row <- function(p12, signal1 = "NA", signal2 = "NA", nsnps = nrow(order),
                          h0 = NA_real_, h1 = NA_real_, h2 = NA_real_, h3 = NA_real_,
                          h4 = NA_real_, top_snp = "NA", top_pp = NA_real_,
                          engine_status = "COLOC_SUSIE_NOT_ESTIMABLE", error = "NA") {
  ratio <- if (is.finite(h3) && h3 > 0) h4 / h3 else if (is.finite(h4) && h4 > 0) Inf else NA_real_
  data.frame(
    analysis_id = task$analysis_id, pair_id = task$pair_id,
    family_role = task$family_role, locus_entry_id = task$locus_entry_id,
    locus_id = task$ld_block_id, CHR = task$CHR, START = task$START, STOP = task$STOP,
    trait1 = task$trait1, trait2 = task$trait2, signal1 = signal1, signal2 = signal2,
    coloc_method = "coloc.susie", p1 = 1e-4, p2 = 1e-4, p12 = p12,
    prior_role = if (p12 == 1e-5) "PRIMARY" else "SENSITIVITY", nsnps = nsnps,
    PP_H0 = h0, PP_H1 = h1, PP_H2 = h2, PP_H3 = h3, PP_H4 = h4,
    PP_H4_over_PP_H3 = ratio, top_shared_variant = top_snp,
    top_shared_variant_PP_H4 = top_pp,
    fine_mapping_qc = if (status1 == "PASS" && status2 == "PASS") "PASS" else "FAILED_QC",
    ld_qc = if (status1 == "PASS" && status2 == "PASS") "PASS" else "FAILED_QC",
    engine_status = engine_status, error = error, claim_limit = coloc_claim,
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

coloc_rows <- list()
coloc_index <- 0L
any_coloc_failure <- FALSE
any_no_estimable <- FALSE
fits_valid <- fit1$ok && fit2$ok && status1 == "PASS" && status2 == "PASS"
has_cs <- fits_valid && length(fit_cs(fit1)) > 0L && length(fit_cs(fit2)) > 0L
if (has_cs) {
  for (p12 in p12_grid) {
    attempt <- tryCatch(
      list(ok = TRUE, value = coloc::coloc.susie(fit1$fit, fit2$fit, p1 = 1e-4, p2 = 1e-4, p12 = p12), error = ""),
      error = function(error) list(ok = FALSE, value = NULL, error = conditionMessage(error))
    )
    if (!attempt$ok) {
      any_coloc_failure <- TRUE
      coloc_index <- coloc_index + 1L
      coloc_rows[[coloc_index]] <- add_coloc_row(
        p12, engine_status = "COLOC_SUSIE_FAILURE", error = attempt$error
      )
      next
    }
    result <- attempt$value
    summaries <- as.data.frame(result$summary)
    variants <- as.data.frame(result$results)
    required_pp <- c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
    if (!nrow(summaries) || !all(required_pp %in% names(summaries))) {
      any_no_estimable <- TRUE
      coloc_index <- coloc_index + 1L
      coloc_rows[[coloc_index]] <- add_coloc_row(
        p12, engine_status = "COLOC_SUSIE_NO_ESTIMABLE_SIGNAL_PAIR"
      )
      next
    }
    for (index in seq_len(nrow(summaries))) {
      pp_column <- if (nrow(summaries) == 1L) "SNP.PP.H4.abf" else paste0("SNP.PP.H4.row", index)
      top_index <- if (pp_column %in% names(variants)) which.max(variants[[pp_column]]) else integer(0)
      coloc_index <- coloc_index + 1L
      coloc_rows[[coloc_index]] <- add_coloc_row(
        p12, signal1 = paste0("L", summaries$idx1[[index]]),
        signal2 = paste0("L", summaries$idx2[[index]]), nsnps = summaries$nsnps[[index]],
        h0 = summaries$PP.H0.abf[[index]], h1 = summaries$PP.H1.abf[[index]],
        h2 = summaries$PP.H2.abf[[index]], h3 = summaries$PP.H3.abf[[index]],
        h4 = summaries$PP.H4.abf[[index]],
        top_snp = if (length(top_index)) variants$snp[[top_index]] else "NA",
        top_pp = if (length(top_index)) variants[[pp_column]][[top_index]] else NA_real_,
        engine_status = "COLOC_SUSIE_COMPLETE"
      )
    }
  }
} else {
  reason <- if (fits_valid) "NO_VALID_CREDIBLE_SET_PAIR" else "INVALID_FINE_MAPPING_DIAGNOSTICS"
  for (p12 in p12_grid) {
    coloc_index <- coloc_index + 1L
    coloc_rows[[coloc_index]] <- add_coloc_row(
      p12, engine_status = reason,
      error = if (fits_valid) "NA" else paste(c(diagnostics$error)[diagnostics$error != "NA"], collapse = ";")
    )
  }
}
if (has_cs && (any_coloc_failure || any_no_estimable)) {
  retained_error <- if (any_coloc_failure) {
    paste(vapply(coloc_rows, function(row) as.character(row$error[[1L]]), character(1L)), collapse = ";")
  } else "NA"
  retained_status <- if (any_coloc_failure) "COLOC_SUSIE_FAILURE" else "COLOC_SUSIE_NO_ESTIMABLE_SIGNAL_PAIR"
  coloc_rows <- lapply(p12_grid, function(p12) {
    add_coloc_row(p12, engine_status = retained_status, error = retained_error)
  })
}
coloc_raw <- do.call(rbind, coloc_rows)

engine_status <- data.frame(
  analysis_id = task$analysis_id, pair_id = task$pair_id,
  locus_entry_id = task$locus_entry_id,
  engine_status = "COMPLETE_SERIALIZED",
  trait1_status = if (status1 == "PASS") "COMPLETE" else "FAILED_QC",
  trait2_status = if (status2 == "PASS") "COMPLETE" else "FAILED_QC",
  coloc_status = if (fits_valid) "SERIALIZED" else "INVALID_INPUT_SERIALIZED",
  error = if (status1 == "PASS" && status2 == "PASS") "NA" else paste(
    diagnostics$error[diagnostics$error != "NA"], collapse = ";"
  ),
  stringsAsFactors = FALSE, check.names = FALSE
)

write_tsv_no_replace(variants1[, raw_variant_fields], "trait1.raw.tsv", raw_variant_fields)
write_tsv_no_replace(variants2[, raw_variant_fields], "trait2.raw.tsv", raw_variant_fields)
write_tsv_no_replace(credible_sets[, credible_fields], "credible_sets.tsv", credible_fields)
write_tsv_no_replace(diagnostics[, diagnostic_fields], "diagnostics.tsv", diagnostic_fields)
write_tsv_no_replace(
  sample_size_diagnostics[, sample_size_fields], "sample_size_diagnostics.tsv", sample_size_fields
)
write_tsv_no_replace(finemap_qc[, finemap_qc_fields], "finemap_qc.tsv", finemap_qc_fields)
write_tsv_no_replace(coloc_raw[, raw_coloc_fields], "coloc.raw.tsv", raw_coloc_fields)
write_tsv_no_replace(engine_status[, engine_status_fields], "engine_status.tsv", engine_status_fields)

# Rehash the exact installed namespace closure after every scientific output
# has been constructed and written. The coordinator performs a third check
# after process exit before it can seal the run.
verify_runtime_lock("POST", require_loaded_closure = TRUE)
runtime_attestation <- data.frame(
  analysis_id = task$analysis_id, locus_entry_id = task$locus_entry_id,
  execution_amendment_sha256 = task$execution_amendment_sha256,
  runtime_package_lock_sha256 = task$runtime_package_lock_sha256,
  precheck_status = "PASS_EXACT_PATHS_AND_BYTES",
  postcheck_status = "PASS_EXACT_LOADED_CLOSURE_PATHS_AND_BYTES",
  loaded_namespace_count = length(loadedNamespaces()),
  stringsAsFactors = FALSE, check.names = FALSE
)
write_tsv_no_replace(
  runtime_attestation[, runtime_attestation_fields],
  "runtime_attestation.tsv", runtime_attestation_fields
)

cat(sprintf(
  paste0(
    "TRACK_B_SUSIE_COLOC_SERIALIZED locus=%s variants=%d trait1=%s trait2=%s ",
    "coloc_rows=%d\n"
  ), task$locus_entry_id, nrow(order), status1, status2, nrow(coloc_raw)
))
rm(R, ld_raw, dataset1, dataset2, fit1, fit2)
invisible(gc(full = TRUE))
