#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 14) {
  stop(paste(
    "usage: 35_run_susie_coloc.R MANIFEST COMPARISON_ID VARIANT_OUT CREDIBLE_SET_OUT",
    "COLOC_OUT SHARED_VARIANT_OUT DIAGNOSTIC_OUT L COVERAGE MIN_ABS_CORR MAXIT P1 P2 P12_CSV"
  ))
}

manifest_path <- args[[1]]
comparison_id <- args[[2]]
variant_out <- args[[3]]
credible_set_out <- args[[4]]
coloc_out <- args[[5]]
shared_variant_out <- args[[6]]
diagnostic_out <- args[[7]]
L <- as.integer(args[[8]])
coverage <- as.numeric(args[[9]])
min_abs_corr <- as.numeric(args[[10]])
maxit <- as.integer(args[[11]])
p1 <- as.numeric(args[[12]])
p2 <- as.numeric(args[[13]])
p12_grid <- as.numeric(strsplit(args[[14]], ",", fixed = TRUE)[[1]])

if (!is.finite(L) || L < 1 || !is.finite(coverage) || coverage <= 0 || coverage >= 1 ||
    !is.finite(min_abs_corr) || min_abs_corr < 0 || min_abs_corr > 1 ||
    !is.finite(maxit) || maxit < 1 || !is.finite(p1) || !is.finite(p2) ||
    any(!is.finite(p12_grid)) || any(p12_grid <= 0) || length(p12_grid) < 1) {
  stop("invalid SuSiE/coloc parameters")
}

suppressPackageStartupMessages(library(susieR))
suppressPackageStartupMessages(library(coloc))
if (as.character(packageVersion("susieR")) != "0.14.2" || as.character(packageVersion("coloc")) != "5.2.3") {
  stop("runtime package version differs from locked workflow")
}

sha256_file <- function(path) {
  sha256sum <- Sys.which("sha256sum")
  shasum <- Sys.which("shasum")
  if (nzchar(sha256sum)) {
    output <- system2(sha256sum, path, stdout = TRUE, stderr = TRUE)
  } else if (nzchar(shasum)) {
    output <- system2(shasum, c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  } else {
    stop("neither sha256sum nor shasum is available")
  }
  if (!length(output)) stop("SHA-256 command produced no output")
  strsplit(trimws(output[[1]]), "[[:space:]]+")[[1]][[1]]
}

verify_file <- function(path, expected, label) {
  if (!file.exists(path) || nchar(expected) != 64 || sha256_file(path) != expected) {
    stop(sprintf("locked %s file is missing or checksum-mismatched", label))
  }
}

write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(data, path, quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")
}

variant_fields <- c(
  "comparison_id", "pair_id", "locus_id", "comparison_type", "dataset_id", "dataset_role",
  "dataset_type", "SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "MAF", "INFO",
  "prior_method", "normalized_prior_weight", "PIP", "credible_set_ids", "max_alpha_component",
  "model_converged", "rss_ld_s", "kriging_allele_switch_outlier"
)
credible_fields <- c(
  "comparison_id", "pair_id", "locus_id", "comparison_type", "dataset_id", "dataset_role",
  "signal_id", "component_index", "lead_snp", "lead_pip", "credible_set_size",
  "credible_set_snps", "requested_coverage", "achieved_coverage", "min_abs_corr",
  "mean_abs_corr", "median_abs_corr", "cs_log10bf", "model_converged"
)
coloc_fields <- c(
  "comparison_id", "pair_id", "locus_id", "comparison_type", "dataset1_id", "dataset2_id",
  "molecular_feature_id", "tissue_cell_context", "coloc_method", "p1", "p2", "p12",
  "prior_role", "signal1", "signal2", "hit1", "hit2", "nsnps", "PP_H0", "PP_H1",
  "PP_H2", "PP_H3", "PP_H4", "PP_H4_over_PP_H3", "top_shared_variant",
  "top_shared_variant_PP_H4", "fine_mapping_qc", "analysis_status",
  "single_signal_fallback_justification", "claim_limit"
)
shared_fields <- c(
  "comparison_id", "pair_id", "locus_id", "comparison_type", "coloc_method", "p12",
  "signal1", "signal2", "SNP", "SNP_PP_H4"
)
diagnostic_fields <- c(
  "comparison_id", "pair_id", "locus_id", "comparison_type", "dataset_id", "dataset_role",
  "variant_count", "model_converged", "niter", "credible_set_count", "max_pip", "rss_ld_s",
  "kriging_allele_switch_outlier_count", "kriging_allele_switch_outliers", "diagnostic_status"
)

empty_frame <- function(fields) {
  as.data.frame(setNames(replicate(length(fields), character(0), simplify = FALSE), fields), stringsAsFactors = FALSE)
}

manifest <- read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
selected <- manifest[manifest$queue_row_id == comparison_id, , drop = FALSE]
if (nrow(selected) != 1) stop("COMPARISON_ID must identify exactly one locked manifest row")
row <- selected[1, ]
if (row$source_search_status != "VERIFIED_ANALYSIS" || row$results_accessed_before_lock != "NO") {
  stop("comparison is not a result-free verified analysis row")
}
for (spec in list(
  c("summary1_path", "summary1_sha256"), c("summary2_path", "summary2_sha256"),
  c("ld_path", "ld_sha256"), c("ld_variant_order_path", "ld_variant_order_sha256")
)) {
  verify_file(row[[spec[[1]]]], row[[spec[[2]]]], spec[[1]])
}
for (prefix in c("dataset1", "dataset2")) {
  if (row[[paste0(prefix, "_prior_method")]] == "POLYFUN_FUNCTIONAL") {
    verify_file(row[[paste0(prefix, "_prior_source_path")]], row[[paste0(prefix, "_prior_source_sha256")]], paste0(prefix, " prior source"))
  }
}

order <- read.delim(row$ld_variant_order_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(names(order), "SNP") || anyDuplicated(order$SNP)) stop("invalid locked variant-order file")
summary1 <- read.delim(row$summary1_path, stringsAsFactors = FALSE, check.names = FALSE)
summary2 <- read.delim(row$summary2_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "MAF", "INFO", "PRIOR_WEIGHT")
if (!all(required %in% names(summary1)) || !all(required %in% names(summary2))) stop("summary input lacks canonical columns")
identity_fields <- c("SNP", "CHR", "BP", "A1", "A2")
if (!identical(summary1[, identity_fields], summary2[, identity_fields]) ||
    !identical(as.character(summary1$SNP), as.character(order$SNP))) {
  stop("summary inputs and LD order are not coordinate/allele identical")
}

ld_raw <- read.delim(row$ld_path, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
R <- as.matrix(ld_raw)
storage.mode(R) <- "double"
if (!identical(rownames(R), as.character(order$SNP)) || !identical(colnames(R), as.character(order$SNP))) {
  stop("signed LD matrix rows/columns differ from locked order")
}
if (any(!is.finite(R)) || max(abs(R - t(R))) > 1e-8 || max(abs(diag(R) - 1)) > 1e-6) {
  stop("signed LD matrix fails finite/symmetry/unit-diagonal QC")
}
minimum_eigenvalue <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
if (!is.finite(minimum_eigenvalue) || minimum_eigenvalue < -1e-6) stop("signed LD matrix is not positive semidefinite")

make_dataset <- function(dat, prefix) {
  dataset <- list(
    beta = as.numeric(dat$BETA), varbeta = as.numeric(dat$SE)^2,
    snp = as.character(dat$SNP), position = as.integer(dat$BP),
    MAF = as.numeric(dat$MAF), type = row[[paste0(prefix, "_type")]],
    N = as.numeric(row[[paste0(prefix, "_N")]]), LD = R
  )
  if (dataset$type == "quant") dataset$sdY <- as.numeric(row[[paste0(prefix, "_sdY")]])
  if (dataset$type == "cc") dataset$s <- as.numeric(row[[paste0(prefix, "_case_fraction")]])
  dataset
}

normalized_prior <- function(dat, prefix) {
  if (row[[paste0(prefix, "_prior_method")]] == "FLAT") return(rep(1 / nrow(dat), nrow(dat)))
  weights <- as.numeric(dat$PRIOR_WEIGHT)
  if (any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) stop("invalid PolyFun prior weights")
  weights / sum(weights)
}

fit_dataset <- function(dataset, dat, prefix) {
  fit_args <- list(
    d = dataset, suffix = prefix, maxit = maxit, repeat_until_convergence = FALSE,
    L = min(L, nrow(dat)), coverage = coverage, min_abs_corr = min_abs_corr,
    estimate_residual_variance = FALSE
  )
  if (row[[paste0(prefix, "_prior_method")]] == "POLYFUN_FUNCTIONAL") {
    fit_args$prior_weights <- normalized_prior(dat, prefix)
  }
  do.call(coloc::runsusie, fit_args)
}

diagnose_dataset <- function(dataset) {
  z <- dataset$beta / sqrt(dataset$varbeta)
  s <- susieR::estimate_s_rss(z, dataset$LD, n = dataset$N)
  kriging <- suppressWarnings(susieR::kriging_rss(z, dataset$LD, n = dataset$N, s = s))
  outlier <- which(kriging$conditional_dist$logLR > 2 & abs(kriging$conditional_dist$z) > 2)
  list(s = as.numeric(s), outlier = outlier)
}

dataset1 <- make_dataset(summary1, "dataset1")
dataset2 <- make_dataset(summary2, "dataset2")
diag1 <- diagnose_dataset(dataset1)
diag2 <- diagnose_dataset(dataset2)
fit1_attempt <- tryCatch(list(ok = TRUE, fit = fit_dataset(dataset1, summary1, "dataset1"), error = ""),
                         error = function(error) list(ok = FALSE, fit = NULL, error = conditionMessage(error)))
fit2_attempt <- tryCatch(list(ok = TRUE, fit = fit_dataset(dataset2, summary2, "dataset2"), error = ""),
                         error = function(error) list(ok = FALSE, fit = NULL, error = conditionMessage(error)))

diagnostic_row <- function(prefix, dat, diagnostic, attempt) {
  fit <- attempt$fit
  cs_count <- if (is.null(fit) || is.null(fit$sets$cs)) 0 else length(fit$sets$cs)
  max_pip <- if (is.null(fit)) NA_real_ else max(fit$pip)
  converged <- !is.null(fit) && isTRUE(fit$converged)
  status <- if (!attempt$ok || !converged) "FINE_MAPPING_FAILED_OR_NONCONVERGED" else if (length(diagnostic$outlier)) "LD_SUMSTAT_INCONSISTENCY_FLAGGED" else "PASS"
  data.frame(
    comparison_id = comparison_id, pair_id = row$pair_id, locus_id = row$locus_id,
    comparison_type = row$comparison_type, dataset_id = row[[paste0(prefix, "_id")]],
    dataset_role = row[[paste0(prefix, "_role")]], variant_count = nrow(dat),
    model_converged = converged, niter = if (is.null(fit)) NA else fit$niter,
    credible_set_count = cs_count, max_pip = max_pip, rss_ld_s = diagnostic$s,
    kriging_allele_switch_outlier_count = length(diagnostic$outlier),
    kriging_allele_switch_outliers = if (length(diagnostic$outlier)) paste(dat$SNP[diagnostic$outlier], collapse = ";") else "NONE",
    diagnostic_status = if (!attempt$ok) paste0(status, ":", attempt$error) else status,
    stringsAsFactors = FALSE
  )
}
diagnostics <- rbind(
  diagnostic_row("dataset1", summary1, diag1, fit1_attempt),
  diagnostic_row("dataset2", summary2, diag2, fit2_attempt)
)
write_tsv(diagnostics[, diagnostic_fields], diagnostic_out)

if (!fit1_attempt$ok || !fit2_attempt$ok || !isTRUE(fit1_attempt$fit$converged) || !isTRUE(fit2_attempt$fit$converged)) {
  write_tsv(empty_frame(variant_fields), variant_out)
  write_tsv(empty_frame(credible_fields), credible_set_out)
  write_tsv(empty_frame(coloc_fields), coloc_out)
  write_tsv(empty_frame(shared_fields), shared_variant_out)
  cat(sprintf("SUSIE_COLOC_COMPARISON_RETAINED_FAILURE comparison=%s\n", comparison_id))
  quit(save = "no", status = 0)
}
fit1 <- fit1_attempt$fit
fit2 <- fit2_attempt$fit

variant_rows <- function(prefix, dat, fit, diagnostic) {
  memberships <- vector("list", nrow(dat))
  if (!is.null(fit$sets$cs)) {
    for (index in seq_along(fit$sets$cs)) {
      component <- fit$sets$cs_index[[index]]
      for (variant in as.integer(fit$sets$cs[[index]])) memberships[[variant]] <- c(memberships[[variant]], paste0("L", component))
    }
  }
  alpha_component <- apply(fit$alpha, 2, which.max)
  data.frame(
    comparison_id = comparison_id, pair_id = row$pair_id, locus_id = row$locus_id,
    comparison_type = row$comparison_type, dataset_id = row[[paste0(prefix, "_id")]],
    dataset_role = row[[paste0(prefix, "_role")]], dataset_type = row[[paste0(prefix, "_type")]],
    SNP = dat$SNP, CHR = dat$CHR, BP = dat$BP, A1 = dat$A1, A2 = dat$A2,
    BETA = dat$BETA, SE = dat$SE, MAF = dat$MAF, INFO = dat$INFO,
    prior_method = row[[paste0(prefix, "_prior_method")]],
    normalized_prior_weight = normalized_prior(dat, prefix), PIP = fit$pip,
    credible_set_ids = vapply(memberships, function(value) if (length(value)) paste(value, collapse = ";") else "NONE", character(1)),
    max_alpha_component = paste0("L", alpha_component), model_converged = fit$converged,
    rss_ld_s = diagnostic$s, kriging_allele_switch_outlier = seq_len(nrow(dat)) %in% diagnostic$outlier,
    stringsAsFactors = FALSE
  )
}
variants <- rbind(
  variant_rows("dataset1", summary1, fit1, diag1),
  variant_rows("dataset2", summary2, fit2, diag2)
)
write_tsv(variants[, variant_fields], variant_out)

credible_rows <- function(prefix, dat, fit) {
  if (is.null(fit$sets$cs) || !length(fit$sets$cs)) return(empty_frame(credible_fields))
  output <- vector("list", length(fit$sets$cs))
  for (index in seq_along(fit$sets$cs)) {
    component <- fit$sets$cs_index[[index]]
    members <- as.integer(fit$sets$cs[[index]])
    lead_index <- which.max(fit$alpha[component, ])
    purity <- fit$sets$purity[index, ]
    output[[index]] <- data.frame(
      comparison_id = comparison_id, pair_id = row$pair_id, locus_id = row$locus_id,
      comparison_type = row$comparison_type, dataset_id = row[[paste0(prefix, "_id")]],
      dataset_role = row[[paste0(prefix, "_role")]], signal_id = paste0("L", component),
      component_index = component, lead_snp = dat$SNP[[lead_index]], lead_pip = fit$pip[[lead_index]],
      credible_set_size = length(members), credible_set_snps = paste(dat$SNP[members], collapse = ";"),
      requested_coverage = fit$sets$requested_coverage, achieved_coverage = fit$sets$coverage[[index]],
      min_abs_corr = purity$min.abs.corr, mean_abs_corr = purity$mean.abs.corr,
      median_abs_corr = purity$median.abs.corr,
      cs_log10bf = if (!is.null(fit$lbf)) fit$lbf[[component]] / log(10) else NA_real_,
      model_converged = fit$converged, stringsAsFactors = FALSE
    )
  }
  do.call(rbind, output)
}
credible_sets <- rbind(
  credible_rows("dataset1", summary1, fit1),
  credible_rows("dataset2", summary2, fit2)
)
write_tsv(credible_sets[, credible_fields], credible_set_out)

fine_mapping_qc <- if (length(diag1$outlier) || length(diag2$outlier)) "LD_SUMSTAT_INCONSISTENCY_FLAGGED" else "PASS"
has_cs <- !is.null(fit1$sets$cs) && length(fit1$sets$cs) && !is.null(fit2$sets$cs) && length(fit2$sets$cs)
coloc_rows <- list()
shared_rows <- list()
row_index <- 0
shared_index <- 0
claim_limit <- "SHARED_SIGNAL_MODEL_NOT_CAUSAL_VARIANT_GENE_PATHWAY_OR_MECHANISM"

add_coloc_row <- function(method, p12, prior_role, signal1, signal2, hit1, hit2, nsnps,
                          h0, h1, h2, h3, h4, top_snp, top_pp, status) {
  ratio <- if (is.finite(h3) && h3 > 0) h4 / h3 else if (is.finite(h4) && h4 > 0) Inf else NA_real_
  data.frame(
    comparison_id = comparison_id, pair_id = row$pair_id, locus_id = row$locus_id,
    comparison_type = row$comparison_type, dataset1_id = row$dataset1_id, dataset2_id = row$dataset2_id,
    molecular_feature_id = row$molecular_feature_id, tissue_cell_context = row$tissue_cell_context,
    coloc_method = method, p1 = p1, p2 = p2, p12 = p12, prior_role = prior_role,
    signal1 = signal1, signal2 = signal2, hit1 = hit1, hit2 = hit2, nsnps = nsnps,
    PP_H0 = h0, PP_H1 = h1, PP_H2 = h2, PP_H3 = h3, PP_H4 = h4,
    PP_H4_over_PP_H3 = ratio, top_shared_variant = top_snp, top_shared_variant_PP_H4 = top_pp,
    fine_mapping_qc = fine_mapping_qc, analysis_status = status,
    single_signal_fallback_justification = row$single_signal_fallback_justification,
    claim_limit = claim_limit, stringsAsFactors = FALSE
  )
}

if (has_cs) {
  for (p12 in p12_grid) {
    coloc_result <- coloc::coloc.susie(fit1, fit2, p1 = p1, p2 = p2, p12 = p12)
    summary <- as.data.frame(coloc_result$summary)
    if (!nrow(summary) || !all(c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf") %in% names(summary))) {
      row_index <- row_index + 1
      coloc_rows[[row_index]] <- add_coloc_row(
        "COLOC_SUSIE", p12, if (abs(p12 - 1e-5) < .Machine$double.eps^0.5) "PRIMARY" else "SENSITIVITY",
        "NA", "NA", "NA", "NA", nrow(summary1), NA, NA, NA, NA, NA, "NA", NA, "COLOC_SUSIE_NO_ESTIMABLE_SIGNAL_PAIR"
      )
      next
    }
    result_table <- as.data.frame(coloc_result$results)
    for (index in seq_len(nrow(summary))) {
      pp_column <- if (nrow(summary) == 1) "SNP.PP.H4.abf" else paste0("SNP.PP.H4.row", index)
      top_index <- if (pp_column %in% names(result_table)) which.max(result_table[[pp_column]]) else integer(0)
      top_snp <- if (length(top_index)) result_table$snp[[top_index]] else "NA"
      top_pp <- if (length(top_index)) result_table[[pp_column]][[top_index]] else NA_real_
      row_index <- row_index + 1
      coloc_rows[[row_index]] <- add_coloc_row(
        "COLOC_SUSIE", p12, if (abs(p12 - 1e-5) < .Machine$double.eps^0.5) "PRIMARY" else "SENSITIVITY",
        paste0("L", summary$idx1[[index]]), paste0("L", summary$idx2[[index]]),
        summary$hit1[[index]], summary$hit2[[index]], summary$nsnps[[index]],
        summary$PP.H0.abf[[index]], summary$PP.H1.abf[[index]], summary$PP.H2.abf[[index]],
        summary$PP.H3.abf[[index]], summary$PP.H4.abf[[index]], top_snp, top_pp,
        if (fine_mapping_qc == "PASS") "COLOC_SUSIE_COMPLETE" else "COLOC_SUSIE_COMPLETE_QC_FLAGGED"
      )
      if (length(top_index)) {
        shared_index <- shared_index + 1
        shared_rows[[shared_index]] <- data.frame(
          comparison_id = comparison_id, pair_id = row$pair_id, locus_id = row$locus_id,
          comparison_type = row$comparison_type, coloc_method = "COLOC_SUSIE", p12 = p12,
          signal1 = paste0("L", summary$idx1[[index]]), signal2 = paste0("L", summary$idx2[[index]]),
          SNP = result_table$snp, SNP_PP_H4 = result_table[[pp_column]], stringsAsFactors = FALSE
        )
      }
    }
  }
} else if (row$single_signal_fallback_justification != "NOT_JUSTIFIED") {
  for (p12 in p12_grid) {
    coloc_result <- coloc::coloc.abf(dataset1, dataset2, p1 = p1, p2 = p2, p12 = p12)
    summary <- coloc_result$summary
    result_table <- as.data.frame(coloc_result$results)
    top_index <- which.max(result_table$SNP.PP.H4)
    row_index <- row_index + 1
    coloc_rows[[row_index]] <- add_coloc_row(
      "COLOC_ABF_FALLBACK", p12, if (abs(p12 - 1e-5) < .Machine$double.eps^0.5) "PRIMARY" else "SENSITIVITY",
      "SINGLE_SIGNAL_ASSUMPTION", "SINGLE_SIGNAL_ASSUMPTION", "NA", "NA", summary[["nsnps"]],
      summary[["PP.H0.abf"]], summary[["PP.H1.abf"]], summary[["PP.H2.abf"]],
      summary[["PP.H3.abf"]], summary[["PP.H4.abf"]], result_table$snp[[top_index]],
      result_table$SNP.PP.H4[[top_index]], if (fine_mapping_qc == "PASS") "COLOC_ABF_FALLBACK_COMPLETE" else "COLOC_ABF_FALLBACK_COMPLETE_QC_FLAGGED"
    )
    shared_index <- shared_index + 1
    shared_rows[[shared_index]] <- data.frame(
      comparison_id = comparison_id, pair_id = row$pair_id, locus_id = row$locus_id,
      comparison_type = row$comparison_type, coloc_method = "COLOC_ABF_FALLBACK", p12 = p12,
      signal1 = "SINGLE_SIGNAL_ASSUMPTION", signal2 = "SINGLE_SIGNAL_ASSUMPTION",
      SNP = result_table$snp, SNP_PP_H4 = result_table$SNP.PP.H4, stringsAsFactors = FALSE
    )
  }
} else {
  for (p12 in p12_grid) {
    row_index <- row_index + 1
    coloc_rows[[row_index]] <- add_coloc_row(
      "NONE", p12, if (abs(p12 - 1e-5) < .Machine$double.eps^0.5) "PRIMARY" else "SENSITIVITY",
      "NA", "NA", "NA", "NA", nrow(summary1), NA, NA, NA, NA, NA, "NA", NA,
      "NO_CREDIBLE_SET_AND_SINGLE_SIGNAL_FALLBACK_NOT_JUSTIFIED"
    )
  }
}

coloc_output <- if (length(coloc_rows)) do.call(rbind, coloc_rows) else empty_frame(coloc_fields)
shared_output <- if (length(shared_rows)) do.call(rbind, shared_rows) else empty_frame(shared_fields)
write_tsv(coloc_output[, coloc_fields], coloc_out)
write_tsv(shared_output[, shared_fields], shared_variant_out)
cat(sprintf(
  "SUSIE_COLOC_COMPARISON_OK comparison=%s variants=%d cs1=%d cs2=%d coloc_rows=%d qc=%s\n",
  comparison_id, nrow(summary1), if (is.null(fit1$sets$cs)) 0 else length(fit1$sets$cs),
  if (is.null(fit2$sets$cs)) 0 else length(fit2$sets$cs), nrow(coloc_output), fine_mapping_qc
))
