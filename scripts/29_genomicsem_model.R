#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
fail <- function(...) stop(paste0(...), call. = FALSE)

args <- commandArgs(trailingOnly = TRUE)
defaults <- list(
  discovery = "results/tables/genomicsem_discovery_odd.rds",
  validation = "results/tables/genomicsem_validation_even.rds",
  inclusion = "results/tables/genomicsem_trait_inclusion.tsv",
  out_dir = "results/tables"
)
if (length(args) %% 2L) fail("Arguments must be --name value pairs")
if (length(args)) for (i in seq(1L, length(args), by = 2L)) {
  key <- gsub("-", "_", sub("^--", "", args[[i]]))
  if (!key %in% names(defaults)) fail("Unknown argument: ", args[[i]])
  defaults[[key]] <- args[[i + 1L]]
}

# lavaan 0.6-19 calls parallel::detectCores() before applying ncpus. Some
# restricted runtimes return NA. Pinning that discovery to one core changes
# scheduling only and makes model estimation deterministic.
if (is.na(parallel::detectCores())) {
  parallel_namespace <- asNamespace("parallel")
  unlockBinding("detectCores", parallel_namespace)
  assign("detectCores", function(all.tests = FALSE, logical = TRUE) 1L,
         envir = parallel_namespace)
  lockBinding("detectCores", parallel_namespace)
}

suppressPackageStartupMessages(library(GenomicSEM))
suppressPackageStartupMessages(library(psych))
if (getRversion() != "4.3.3" || as.character(packageVersion("GenomicSEM")) != "0.0.5" ||
    as.character(packageVersion("psych")) != "2.5.6") {
  fail("Pinned R 4.3.3, GenomicSEM 0.0.5, and psych 2.5.6 are required")
}

odd <- readRDS(defaults$discovery)
even <- readRDS(defaults$validation)
inclusion <- read.delim(defaults$inclusion, check.names = FALSE)
traits <- inclusion$trait_id[inclusion$include_genomic_sem == TRUE]
if (length(traits) != 42L || !identical(colnames(odd$S), traits) ||
    !identical(colnames(even$S), traits)) fail("Split covariance traits differ from the 42-trait QC ledger")
for (object in list(odd, even)) {
  if (!identical(dim(object$S), c(42L, 42L)) || !identical(dim(object$V), c(903L, 903L)) ||
      any(!is.finite(object$S)) || any(!is.finite(object$V))) fail("Invalid split covariance object")
}

odd_rg <- odd$S_Stand
smoothed_rg <- psych::cor.smooth(odd_rg)
smooth_frobenius <- sqrt(sum((smoothed_rg - odd_rg)^2))
eigenvalues <- eigen(smoothed_rg, symmetric = TRUE, only.values = TRUE)$values

efa_rows <- list()
efa_fits <- list()
for (method in c("ml", "minres")) {
  for (factors in 1:10) {
    rotation <- if (factors == 1L) "none" else "oblimin"
    fit <- suppressWarnings(psych::fa(
      smoothed_rg, nfactors = factors, fm = method, rotate = rotation,
      SMC = FALSE, warnings = FALSE, max.iter = 1000
    ))
    loadings <- unclass(fit$loadings)
    primary <- max.col(abs(loadings))
    maximum <- apply(abs(loadings), 1L, max)
    efa_rows[[length(efa_rows) + 1L]] <- data.frame(
      method = method,
      factors = factors,
      rms = fit$rms,
      crms = fit$crms,
      off_diagonal_fit = fit$fit.off,
      heywood_communality_gt_1 = sum(fit$communality > 1),
      minimum_uniqueness = min(fit$uniquenesses),
      indicators_below_abs_0_3 = sum(maximum < 0.3),
      smallest_primary_factor_size = min(tabulate(primary, nbins = factors)),
      smoothing_frobenius = smooth_frobenius
    )
    efa_fits[[paste(method, factors, sep = "_")]] <- fit
  }
}
efa_table <- do.call(rbind, efa_rows)

candidate_models <- list()
candidate_meta <- list()
for (factors in 3:7) {
  fit <- efa_fits[[paste("ml", factors, sep = "_")]]
  loadings <- unclass(fit$loadings)
  primary <- max.col(abs(loadings))
  primary_sizes <- tabulate(primary, nbins = factors)
  if (all(primary_sizes >= 3L)) {
    lines <- vapply(seq_len(factors), function(index) {
      paste0("F", index, " =~ ", paste(rownames(loadings)[primary == index], collapse = " + "))
    }, character(1))
    id <- paste0("ml_k", factors, "_primary")
    candidate_models[[id]] <- paste(lines, collapse = "\n")
    candidate_meta[[id]] <- list(factors = factors, indicators = nrow(loadings),
                                 cross_loaders = 0L, boundary = "none")
  }

  salient <- abs(loadings) >= 0.30
  salient_sizes <- colSums(salient)
  if (all(salient_sizes >= 3L)) {
    lines <- vapply(seq_len(factors), function(index) {
      paste0("F", index, " =~ ", paste(rownames(loadings)[salient[, index]], collapse = " + "))
    }, character(1))
    id <- paste0("ml_k", factors, "_salient_crossload")
    candidate_models[[id]] <- paste(lines, collapse = "\n")
    candidate_meta[[id]] <- list(
      factors = factors, indicators = sum(rowSums(salient) > 0L),
      cross_loaders = sum(rowSums(salient) > 1L), boundary = "none"
    )
  }
}

# This candidate is determined only from odd-chromosome discovery results.
# Boundary residuals are fixed for discovery Heywood traits before the even
# chromosome object is inspected.
base_id <- "ml_k6_salient_crossload"
if (base_id %in% names(candidate_models)) {
  discovery_fit <- GenomicSEM::usermodel(
    covstruc = odd, estimation = "DWLS", model = candidate_models[[base_id]],
    std.lv = TRUE, fix_resid = FALSE
  )
  discovery_results <- discovery_fit$results
  discovery_heywood <- discovery_results$lhs[
    discovery_results$op == "~~" &
      discovery_results$lhs == discovery_results$rhs &
      discovery_results$lhs %in% traits &
      discovery_results$Unstand_Est < 0
  ]
  if (length(discovery_heywood)) {
    id <- paste0(base_id, "_discovery_boundary")
    constraints <- paste0(discovery_heywood, " ~~ 0.0001*", discovery_heywood)
    candidate_models[[id]] <- paste(c(candidate_models[[base_id]], constraints), collapse = "\n")
    candidate_meta[[id]] <- candidate_meta[[base_id]]
    candidate_meta[[id]]$boundary <- paste(discovery_heywood, collapse = ",")
  }
}

fit_rows <- list()
loading_rows <- list()
syntax_rows <- list()
for (id in names(candidate_models)) {
  model <- candidate_models[[id]]
  captured_warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      GenomicSEM::usermodel(
        covstruc = even, estimation = "DWLS", model = model,
        std.lv = TRUE, fix_resid = FALSE
      ),
      warning = function(warning) {
        captured_warnings <<- c(captured_warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) error
  )
  meta <- candidate_meta[[id]]
  if (inherits(fit, "error")) {
    fit_rows[[length(fit_rows) + 1L]] <- data.frame(
      model_id = id, factors = meta$factors, indicators = meta$indicators,
      discovery_cross_loaders = meta$cross_loaders,
      discovery_boundary_residuals = meta$boundary,
      chisq = NA, df = NA, p_chisq = NA, AIC = NA, BIC = NA, CFI = NA, SRMR = NA,
      negative_observed_residuals = NA, warning_count = length(captured_warnings),
      validation_status = "FIT_ERROR", validation_reason = conditionMessage(fit)
    )
    next
  }
  model_fit <- fit$modelfit[1, ]
  result <- fit$results
  negative <- result$lhs[
    result$op == "~~" & result$lhs == result$rhs & result$lhs %in% traits &
      result$Unstand_Est < -1e-10
  ]
  cfi <- as.numeric(model_fit$CFI)
  srmr <- as.numeric(model_fit$SRMR)
  passes <- is.finite(cfi) && is.finite(srmr) && cfi >= 0.90 && srmr <= 0.10 && !length(negative)
  reasons <- c()
  if (!is.finite(cfi) || cfi < 0.90) reasons <- c(reasons, "CFI_below_0.90")
  if (!is.finite(srmr) || srmr > 0.10) reasons <- c(reasons, "SRMR_above_0.10")
  if (length(negative)) reasons <- c(reasons, paste0("negative_residuals:", paste(negative, collapse = ",")))
  if (any(grepl("difference greater than .025", captured_warnings, fixed = TRUE))) {
    reasons <- c(reasons, "smoothing_Z_change_above_0.025_warning")
  }
  fit_rows[[length(fit_rows) + 1L]] <- data.frame(
    model_id = id, factors = meta$factors, indicators = meta$indicators,
    discovery_cross_loaders = meta$cross_loaders,
    discovery_boundary_residuals = meta$boundary,
    chisq = as.numeric(model_fit$chisq), df = as.numeric(model_fit$df),
    p_chisq = as.numeric(model_fit$p_chisq), AIC = as.numeric(model_fit$AIC), BIC = NA,
    CFI = cfi, SRMR = srmr, negative_observed_residuals = length(negative),
    warning_count = length(captured_warnings),
    validation_status = if (passes) "VALIDATED" else "NOT_VALIDATED",
    validation_reason = if (passes) "pass" else paste(reasons, collapse = ";")
  )
  parameter_rows <- result[result$op == "=~", , drop = FALSE]
  if (nrow(parameter_rows)) {
    loading_rows[[length(loading_rows) + 1L]] <- data.frame(
      model_id = id, validation_status = if (passes) "VALIDATED" else "NOT_VALIDATED",
      factor = parameter_rows$lhs, trait_id = parameter_rows$rhs,
      estimate = parameter_rows$Unstand_Est, se = parameter_rows$Unstand_SE,
      standardized_loading = parameter_rows$STD_All,
      p = parameter_rows$p_value
    )
  }
  syntax_rows[[length(syntax_rows) + 1L]] <- data.frame(
    model_id = id, model = gsub("[\r\n\t]+", "; ", model),
    warnings = gsub("[\r\n\t]+", " ", paste(unique(captured_warnings), collapse = " | "))
  )
}

fit_table <- do.call(rbind, fit_rows)
loading_table <- if (length(loading_rows)) do.call(rbind, loading_rows) else data.frame()
syntax_table <- if (length(syntax_rows)) do.call(rbind, syntax_rows) else data.frame()
dir.create(defaults$out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(efa_table, file.path(defaults$out_dir, "genomic_sem_efa_models.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write.table(fit_table, file.path(defaults$out_dir, "genomic_sem_model_fit.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write.table(loading_table, file.path(defaults$out_dir, "genomic_sem_factor_loadings.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write.table(syntax_table, file.path(defaults$out_dir, "genomic_sem_model_syntax.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write.table(data.frame(
  metric = c("odd_even_off_diagonal_pearson", "odd_even_off_diagonal_spearman",
             "odd_even_median_absolute_rg_difference", "odd_even_sign_agreement",
             "odd_smoothing_frobenius", "odd_smoothed_eigenvalues"),
  value = c(
    cor(odd$S_Stand[lower.tri(odd$S_Stand)], even$S_Stand[lower.tri(even$S_Stand)]),
    cor(odd$S_Stand[lower.tri(odd$S_Stand)], even$S_Stand[lower.tri(even$S_Stand)], method = "spearman"),
    median(abs(odd$S_Stand[lower.tri(odd$S_Stand)] - even$S_Stand[lower.tri(even$S_Stand)])),
    mean(sign(odd$S_Stand[lower.tri(odd$S_Stand)]) == sign(even$S_Stand[lower.tri(even$S_Stand)])),
    smooth_frobenius, paste(eigenvalues, collapse = ",")
  )
), file.path(defaults$out_dir, "genomic_sem_split_diagnostics.tsv"),
sep = "\t", row.names = FALSE, quote = FALSE)

validated <- fit_table$model_id[fit_table$validation_status == "VALIDATED"]
if (length(validated)) {
  cat("Validated held-out GenomicSEM model(s): ", paste(validated, collapse = ", "), "\n", sep = "")
} else {
  cat("No candidate GenomicSEM model passed held-out CFI/SRMR and residual admissibility criteria\n")
}
