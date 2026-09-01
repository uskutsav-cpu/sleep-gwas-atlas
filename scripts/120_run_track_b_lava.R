#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(LAVA))

if (as.character(packageVersion("LAVA")) != "0.1.5") stop("LAVA 0.1.5 is required")

values <- list(
  input_info = "results/track_b/lava_input_info.tsv",
  sample_overlap = "results/track_b/lava_sample_overlap.txt",
  pair_manifest = "results/track_b/lava_pair_manifest.tsv",
  conditional_manifest = "results/track_b/local_conditional_manifest.tsv",
  runtime_policy = "results/track_b/lava_runtime_policy.tsv",
  checkpoint_dir = "results/track_b/checkpoints/lava",
  out_dir = "results/track_b/local", locus_start = 1L, locus_end = NA_integer_,
  finalize_only = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
i <- 1L
mapping <- c(
  "--input-info" = "input_info", "--sample-overlap" = "sample_overlap",
  "--pair-manifest" = "pair_manifest", "--conditional-manifest" = "conditional_manifest",
  "--runtime-policy" = "runtime_policy", "--checkpoint-dir" = "checkpoint_dir",
  "--out-dir" = "out_dir", "--locus-start" = "locus_start", "--locus-end" = "locus_end"
)
while (i <= length(args)) {
  if (args[[i]] == "--finalize-only") {
    values$finalize_only <- TRUE; i <- i + 1L; next
  }
  if (!(args[[i]] %in% names(mapping)) || i == length(args)) stop(paste("Unknown/incomplete argument:", args[[i]]))
  values[[unname(mapping[[args[[i]]]])]] <- args[[i + 1L]]; i <- i + 2L
}
values$locus_start <- as.integer(values$locus_start)
if (!is.na(values$locus_end)) values$locus_end <- as.integer(values$locus_end)

required <- c(values$input_info, values$sample_overlap, values$pair_manifest, values$conditional_manifest, values$runtime_policy)
if (!all(file.exists(required))) stop(paste("Missing Track B LAVA input:", paste(required[!file.exists(required)], collapse = ", ")))

policy_table <- read.delim(values$runtime_policy, stringsAsFactors = FALSE, check.names = FALSE)
policy <- setNames(as.list(policy_table$value), policy_table$key)
num <- function(key) as.numeric(policy[[key]])
integer <- function(key) as.integer(policy[[key]])
logical_value <- function(key) identical(tolower(policy[[key]]), "true")
if (policy$analysis_id != "track-b-v1.0-local" || policy$lava_version != "0.1.5") stop("Track B runtime policy drifted")

input_info <- read.delim(values$input_info, stringsAsFactors = FALSE, check.names = FALSE)
pairs <- read.delim(values$pair_manifest, stringsAsFactors = FALSE, check.names = FALSE)
conditioners <- read.delim(values$conditional_manifest, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(input_info) != integer("expected_traits") || length(unique(input_info$phenotype)) != integer("expected_traits")) stop("Input is not the exact eight-trait family")
if (nrow(pairs) != integer("expected_pairs") || !identical(pairs$pair_id, c("A", "B", "CONTROL"))) stop("Pair manifest is not A/B/CONTROL")
if (any(conditioners$selection_timing != "BEFORE_LOCAL_RESULT_ACCESS")) stop("Conditional covariates were not frozen before results")

contract <- suppressWarnings(system2("python3", c("scripts/119_track_b_lava_contract.py", "--run-fingerprint"), stdout = TRUE, stderr = TRUE))
if (!is.null(attr(contract, "status")) && attr(contract, "status") != 0L) stop(paste("Track B LAVA contract failed:", paste(contract, collapse = " | ")))
fingerprints <- contract[grepl("^[0-9a-f]{64}$", contract)]
if (length(fingerprints) != 1L) stop("Track B LAVA contract did not return one fingerprint")
fingerprint <- fingerprints[[1L]]

reference_prefix <- policy$reference_prefix
reference_files <- unlist(lapply(1:22, function(chr) paste0(reference_prefix, "_chr", chr, c(".info", ".bcor"))))
if (!all(file.exists(reference_files))) stop(paste("Complete LAVA reference missing", sum(!file.exists(reference_files)), "chromosome files"))
locus_file <- "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
loci <- read.loci(locus_file)
if (nrow(loci) != integer("expected_loci") || anyDuplicated(loci$LOC)) stop("Locus file is not the exact 2,495 family")
if (is.na(values$locus_end)) values$locus_end <- nrow(loci)
if (values$locus_start < 1L || values$locus_end > nrow(loci) || values$locus_start > values$locus_end) stop("Invalid locus range")
dir.create(values$checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(values$out_dir, recursive = TRUE, showWarnings = FALSE)

empty_univ <- function(locus_row, status, message = NA_character_, n_snps = NA_integer_, K = NA_integer_) {
  data.frame(LOC = locus_row$LOC, CHR = locus_row$CHR, START = locus_row$START, STOP = locus_row$STOP,
    phen = input_info$phenotype, h2.obs = NA_real_, h2.latent = NA_real_, ascertained = NA,
    p = NA_real_, analysis_status = status, error = message, n_snps = n_snps, K = K,
    stringsAsFactors = FALSE)
}
empty_bivar <- function() data.frame(
  LOC = numeric(0), CHR = numeric(0), START = numeric(0), STOP = numeric(0),
  pair_id = character(0), trait1 = character(0), trait2 = character(0),
  discovery_rg = numeric(0), discovery_SE = numeric(0), discovery_P = numeric(0), discovery_FDR = numeric(0),
  rho = numeric(0), rho.lower = numeric(0), rho.upper = numeric(0), r2 = numeric(0),
  r2.lower = numeric(0), r2.upper = numeric(0), p = numeric(0),
  analysis_status = character(0), error = character(0), stringsAsFactors = FALSE)

process_one <- function(index, input) {
  locus_row <- loci[index, , drop = FALSE]; started <- Sys.time()
  locus <- tryCatch(process.locus(locus_row, input, min.K = integer("min_K"),
    prune.thresh = num("prune_threshold"), max.prop.K = num("max_proportion_K"),
    drop.failed = TRUE, max.block.size = integer("max_block_size"), cap.estimates = logical_value("cap_estimates")),
    error = function(error) structure(list(message = conditionMessage(error)), class = "process_error"))
  if (inherits(locus, "process_error") || is.null(locus)) {
    message <- if (inherits(locus, "process_error")) locus$message else "process.locus returned NULL"
    univ <- empty_univ(locus_row, "LOCUS_PROCESS_FAILED", message); bivar <- empty_bivar()
    status <- "PROCESS_FAILED"; n_snps <- K <- NA_integer_
  } else {
    n_snps <- as.integer(locus$n.snps); K <- as.integer(locus$K)
    result <- tryCatch(run.univ(locus), error = function(error) structure(list(message = conditionMessage(error)), class = "univ_error"))
    if (inherits(result, "univ_error")) {
      univ <- empty_univ(locus_row, "UNIVARIATE_FAILED", result$message, n_snps, K)
      bivar <- empty_bivar(); status <- "UNIVARIATE_FAILED"
    } else {
      univ <- empty_univ(locus_row, "PHENOTYPE_DROPPED", "phenotype removed during locus processing", n_snps, K)
      matched <- match(result$phen, univ$phen)
      univ$h2.obs[matched] <- result$h2.obs; univ$h2.latent[matched] <- result$h2.latent
      univ$ascertained[matched] <- result$ascertained; univ$p[matched] <- result$p
      univ$analysis_status[matched] <- "TESTED"; univ$error[matched] <- NA_character_
      eligible <- univ$phen[univ$analysis_status == "TESTED" & is.finite(univ$p) & univ$p <= num("univariate_p_threshold")]
      eligible_pairs <- pairs[pairs$trait1 %in% eligible & pairs$trait2 %in% eligible, , drop = FALSE]
      rows <- vector("list", nrow(eligible_pairs))
      if (nrow(eligible_pairs)) for (j in seq_len(nrow(eligible_pairs))) {
        pair <- eligible_pairs[j, , drop = FALSE]; set.seed(integer("random_seed") + index * 1000L + j)
        local <- tryCatch(run.bivar(locus, phenos = c(pair$trait1, pair$trait2)),
          error = function(error) structure(list(message = conditionMessage(error)), class = "bivar_error"))
        base <- data.frame(LOC = locus_row$LOC, CHR = locus_row$CHR, START = locus_row$START, STOP = locus_row$STOP,
          pair_id = pair$pair_id, trait1 = pair$trait1, trait2 = pair$trait2,
          discovery_rg = as.numeric(pair$discovery_rg), discovery_SE = as.numeric(pair$discovery_SE),
          discovery_P = as.numeric(pair$discovery_P), discovery_FDR = as.numeric(pair$discovery_FDR), stringsAsFactors = FALSE)
        if (inherits(local, "bivar_error") || is.null(local) || nrow(local) != 1L) {
          message <- if (inherits(local, "bivar_error")) local$message else "run.bivar returned no unique result"
          rows[[j]] <- cbind(base, rho = NA_real_, rho.lower = NA_real_, rho.upper = NA_real_, r2 = NA_real_,
            r2.lower = NA_real_, r2.upper = NA_real_, p = NA_real_, analysis_status = "BIVARIATE_FAILED", error = message)
        } else rows[[j]] <- cbind(base, rho = local$rho, rho.lower = local$rho.lower, rho.upper = local$rho.upper,
          r2 = local$r2, r2.lower = local$r2.lower, r2.upper = local$r2.upper, p = local$p,
          analysis_status = "TESTED", error = NA_character_)
      }
      bivar <- if (length(rows)) do.call(rbind, rows) else empty_bivar(); status <- "PROCESSED"
    }
  }
  status_row <- data.frame(LOC = locus_row$LOC, CHR = locus_row$CHR, START = locus_row$START, STOP = locus_row$STOP,
    status = status, n_snps = n_snps, K = K, univariate_tested = sum(univ$analysis_status == "TESTED"),
    eligible_bivariate_pairs = nrow(bivar), elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    stringsAsFactors = FALSE)
  list(fingerprint = fingerprint, status = status_row, univ = univ, bivar = bivar)
}

input <- NULL
if (!values$finalize_only) {
  input <- process.input(input.info.file = values$input_info, sample.overlap.file = values$sample_overlap,
    ref.prefix = reference_prefix, phenos = input_info$phenotype)
  for (index in seq.int(values$locus_start, values$locus_end)) {
    checkpoint <- file.path(values$checkpoint_dir, sprintf("locus_%04d.rds", index))
    if (file.exists(checkpoint)) {
      existing <- readRDS(checkpoint)
      if (!identical(existing$fingerprint, fingerprint)) stop(paste("Stale checkpoint:", checkpoint))
      next
    }
    result <- process_one(index, input); temporary <- paste0(checkpoint, ".tmp")
    saveRDS(result, temporary, version = 3)
    if (!file.rename(temporary, checkpoint)) stop(paste("Could not publish checkpoint:", checkpoint))
  }
}

checkpoint_files <- file.path(values$checkpoint_dir, sprintf("locus_%04d.rds", seq_len(nrow(loci))))
if (!all(file.exists(checkpoint_files))) {
  message(sprintf("Track B LAVA range complete; %d/%d checkpoints exist; no canonical outputs published", sum(file.exists(checkpoint_files)), nrow(loci)))
  quit(save = "no", status = 3L)
}
results <- lapply(checkpoint_files, readRDS)
if (any(!vapply(results, function(x) identical(x$fingerprint, fingerprint), logical(1)))) stop("Stale checkpoint fingerprint")
status <- do.call(rbind, lapply(results, `[[`, "status")); univ <- do.call(rbind, lapply(results, `[[`, "univ"))
bivar <- do.call(rbind, lapply(results, `[[`, "bivar"))
if (nrow(status) != integer("expected_loci") || nrow(univ) != integer("planned_univariate_tests")) stop("Incomplete checkpoint family")
univ$univariate_test_family_n <- integer("planned_univariate_tests"); univ$univariate_p_threshold <- num("univariate_p_threshold")
univ$p_bonferroni <- pmin(1, univ$p * integer("planned_univariate_tests")); univ$p_fdr <- NA_real_
tested_univ <- which(univ$analysis_status == "TESTED" & is.finite(univ$p)); univ$p_fdr[tested_univ] <- p.adjust(univ$p[tested_univ], method = "BH")
bivar$bivariate_test_family_n <- rep.int(nrow(bivar), nrow(bivar)); bivar$p_fdr <- rep(NA_real_, nrow(bivar))
tested_bivar <- which(bivar$analysis_status == "TESTED" & is.finite(bivar$p)); bivar$p_fdr[tested_bivar] <- p.adjust(bivar$p[tested_bivar], method = "BH")
bivar$fdr_significant <- !is.na(bivar$p_fdr) & bivar$p_fdr <= num("bivariate_fdr_alpha")

empty_conditional <- function() data.frame(
  LOC = numeric(0), CHR = numeric(0), START = numeric(0), STOP = numeric(0), pair_id = character(0),
  trait1 = character(0), trait2 = character(0), covariates = character(0), pcor = numeric(0),
  ci.lower = numeric(0), ci.upper = numeric(0), p = numeric(0), r2.trait1_z = numeric(0),
  r2.trait2_z = numeric(0), analysis_status = character(0), error = character(0), stringsAsFactors = FALSE)
conditional_rows <- list(); targets <- bivar[bivar$fdr_significant & bivar$pair_id %in% c("A", "B"), , drop = FALSE]
if (nrow(targets) && is.null(input)) input <- process.input(input.info.file = values$input_info,
  sample.overlap.file = values$sample_overlap, ref.prefix = reference_prefix, phenos = input_info$phenotype)
if (nrow(targets)) for (j in seq_len(nrow(targets))) {
  target <- targets[j, , drop = FALSE]; covars <- conditioners$covariate[conditioners$pair_id == target$pair_id & conditioners$covariate != "NONE"]
  base <- data.frame(LOC = target$LOC, CHR = target$CHR, START = target$START, STOP = target$STOP,
    pair_id = target$pair_id, trait1 = target$trait1, trait2 = target$trait2, covariates = paste(covars, collapse = ";"), stringsAsFactors = FALSE)
  local_univ <- univ[univ$LOC == target$LOC & univ$phen %in% covars, , drop = FALSE]
  eligible <- length(covars) > 0L && nrow(local_univ) == length(covars) && all(local_univ$analysis_status == "TESTED" & is.finite(local_univ$p) & local_univ$p <= num("univariate_p_threshold"))
  if (!eligible) {
    conditional_rows[[j]] <- cbind(base, pcor = NA_real_, ci.lower = NA_real_, ci.upper = NA_real_, p = NA_real_,
      r2.trait1_z = NA_real_, r2.trait2_z = NA_real_, analysis_status = "CONDITIONER_LOCAL_H2_INELIGIBLE", error = "At least one predeclared conditioner did not pass the same-locus univariate gate")
    next
  }
  index <- match(target$LOC, loci$LOC); locus <- tryCatch(process.locus(loci[index, , drop = FALSE], input,
    min.K = integer("min_K"), prune.thresh = num("prune_threshold"), max.prop.K = num("max_proportion_K"),
    drop.failed = TRUE, max.block.size = integer("max_block_size"), cap.estimates = logical_value("cap_estimates")),
    error = function(error) structure(list(message = conditionMessage(error)), class = "conditional_error"))
  set.seed(integer("random_seed") + index * 1000L + 500L + j)
  partial <- if (inherits(locus, "conditional_error") || is.null(locus)) locus else tryCatch(
    run.pcor(locus, target = c(target$trait1, target$trait2), phenos = c(target$trait1, target$trait2, covars), max.r2 = num("conditional_max_r2")),
    error = function(error) structure(list(message = conditionMessage(error)), class = "conditional_error"))
  if (inherits(partial, "conditional_error") || is.null(partial) || nrow(partial) != 1L) {
    message <- if (inherits(partial, "conditional_error")) partial$message else "run.pcor returned no unique result"
    conditional_rows[[j]] <- cbind(base, pcor = NA_real_, ci.lower = NA_real_, ci.upper = NA_real_, p = NA_real_,
      r2.trait1_z = NA_real_, r2.trait2_z = NA_real_, analysis_status = "CONDITIONAL_FAILED", error = message)
  } else conditional_rows[[j]] <- cbind(base, pcor = partial$pcor, ci.lower = partial$ci.lower, ci.upper = partial$ci.upper,
    p = partial$p, r2.trait1_z = partial[["r2.phen1_z"]], r2.trait2_z = partial[["r2.phen2_z"]],
    analysis_status = "TESTED", error = NA_character_)
}
conditional <- if (length(conditional_rows)) do.call(rbind, conditional_rows) else empty_conditional()
conditional$conditional_test_family_n <- rep.int(sum(conditional$analysis_status == "TESTED"), nrow(conditional))
conditional$p_fdr <- rep(NA_real_, nrow(conditional)); tested_cond <- which(conditional$analysis_status == "TESTED" & is.finite(conditional$p))
conditional$p_fdr[tested_cond] <- p.adjust(conditional$p[tested_cond], method = "BH")
conditional$fdr_significant <- !is.na(conditional$p_fdr) & conditional$p_fdr <= num("conditional_fdr_alpha")
status$analysis_fingerprint <- fingerprint; status$lava_version <- as.character(packageVersion("LAVA")); status$reference_prefix <- reference_prefix

write_atomic <- function(table, path) {
  temporary <- paste0(path, ".tmp"); write.table(table, temporary, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  if (!file.rename(temporary, path)) stop(paste("Could not publish:", path))
}
write_atomic(status, file.path(values$out_dir, "lava_locus_status.tsv")); write_atomic(univ, file.path(values$out_dir, "lava_univariate.tsv"))
write_atomic(bivar, file.path(values$out_dir, "lava_bivariate.tsv")); write_atomic(conditional, file.path(values$out_dir, "lava_conditional.tsv"))
validation <- system2("python3", c("scripts/121_validate_track_b_lava.py", "--seal-results"), stdout = TRUE, stderr = TRUE)
if (!is.null(attr(validation, "status")) && attr(validation, "status") != 0L) stop(paste("Track B LAVA validation failed:", paste(validation, collapse = " | ")))
message(sprintf("Published Track B LAVA: %d loci, %d univariate, %d bivariate, %d conditional rows", nrow(status), nrow(univ), nrow(bivar), nrow(conditional)))
