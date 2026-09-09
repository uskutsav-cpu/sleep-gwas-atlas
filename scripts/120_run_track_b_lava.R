#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(LAVA))

if (as.character(packageVersion("LAVA")) != "0.1.5") stop("LAVA 0.1.5 is required")
runtime_check <- suppressWarnings(system2(".r-env/bin/Rscript", "scripts/133_validate_track_b_lava_runtime.R", stdout = TRUE, stderr = TRUE))
if (!is.null(attr(runtime_check, "status")) && attr(runtime_check, "status") != 0L) stop(paste("Pinned LAVA implementation failed verification:", paste(runtime_check, collapse = " | ")))
if (sum(grepl("^TRACK_B_LAVA_RUNTIME\\tstatus=PASS\\t", runtime_check)) != 1L) stop("Pinned LAVA runtime did not emit one semantic PASS marker")

values <- list(
  input_info = "results/track_b/lava_input_info.tsv",
  sample_overlap = "results/track_b/lava_sample_overlap.txt",
  pair_manifest = "results/track_b/lava_pair_manifest.tsv",
  conditional_manifest = "results/track_b/local_conditional_manifest.tsv",
  runtime_policy = "results/track_b/lava_runtime_policy.tsv",
  checkpoint_dir = "results/track_b/checkpoints/lava",
  out_dir = NA_character_, phase = NA_character_, locus_index = NA_integer_,
  fingerprint = NA_character_, worker_output = NA_character_
)

args <- commandArgs(trailingOnly = TRUE)
i <- 1L
mapping <- c("--phase" = "phase", "--locus-index" = "locus_index", "--fingerprint" = "fingerprint",
  "--worker-output" = "worker_output")
while (i <= length(args)) {
  if (!(args[[i]] %in% names(mapping)) || i == length(args)) stop(paste("Unknown/incomplete argument:", args[[i]]))
  values[[unname(mapping[[args[[i]]]])]] <- args[[i + 1L]]
  i <- i + 2L
}
values$locus_index <- as.integer(values$locus_index)
allowed_phases <- c("discovery", "aggregate-discovery", "conditional", "finalize")
if (is.na(values$phase) || !(values$phase %in% allowed_phases)) stop("--phase must be discovery, aggregate-discovery, conditional, or finalize")

required <- c(values$input_info, values$sample_overlap, values$pair_manifest, values$conditional_manifest, values$runtime_policy)
if (!all(file.exists(required))) stop(paste("Missing Track B LAVA input:", paste(required[!file.exists(required)], collapse = ", ")))

policy_table <- read.delim(values$runtime_policy, stringsAsFactors = FALSE, check.names = FALSE)
policy <- setNames(as.list(policy_table$value), policy_table$key)
num <- function(key) as.numeric(policy[[key]])
integer <- function(key) as.integer(policy[[key]])
logical_value <- function(key) identical(tolower(policy[[key]]), "true")
if (policy$analysis_id != "track-b-v1.0-local" || policy$lava_version != "0.1.5") stop("Track B runtime policy drifted")
if (policy$ram_execution_unit != "WHOLE_PREDECLARED_LAVA_LOCUS" || integer("maximum_loci_per_process") != 1L) stop("Track B RAM-aware execution unit drifted")

input_info <- read.delim(values$input_info, stringsAsFactors = FALSE, check.names = FALSE)
pairs <- read.delim(values$pair_manifest, stringsAsFactors = FALSE, check.names = FALSE)
conditioners <- read.delim(values$conditional_manifest, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(input_info) != integer("expected_traits") || length(unique(input_info$phenotype)) != integer("expected_traits")) stop("Input is not the exact eight-trait family")
if (nrow(pairs) != integer("expected_pairs") || !identical(pairs$pair_id, c("A", "B", "CONTROL"))) stop("Pair manifest is not A/B/CONTROL")
if (any(conditioners$selection_timing != "BEFORE_LOCAL_RESULT_ACCESS")) stop("Conditional covariates were not frozen before results")
if (!identical(conditioners$conditional_model_id, c("A_BMI_ONLY", "A_SLEEP_APNEA_ONLY", "B_MDD_ONLY", "NONE"))) stop("Conditional model family drifted")

if (is.na(values$fingerprint)) {
  contract <- suppressWarnings(system2("python3", c("scripts/119_track_b_lava_contract.py", "--execution-fingerprint"), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(contract, "status")) && attr(contract, "status") != 0L) stop(paste("Track B LAVA contract failed:", paste(contract, collapse = " | ")))
  fingerprints <- contract[grepl("^[0-9a-f]{64}$", contract)]
  if (length(fingerprints) != 1L) stop("Track B LAVA contract did not return one fingerprint")
  values$fingerprint <- fingerprints[[1L]]
}
if (!grepl("^[0-9a-f]{64}$", values$fingerprint)) stop("Invalid Track B LAVA execution fingerprint")
fingerprint <- values$fingerprint
verify_current_fingerprint <- function() {
  contract <- suppressWarnings(system2("python3", c("scripts/119_track_b_lava_contract.py", "--quick-execution-fingerprint"), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(contract, "status")) && attr(contract, "status") != 0L) stop(paste("Track B LAVA quick contract failed:", paste(contract, collapse = " | ")))
  observed <- contract[grepl("^[0-9a-f]{64}$", contract)]
  if (length(observed) != 1L || !identical(observed[[1L]], fingerprint)) stop("Live Track B inputs/runtime differ from the supplied fingerprint")
  invisible(TRUE)
}
verify_current_fingerprint()
values$checkpoint_dir <- file.path(values$checkpoint_dir, fingerprint)
if (is.na(values$worker_output) || !nzchar(values$worker_output)) stop("--worker-output is required")
worker_parent <- normalizePath(dirname(values$worker_output), mustWork = TRUE)
checkpoint_root <- normalizePath(values$checkpoint_dir, mustWork = TRUE)
if (!startsWith(paste0(worker_parent, "/"), paste0(checkpoint_root, "/")) ||
    basename(dirname(worker_parent)) != ".attempts") {
  stop("Worker output must be inside the current fingerprint's supervised attempt directory")
}
values$worker_output <- file.path(worker_parent, basename(values$worker_output))

result_lock <- "results/track_b/local/lava_results.provenance.json"
if (file.exists(result_lock)) {
  verification <- system2("python3", c("scripts/119_track_b_lava_contract.py", "--verify-results"), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(verification, "status")) && attr(verification, "status") != 0L) stop(paste("Existing Track B LAVA result seal is invalid:", paste(verification, collapse = " | ")))
  if (values$phase != "finalize") stop("Track B LAVA results are already sealed; no non-final worker may run")
}

reference_prefix <- policy$reference_prefix
reference_files <- unlist(lapply(1:22, function(chr) paste0(reference_prefix, "_chr", chr, c(".info", ".bcor"))))
if (!all(file.exists(reference_files))) stop(paste("Complete LAVA reference missing", sum(!file.exists(reference_files)), "chromosome files"))
locus_file <- "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
loci <- read.loci(locus_file)
if (nrow(loci) != integer("expected_loci") || anyDuplicated(loci$LOC)) stop("Locus file is not the exact 2,495 family")
if (values$phase %in% c("discovery", "conditional") && (is.na(values$locus_index) || values$locus_index < 1L || values$locus_index > nrow(loci))) stop("A valid --locus-index is required for a locus worker")
if (!dir.exists(values$checkpoint_dir)) stop("Fingerprint checkpoint namespace was not created by the supervisor")

finite_or_na <- function(value) {
  if (length(value) == 1L && is.finite(value)) as.numeric(value) else NA_real_
}

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
  local_covariance = numeric(0), rho = numeric(0), rho.lower = numeric(0), rho.upper = numeric(0), r2 = numeric(0),
  r2.lower = numeric(0), r2.upper = numeric(0), p = numeric(0),
  analysis_status = character(0), error = character(0), stringsAsFactors = FALSE)

empty_conditional <- function() data.frame(
  LOC = numeric(0), CHR = numeric(0), START = numeric(0), STOP = numeric(0), pair_id = character(0),
  trait1 = character(0), trait2 = character(0), conditional_model_id = character(0),
  covariates = character(0), pcor = numeric(0), ci.lower = numeric(0), ci.upper = numeric(0), p = numeric(0),
  r2.trait1_z = numeric(0), r2.trait2_z = numeric(0), analysis_status = character(0), error = character(0),
  stringsAsFactors = FALSE)

checkpoint_path <- function(phase, index) {
  file.path(values$checkpoint_dir, phase, sprintf("locus_%04d", index), "result.rds")
}

validate_checkpoint <- function(value, phase, index) {
  if (!is.list(value) || !identical(value$fingerprint, fingerprint) || !identical(value$phase, phase) ||
      !identical(as.integer(value$locus_index), as.integer(index)) ||
      !identical(as.character(value$LOC), as.character(loci$LOC[index]))) {
    stop(paste("Stale or malformed", phase, "checkpoint for locus index", index))
  }
  value
}

publish_rds_no_replace <- function(value, path) {
  if (file.exists(path)) stop(paste("Checkpoint already exists; overwrite is forbidden:", path))
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(value, temporary, version = 3)
  if (!file.link(temporary, path)) stop(paste("Could not atomically publish checkpoint without replacement:", path))
  invisible(path)
}

write_table_no_replace <- function(table, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  write.table(table, temporary, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  if (file.exists(path)) {
    if (!identical(unname(tools::md5sum(temporary)), unname(tools::md5sum(path)))) stop(paste("Existing output differs; overwrite is forbidden:", path))
    return(invisible(path))
  }
  if (!file.link(temporary, path)) stop(paste("Could not atomically publish output without replacement:", path))
  invisible(path)
}

worker_marker <- function(phase, index, locus, chromosome, n_snps, pair, qc, construction_complete, output) {
  cat(paste("TRACK_B_LAVA_WORKER", paste0("phase=", phase), paste0("index=", index),
    paste0("locus=", locus), paste0("chromosome=", chromosome), paste0("n_snps=", n_snps),
    paste0("pair=", pair), paste0("qc=", qc), paste0("construction_complete=", construction_complete),
    paste0("output=", output), sep = "\t"), "\n", sep = "")
}

univariate_result_error <- function(result, locus) {
  required <- c("phen", "h2.obs", "p")
  if (!is.data.frame(result) || !all(required %in% names(result)) || !nrow(result)) return("run.univ returned no valid data-frame result")
  phen <- as.character(result$phen)
  if (anyNA(phen) || anyDuplicated(phen) || any(!(phen %in% locus$phenos))) return("run.univ returned missing, duplicate, or unknown phenotypes")
  if (any(!is.finite(result$h2.obs) | result$h2.obs < 0) || any(!is.finite(result$p) | result$p < 0 | result$p > 1)) return("run.univ returned a non-finite/out-of-range h2 or p-value")
  binary <- as.logical(locus$binary[phen])
  if (any(binary)) {
    if (!("h2.latent" %in% names(result)) || any(!is.finite(result$h2.latent[binary]) | result$h2.latent[binary] < 0)) return("run.univ returned an invalid binary-trait latent-scale h2")
  }
  if ("h2.latent" %in% names(result) && any(!binary) && any(!is.na(result$h2.latent[!binary]))) return("run.univ returned a latent-scale h2 for a continuous trait")
  NA_character_
}

scoped_input <- function(locus_row) {
  chromosome <- as.integer(locus_row$CHR)
  chromosome_prefix <- paste0(reference_prefix, "_chr", chromosome)
  chromosome_input_info <- file.path("results/track_b/lava_chromosome_inputs",
    sprintf("chr%02d", chromosome), "lava_input_info.tsv")
  if (!file.exists(chromosome_input_info)) stop(paste("Verified chromosome input-info is missing:", chromosome_input_info))
  process.input(input.info.file = chromosome_input_info, sample.overlap.file = values$sample_overlap,
    ref.prefix = chromosome_prefix, phenos = input_info$phenotype)
}

process_discovery_locus <- function(index, input) {
  locus_row <- loci[index, , drop = FALSE]
  started <- Sys.time()
  locus <- tryCatch(process.locus(locus_row, input, min.K = integer("min_K"),
    prune.thresh = num("prune_threshold"), max.prop.K = num("max_proportion_K"),
    drop.failed = TRUE, max.block.size = integer("max_block_size"), cap.estimates = logical_value("cap_estimates")),
    error = function(error) structure(list(message = conditionMessage(error)), class = "process_error"))
  if (inherits(locus, "process_error") || is.null(locus)) {
    message_text <- if (inherits(locus, "process_error")) locus$message else "process.locus returned NULL"
    univ <- empty_univ(locus_row, "LOCUS_PROCESS_FAILED", message_text)
    bivar <- empty_bivar()
    status_value <- "PROCESS_FAILED"
    n_snps <- K <- NA_integer_
  } else {
    n_snps <- as.integer(locus$n.snps)
    K <- as.integer(locus$K)
    result <- tryCatch(run.univ(locus), error = function(error) structure(list(message = conditionMessage(error)), class = "univ_error"))
    result_error <- if (inherits(result, "univ_error")) result$message else univariate_result_error(result, locus)
    if (!is.na(result_error)) {
      univ <- empty_univ(locus_row, "UNIVARIATE_FAILED", result_error, n_snps, K)
      bivar <- empty_bivar()
      status_value <- "UNIVARIATE_FAILED"
    } else {
      univ <- empty_univ(locus_row, "PHENOTYPE_DROPPED", "phenotype removed during locus processing", n_snps, K)
      matched <- match(result$phen, univ$phen)
      univ$h2.obs[matched] <- result$h2.obs
      if ("h2.latent" %in% names(result)) univ$h2.latent[matched] <- result$h2.latent
      if ("ascertained" %in% names(result)) univ$ascertained[matched] <- result$ascertained
      univ$p[matched] <- result$p
      univ$analysis_status[matched] <- "TESTED"
      univ$error[matched] <- NA_character_
      eligible <- univ$phen[univ$analysis_status == "TESTED" & is.finite(univ$p) & univ$p <= num("univariate_p_threshold")]
      eligible_pairs <- pairs[pairs$trait1 %in% eligible & pairs$trait2 %in% eligible, , drop = FALSE]
      rows <- vector("list", nrow(eligible_pairs))
      if (nrow(eligible_pairs)) for (j in seq_len(nrow(eligible_pairs))) {
        pair <- eligible_pairs[j, , drop = FALSE]
        set.seed(integer("random_seed") + index * 1000L + j)
        local <- tryCatch(run.bivar(locus, phenos = c(pair$trait1, pair$trait2)),
          error = function(error) structure(list(message = conditionMessage(error)), class = "bivar_error"))
        base <- data.frame(LOC = locus_row$LOC, CHR = locus_row$CHR, START = locus_row$START, STOP = locus_row$STOP,
          pair_id = pair$pair_id, trait1 = pair$trait1, trait2 = pair$trait2,
          discovery_rg = as.numeric(pair$discovery_rg), discovery_SE = as.numeric(pair$discovery_SE),
          discovery_P = as.numeric(pair$discovery_P), discovery_FDR = as.numeric(pair$discovery_FDR), stringsAsFactors = FALSE)
        covariance <- finite_or_na(locus$omega[pair$trait1, pair$trait2])
        unique_result <- !inherits(local, "bivar_error") && !is.null(local) && nrow(local) == 1L
        valid_result <- unique_result && is.finite(covariance) && all(is.finite(unlist(local[1L, c("rho", "rho.lower", "rho.upper", "r2", "r2.lower", "r2.upper", "p")])))
        if (!valid_result) {
          message_text <- if (inherits(local, "bivar_error")) local$message else if (!unique_result) "run.bivar returned no unique result" else "run.bivar returned an incomplete covariance, estimate, confidence interval, r2, or p-value"
          rows[[j]] <- cbind(base, local_covariance = covariance,
            rho = if (unique_result) finite_or_na(local$rho) else NA_real_,
            rho.lower = if (unique_result) finite_or_na(local$rho.lower) else NA_real_,
            rho.upper = if (unique_result) finite_or_na(local$rho.upper) else NA_real_,
            r2 = if (unique_result) finite_or_na(local$r2) else NA_real_,
            r2.lower = if (unique_result) finite_or_na(local$r2.lower) else NA_real_,
            r2.upper = if (unique_result) finite_or_na(local$r2.upper) else NA_real_,
            p = if (unique_result) finite_or_na(local$p) else NA_real_,
            analysis_status = "BIVARIATE_FAILED", error = message_text)
        } else {
          rows[[j]] <- cbind(base, local_covariance = covariance,
            rho = local$rho, rho.lower = local$rho.lower, rho.upper = local$rho.upper,
            r2 = local$r2, r2.lower = local$r2.lower, r2.upper = local$r2.upper, p = local$p,
            analysis_status = "TESTED", error = NA_character_)
        }
      }
      bivar <- if (length(rows)) do.call(rbind, rows) else empty_bivar()
      status_value <- "PROCESSED"
    }
  }
  status_row <- data.frame(LOC = locus_row$LOC, CHR = locus_row$CHR, START = locus_row$START, STOP = locus_row$STOP,
    status = status_value, n_snps = n_snps, K = K, univariate_tested = sum(univ$analysis_status == "TESTED"),
    eligible_bivariate_pairs = nrow(bivar), elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    stringsAsFactors = FALSE)
  list(fingerprint = fingerprint, phase = "discovery", locus_index = as.integer(index),
    LOC = as.character(locus_row$LOC), CHR = as.integer(locus_row$CHR), status = status_row, univ = univ, bivar = bivar)
}

if (values$phase == "discovery") {
  index <- values$locus_index
  path <- checkpoint_path("discovery", index)
  if (file.exists(path)) stop("Supervisor invoked a discovery worker whose ready bundle already exists")
  input <- scoped_input(loci[index, , drop = FALSE])
  output <- process_discovery_locus(index, input)
  publish_rds_no_replace(output, values$worker_output)
  marker <- list(locus = output$LOC, chromosome = output$CHR, n_snps = output$status$n_snps,
    pair = if (nrow(output$bivar)) paste(output$bivar$pair_id, collapse = ";") else "NONE", qc = output$status$status)
  rm(input, output)
  invisible(gc(full = TRUE))
  verify_current_fingerprint()
  worker_marker("discovery", index, marker$locus, marker$chromosome, marker$n_snps, marker$pair, marker$qc,
    if (marker$qc == "PROCESS_FAILED") "FALSE" else "TRUE", values$worker_output)
  quit(save = "no", status = 0L)
}

discovery_family_path <- file.path(values$checkpoint_dir, "aggregate-discovery", "unit_all", "result.rds")
if (values$phase == "conditional") {
  if (!file.exists(discovery_family_path)) stop("Discovery-family BH checkpoint is missing; run aggregate-discovery first")
  sealed_discovery <- readRDS(discovery_family_path)
  if (!identical(sealed_discovery$fingerprint, fingerprint) || !identical(sealed_discovery$phase, "aggregate-discovery")) stop("Stale discovery-family BH checkpoint")
  status <- sealed_discovery$status
  univ <- sealed_discovery$univ
  bivar <- sealed_discovery$bivar
} else {
  discovery_paths <- file.path(values$checkpoint_dir, "discovery", sprintf("locus_%04d", seq_len(nrow(loci))), "result.rds")
  if (!all(file.exists(discovery_paths))) stop(sprintf("Incomplete discovery checkpoint family: %d/%d", sum(file.exists(discovery_paths)), nrow(loci)))
  discovery_results <- lapply(seq_along(discovery_paths), function(index) validate_checkpoint(readRDS(discovery_paths[[index]]), "discovery", index))
  status <- do.call(rbind, lapply(discovery_results, `[[`, "status"))
  univ <- do.call(rbind, lapply(discovery_results, `[[`, "univ"))
  bivar <- do.call(rbind, lapply(discovery_results, `[[`, "bivar"))
  if (nrow(status) != integer("expected_loci") || nrow(univ) != integer("planned_univariate_tests")) stop("Incomplete discovery checkpoint contents")
  univ$univariate_test_family_n <- integer("planned_univariate_tests")
  univ$univariate_p_threshold <- num("univariate_p_threshold")
  univ$p_bonferroni <- pmin(1, univ$p * integer("planned_univariate_tests"))
  univ$p_fdr <- NA_real_
  tested_univ <- which(univ$analysis_status == "TESTED" & is.finite(univ$p))
  univ_family_p <- rep.int(1, nrow(univ))
  univ_family_p[tested_univ] <- univ$p[tested_univ]
  univ$p_fdr[tested_univ] <- p.adjust(univ_family_p, method = "BH")[tested_univ]
  bivar$bivariate_test_family_n <- rep.int(nrow(bivar), nrow(bivar))
  bivar$p_fdr <- rep(NA_real_, nrow(bivar))
  tested_bivar <- which(bivar$analysis_status == "TESTED" & is.finite(bivar$p))
  bivar_family_p <- rep.int(1, nrow(bivar))
  bivar_family_p[tested_bivar] <- bivar$p[tested_bivar]
  bivar$p_fdr[tested_bivar] <- p.adjust(bivar_family_p, method = "BH")[tested_bivar]
  bivar$fdr_significant <- !is.na(bivar$p_fdr) & bivar$p_fdr <= num("bivariate_fdr_alpha")
  status$analysis_fingerprint <- fingerprint
  status$lava_version <- as.character(packageVersion("LAVA"))
  status$reference_prefix <- reference_prefix
  discovery_family <- list(fingerprint = fingerprint, phase = "aggregate-discovery", status = status, univ = univ, bivar = bivar)
}

if (values$phase == "aggregate-discovery") {
  if (file.exists(discovery_family_path)) stop("Supervisor invoked aggregation after its ready bundle already existed")
  candidate_rows <- list()
  candidate_index <- 0L
  for (index in seq_len(nrow(loci))) {
    locus_row <- loci[index, , drop = FALSE]
    targets <- bivar[bivar$LOC == locus_row$LOC & bivar$fdr_significant & bivar$pair_id %in% c("A", "B"), , drop = FALSE]
    eligible_models <- 0L
    eligible_pairs <- character(0)
    if (nrow(targets)) for (target_index in seq_len(nrow(targets))) {
      target <- targets[target_index, , drop = FALSE]
      models <- conditioners[conditioners$pair_id == target$pair_id & conditioners$conditional_model_id != "NONE", , drop = FALSE]
      if (nrow(models)) for (model_index in seq_len(nrow(models))) {
        model <- models[model_index, , drop = FALSE]
        covars <- strsplit(model$covariates, ";", fixed = TRUE)[[1L]]
        local_univ <- univ[univ$LOC == locus_row$LOC & univ$phen %in% covars, , drop = FALSE]
        eligible <- length(covars) > 0L && nrow(local_univ) == length(covars) &&
          all(local_univ$analysis_status == "TESTED" & is.finite(local_univ$p) & local_univ$p <= num("univariate_p_threshold"))
        if (eligible) {
          eligible_models <- eligible_models + 1L
          eligible_pairs <- unique(c(eligible_pairs, as.character(target$pair_id)))
        }
      }
    }
    if (eligible_models > 0L) {
      candidate_index <- candidate_index + 1L
      candidate_rows[[candidate_index]] <- data.frame(
        locus_index = index, locus = locus_row$LOC, chromosome = locus_row$CHR,
        start = locus_row$START, stop = locus_row$STOP,
        n_snps = status$n_snps[match(locus_row$LOC, status$LOC)],
        eligible_models = eligible_models, eligible_pairs = paste(eligible_pairs, collapse = ";"),
        stringsAsFactors = FALSE)
    }
  }
  candidates <- if (length(candidate_rows)) do.call(rbind, candidate_rows) else data.frame(
    locus_index = integer(0), locus = integer(0), chromosome = integer(0),
    start = integer(0), stop = integer(0), n_snps = integer(0),
    eligible_models = integer(0), eligible_pairs = character(0), stringsAsFactors = FALSE)
  candidate_path <- file.path(dirname(values$worker_output), "conditional_candidates.tsv")
  write_table_no_replace(candidates, candidate_path)
  publish_rds_no_replace(discovery_family, values$worker_output)
  total_snps <- sum(status$n_snps, na.rm = TRUE)
  rm(discovery_results, discovery_family, status, univ, bivar, candidate_rows, candidates)
  invisible(gc(full = TRUE))
  verify_current_fingerprint()
  worker_marker("aggregate-discovery", 0L, "ALL", "ALL", total_snps, "A;B;CONTROL", "FULL_FAMILY_BH_COMPLETE", "TRUE", values$worker_output)
  quit(save = "no", status = 0L)
}

if (values$phase != "conditional") {
  if (!file.exists(discovery_family_path)) stop("Discovery-family BH checkpoint is missing; run aggregate-discovery first")
  sealed_discovery <- readRDS(discovery_family_path)
  if (!identical(sealed_discovery$fingerprint, fingerprint) || !identical(sealed_discovery$phase, "aggregate-discovery") ||
      !identical(sealed_discovery$status, status) || !identical(sealed_discovery$univ, univ) || !identical(sealed_discovery$bivar, bivar)) {
    stop("Discovery-family BH checkpoint differs from the complete immutable locus family")
  }
}

process_conditional_locus <- function(index, family) {
  locus_row <- loci[index, , drop = FALSE]
  targets <- family$bivar[family$bivar$LOC == locus_row$LOC & family$bivar$fdr_significant & family$bivar$pair_id %in% c("A", "B"), , drop = FALSE]
  conditional_rows <- list()
  conditional_index <- 0L
  locus <- NULL
  input <- NULL
  locus_process_attempted <- FALSE
  locus_process_ok <- FALSE
  if (nrow(targets)) for (j in seq_len(nrow(targets))) {
    target <- targets[j, , drop = FALSE]
    models <- conditioners[conditioners$pair_id == target$pair_id & conditioners$conditional_model_id != "NONE", , drop = FALSE]
    if (!nrow(models)) stop(paste("No frozen conditional model for pair", target$pair_id))
    for (model_index in seq_len(nrow(models))) {
      conditional_index <- conditional_index + 1L
      model <- models[model_index, , drop = FALSE]
      covars <- strsplit(model$covariates, ";", fixed = TRUE)[[1L]]
      base <- data.frame(LOC = target$LOC, CHR = target$CHR, START = target$START, STOP = target$STOP,
        pair_id = target$pair_id, trait1 = target$trait1, trait2 = target$trait2,
        conditional_model_id = model$conditional_model_id, covariates = model$covariates, stringsAsFactors = FALSE)
      local_univ <- family$univ[family$univ$LOC == target$LOC & family$univ$phen %in% covars, , drop = FALSE]
      eligible <- length(covars) > 0L && nrow(local_univ) == length(covars) &&
        all(local_univ$analysis_status == "TESTED" & is.finite(local_univ$p) & local_univ$p <= num("univariate_p_threshold"))
      if (!eligible) {
        conditional_rows[[conditional_index]] <- cbind(base, pcor = NA_real_, ci.lower = NA_real_, ci.upper = NA_real_, p = NA_real_,
          r2.trait1_z = NA_real_, r2.trait2_z = NA_real_, analysis_status = "CONDITIONER_LOCAL_H2_INELIGIBLE",
          error = "At least one predeclared conditioner did not pass the same-locus univariate gate")
        next
      }
      if (is.null(input)) input <- scoped_input(locus_row)
      if (is.null(locus)) {
        locus_process_attempted <- TRUE
        locus <- tryCatch(process.locus(locus_row, input, min.K = integer("min_K"),
          prune.thresh = num("prune_threshold"), max.prop.K = num("max_proportion_K"),
          drop.failed = TRUE, max.block.size = integer("max_block_size"), cap.estimates = logical_value("cap_estimates")),
          error = function(error) structure(list(message = conditionMessage(error)), class = "conditional_error"))
        locus_process_ok <- !inherits(locus, "conditional_error") && !is.null(locus)
      }
      model_seed <- match(model$conditional_model_id, c("A_BMI_ONLY", "A_SLEEP_APNEA_ONLY", "B_MDD_ONLY"))
      set.seed(integer("random_seed") + index * 1000L + 500L + model_seed)
      partial <- if (inherits(locus, "conditional_error") || is.null(locus)) locus else tryCatch(
        run.pcor(locus, target = c(target$trait1, target$trait2), phenos = c(target$trait1, target$trait2, covars), max.r2 = num("conditional_max_r2")),
        error = function(error) structure(list(message = conditionMessage(error)), class = "conditional_error"))
      unique_partial <- !inherits(partial, "conditional_error") && !is.null(partial) && nrow(partial) == 1L
      partial_values <- if (unique_partial) unlist(partial[1L, c("pcor", "ci.lower", "ci.upper", "p", "r2.phen1_z", "r2.phen2_z")]) else rep(NA_real_, 6L)
      valid_partial <- unique_partial && all(is.finite(partial_values)) && partial[["r2.phen1_z"]] < num("conditional_max_r2") && partial[["r2.phen2_z"]] < num("conditional_max_r2")
      unstable_max_r2 <- unique_partial && is.finite(partial[["pcor"]]) && is.finite(partial[["ci.lower"]]) && is.finite(partial[["ci.upper"]]) &&
        is.finite(partial[["r2.phen1_z"]]) && is.finite(partial[["r2.phen2_z"]]) &&
        (partial[["r2.phen1_z"]] >= num("conditional_max_r2") || partial[["r2.phen2_z"]] >= num("conditional_max_r2"))
      if (!valid_partial) {
        message_text <- if (inherits(partial, "conditional_error")) partial$message else if (!unique_partial) "run.pcor returned no unique result" else if (unstable_max_r2) "run.pcor exceeded the predeclared conditional max-r2 stability threshold" else "run.pcor returned an incomplete partial estimate, confidence interval, r2, or p-value"
        terminal_status <- if (unstable_max_r2) "CONDITIONAL_UNSTABLE_MAX_R2" else "CONDITIONAL_FAILED"
        conditional_rows[[conditional_index]] <- cbind(base,
          pcor = if (unique_partial) finite_or_na(partial$pcor) else NA_real_,
          ci.lower = if (unique_partial) finite_or_na(partial$ci.lower) else NA_real_,
          ci.upper = if (unique_partial) finite_or_na(partial$ci.upper) else NA_real_,
          p = if (unique_partial && !unstable_max_r2) finite_or_na(partial$p) else NA_real_,
          r2.trait1_z = if (unique_partial) finite_or_na(partial[["r2.phen1_z"]]) else NA_real_,
          r2.trait2_z = if (unique_partial) finite_or_na(partial[["r2.phen2_z"]]) else NA_real_,
          analysis_status = terminal_status, error = message_text)
      } else {
        conditional_rows[[conditional_index]] <- cbind(base, pcor = partial$pcor, ci.lower = partial$ci.lower, ci.upper = partial$ci.upper,
          p = partial$p, r2.trait1_z = partial[["r2.phen1_z"]], r2.trait2_z = partial[["r2.phen2_z"]],
          analysis_status = "TESTED", error = NA_character_)
      }
    }
  }
  conditional <- if (length(conditional_rows)) do.call(rbind, conditional_rows) else empty_conditional()
  state <- if (!nrow(targets)) "NOT_BIVARIATE_FDR_ELIGIBLE" else if (!nrow(conditional)) "NO_FROZEN_MODEL" else if (!locus_process_attempted) "CONDITIONER_LOCAL_H2_INELIGIBLE" else if (locus_process_ok) "CONDITIONAL_LOCUS_PROCESSED" else "CONDITIONAL_LOCUS_PROCESS_FAILED"
  n_snps <- family$status$n_snps[match(locus_row$LOC, family$status$LOC)]
  if (!is.null(locus)) rm(locus)
  if (!is.null(input)) rm(input)
  invisible(gc(full = TRUE))
  list(fingerprint = fingerprint, phase = "conditional", locus_index = as.integer(index),
    LOC = as.character(locus_row$LOC), CHR = as.integer(locus_row$CHR), n_snps = as.integer(n_snps),
    pair = if (nrow(targets)) paste(targets$pair_id, collapse = ";") else "NONE", state = state, conditional = conditional)
}

if (values$phase == "conditional") {
  index <- values$locus_index
  path <- checkpoint_path("conditional", index)
  if (file.exists(path)) stop("Supervisor invoked a conditional worker whose ready bundle already exists")
  output <- process_conditional_locus(index, sealed_discovery)
  publish_rds_no_replace(output, values$worker_output)
  marker <- list(locus = output$LOC, chromosome = output$CHR, n_snps = output$n_snps, pair = output$pair, qc = output$state)
  rm(output, sealed_discovery, status, univ, bivar)
  invisible(gc(full = TRUE))
  verify_current_fingerprint()
  worker_marker("conditional", index, marker$locus, marker$chromosome, marker$n_snps, marker$pair, marker$qc,
    if (marker$qc == "CONDITIONAL_LOCUS_PROCESSED") "TRUE" else "FALSE", values$worker_output)
  quit(save = "no", status = 0L)
}

conditional_paths <- file.path(values$checkpoint_dir, "conditional", sprintf("locus_%04d", seq_len(nrow(loci))), "result.rds")
if (!all(file.exists(conditional_paths))) stop(sprintf("Incomplete conditional checkpoint family: %d/%d", sum(file.exists(conditional_paths)), nrow(loci)))
conditional_results <- lapply(seq_along(conditional_paths), function(index) validate_checkpoint(readRDS(conditional_paths[[index]]), "conditional", index))
conditional_parts <- lapply(conditional_results, `[[`, "conditional")
conditional <- if (any(vapply(conditional_parts, nrow, integer(1)) > 0L)) do.call(rbind, conditional_parts) else empty_conditional()
eligible_cond <- which(conditional$analysis_status != "CONDITIONER_LOCAL_H2_INELIGIBLE")
conditional$conditional_test_family_n <- rep.int(length(eligible_cond), nrow(conditional))
conditional$p_fdr <- rep(NA_real_, nrow(conditional))
tested_cond <- which(conditional$analysis_status == "TESTED" & is.finite(conditional$p))
conditional_family_p <- rep.int(1, length(eligible_cond))
names(conditional_family_p) <- eligible_cond
conditional_family_p[as.character(tested_cond)] <- conditional$p[tested_cond]
conditional_adjusted <- p.adjust(conditional_family_p, method = "BH")
conditional$p_fdr[tested_cond] <- conditional_adjusted[as.character(tested_cond)]
conditional$fdr_significant <- !is.na(conditional$p_fdr) & conditional$p_fdr <= num("conditional_fdr_alpha")

attempt_dir <- dirname(values$worker_output)
write_table_no_replace(status, file.path(attempt_dir, "lava_locus_status.tsv"))
write_table_no_replace(univ, file.path(attempt_dir, "lava_univariate.tsv"))
write_table_no_replace(bivar, file.path(attempt_dir, "lava_bivariate.tsv"))
write_table_no_replace(conditional, file.path(attempt_dir, "lava_conditional.tsv"))
validation <- system2("python3", c("scripts/121_validate_track_b_lava.py", "--validate-staged-results",
  "--staging-dir", attempt_dir, "--fingerprint", fingerprint), stdout = TRUE, stderr = TRUE)
if (!is.null(attr(validation, "status")) && attr(validation, "status") != 0L) stop(paste("Track B LAVA validation failed:", paste(validation, collapse = " | ")))
if (sum(grepl("^TRACK_B_LAVA_STAGED_RESULTS_SEMANTICALLY_VALIDATED", validation)) != 1L) stop("Staged result validator did not emit one PASS marker")
staging_validation <- file.path(attempt_dir, "staging_validation.txt")
temporary_validation <- tempfile(pattern = ".staging-validation.", tmpdir = attempt_dir)
writeLines(validation, temporary_validation, useBytes = TRUE)
if (!file.link(temporary_validation, staging_validation)) {
  unlink(temporary_validation)
  stop("Could not atomically publish staged validation attestation")
}
unlink(temporary_validation)
total_snps <- sum(status$n_snps, na.rm = TRUE)
final_payload <- list(fingerprint = fingerprint, phase = "finalize", locus_index = 0L,
  LOC = "ALL", CHR = "ALL", total_snps = total_snps,
  state = "STAGED_SEMANTIC_VALIDATION_COMPLETE")
publish_rds_no_replace(final_payload, values$worker_output)
rm(discovery_results, discovery_family, sealed_discovery, conditional_results, conditional_parts, status, univ, bivar, conditional)
invisible(gc(full = TRUE))
verify_current_fingerprint()
worker_marker("finalize", 0L, "ALL", "ALL", total_snps, "A;B;CONTROL", "STAGED_SEMANTIC_VALIDATION_COMPLETE", "TRUE", values$worker_output)
