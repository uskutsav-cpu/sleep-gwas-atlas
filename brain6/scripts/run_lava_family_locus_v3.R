#!/usr/bin/env Rscript
# One process per frozen locus. Compute unique trait-only h2 tests, then gate
# the five pair-specific bivariate tests with those shared trait-locus results.
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(LAVA))

fail <- function(message) stop(message, call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) fail("Usage: Rscript run_lava_family_locus_v3.R CONFIG.json")
cfg <- jsonlite::fromJSON(args[[1L]], simplifyVector = FALSE)
need <- function(test, message) if (!isTRUE(test)) fail(message)
need(as.character(packageVersion("LAVA")) == "0.1.5", "LAVA version drift")
need(isTRUE(cfg$runtime_validator_passed), "Coordinator did not pass the pinned LAVA runtime check")
need(grepl("^[0-9a-f]{64}$", cfg$runtime_validation_output_sha256),
     "Pinned LAVA runtime validation receipt is missing")
need(cfg$schema_version == 1L && cfg$analysis_id == "brain6-lava-local-rg-v3",
     "Unexpected locus-worker identity")
need(length(cfg$pairs) == 5L && length(unique(vapply(cfg$pairs, `[[`, "", "pair_id"))) == 5L,
     "Worker must receive the exact five-pair family")
need(length(cfg$trait_ids) == 7L && length(unique(cfg$trait_ids)) == 7L,
     "Worker must receive the exact seven-trait univariate family")
need(length(cfg$locus_id) == 1L && !is.na(cfg$locus_id), "Worker must receive exactly one locus")
need(file.exists(cfg$loci_file) && dir.exists(cfg$output_dir), "Locus input/output path is missing")

# Numerical-stability patch for the pinned LAVA 0.1.5 blockwise LD path.
# Symmetrize only dimension-scaled floating-point roundoff; fail on material asymmetry.
install_lava_roundoff_patch <- function() {
  ns <- asNamespace("LAVA")
  original <- get("decompose.ld", envir = ns, inherits = FALSE)
  body_lines <- deparse(body(original), width.cutoff = 500L)
  hit <- which(trimws(body_lines) == "M = t(R.base) %*% ld$ld %*% R.base")
  need(length(hit) == 1L, "Pinned LAVA blockwise decomposition source changed; refusing patch")
  insert <- c(
    "            asymmetry <- max(abs(M - t(M)))",
    "            tolerance <- 100 * .Machine$double.eps * max(1, nrow(M)) * max(1, max(abs(M)))",
    "            if (!is.finite(asymmetry) || asymmetry > tolerance) stop('Material asymmetry in block-reduced LD matrix; refusing numerical symmetrization', call. = FALSE)",
    "            if (asymmetry > 0) {",
    "                cat(sprintf('BRAIN6_LAVA_SYMMETRIZED_ROUNDOFF n=%d asymmetry=%.17g tolerance=%.17g\\n', nrow(M), asymmetry, tolerance))",
    "                M <- (M + t(M)) / 2",
    "            }"
  )
  body_lines <- append(body_lines, insert, after = hit)
  patched <- original
  body(patched) <- parse(text = paste(body_lines, collapse = "\n"))[[1L]]
  environment(patched) <- ns
  unlockBinding("decompose.ld", ns)
  assign("decompose.ld", patched, envir = ns)
  lockBinding("decompose.ld", ns)
  cat("BRAIN6_LAVA_NUMERICAL_PATCH id=LAVA015_BLOCK_REDUCED_SYMMETRY_V1 tolerance=100*eps*n*scale\n")
}
install_lava_roundoff_patch()

set.seed(as.integer(cfg$random_seed))
loci <- LAVA::read.loci(cfg$loci_file)
locus_index <- match(as.character(cfg$locus_id), as.character(loci$LOC))
need(!is.na(locus_index) && sum(as.character(loci$LOC) == as.character(cfg$locus_id)) == 1L,
     "Requested locus is missing or duplicated")
locus_row <- loci[locus_index, , drop = FALSE]
policy <- cfg$execution_policy
template_pair <- cfg$pairs[[1L]]
trait_ids <- as.character(unlist(cfg$trait_ids, use.names = FALSE))
# The family runner verifies the source locus receipt before creating this
# config. Passing its seven canonical rows avoids recomputing unique univariate
# tests in the second, pairwise stage.
need(length(cfg$canonical_univariate_rows) == 7L,
     "Pairwise stage requires all seven verified canonical univariate rows")
univariate <- rbindlist(lapply(cfg$canonical_univariate_rows, function(row) {
  as_num <- function(value) if (is.null(value) || length(value) == 0L || is.na(value)) NA_real_ else as.numeric(value)
  data.table(phen = as.character(row$phen), locus_id = as.character(row$locus_id),
    status = as.character(row$status), h2.obs = as_num(row$h2.obs),
    h2.latent = as_num(row$h2.latent), p = as_num(row$p), reason = as.character(row$reason))
}), use.names = TRUE, fill = TRUE)
need(nrow(univariate) == 7L && setequal(univariate$phen, trait_ids) &&
       all(univariate$locus_id == as.character(cfg$locus_id)) &&
       !anyDuplicated(univariate$phen) &&
       all(univariate$status %in% c("TESTED", "NOT_RUN", "FAILED")),
     "Canonical univariate rows are incomplete or have invalid identities/statuses")
tested_rows <- univariate[status == "TESTED"]
need(all(is.finite(tested_rows$p)) && all(tested_rows$p >= 0 & tested_rows$p <= 1),
     "Canonical tested p-values must be finite and in [0,1]")

# Defer reference and summary-statistic loading until a pair has passed the
# canonical gate. Underpowered loci require no reference load and no reread.
lava_ns <- asNamespace("LAVA")
check_reference <- get("check.reference", envir = lava_ns, inherits = FALSE)
load_reference <- get("load.reference", envir = lava_ns, inherits = FALSE)
read_sumstats <- get("read.sumstats.file", envir = lava_ns, inherits = FALSE)
harmonize_snps <- get("harmonize.snps", envir = lava_ns, inherits = FALSE)
align_sumstats <- get("align", envir = lava_ns, inherits = FALSE)
process_overlap <- get("process.sample.overlap", envir = lava_ns, inherits = FALSE)
cache <- new.env(parent = emptyenv())
cache$reference <- NULL
cache$input_info <- NULL
cache$sumstats <- list()

prepare_input <- function(phenotypes, sample_overlap_file = NULL) {
  if (is.null(cache$reference))
    cache$reference <- load_reference(check_reference(template_pair$reference_prefix))
  if (is.null(cache$input_info)) {
    input_dir <- dirname(template_pair$input_info)
    cache$input_info <- fread(template_pair$input_info)
    need(all(c("phenotype", "cases", "controls", "filename") %in% names(cache$input_info)) &&
           all(trait_ids %in% cache$input_info$phenotype),
         "Locus input info does not cover the frozen seven-trait family")
    cache$input_info$filename <- file.path(input_dir, cache$input_info$filename)
    cache$input_info$N <- cache$input_info$cases + cache$input_info$controls
    cache$input_info$prop_cases <- cache$input_info$cases / cache$input_info$N
    cache$input_info$binary <- !is.na(cache$input_info$prop_cases) & (cache$input_info$prop_cases != 1)
  }
  input <- new.env(parent = globalenv())
  input$info <- cache$input_info[match(phenotypes, cache$input_info$phenotype), , drop = FALSE]
  input$P <- length(phenotypes)
  input$reference <- cache$reference
  input$sample.overlap <- if (input$P > 1L && !is.null(sample_overlap_file))
    process_overlap(sample_overlap_file, phenotypes) else NULL
  input$sum.stats <- lapply(phenotypes, function(trait) {
    if (is.null(cache$sumstats[[trait]])) {
      path <- input$info[match(trait, input$info$phenotype), filename][[1L]]
      need(file.exists(path), paste("Trait sumstats are missing:", trait))
      cache$sumstats[[trait]] <- read_sumstats(path, trait)
    }
    data.table::copy(cache$sumstats[[trait]])
  })
  names(input$sum.stats) <- phenotypes
  harmonize_snps(input)
  align_sumstats(input)
  input
}

run_pair <- function(pair) {
  pair_id <- as.character(pair$pair_id)
  phenotypes <- as.character(unlist(pair$phenotypes, use.names = FALSE))
  tryCatch({
    need(length(phenotypes) == 2L && all(phenotypes %in% univariate$phen),
         paste("Invalid phenotype pair:", pair_id))
    if (as.integer(pair$shared_reference_variants) < as.integer(policy$locus_processing$min_K)) {
      return(data.table(pair_id = pair_id, locus_id = cfg$locus_id, status = "NO_OVERLAP",
        p = NA_real_, local_rg = NA_real_, reason = "FEWER_THAN_MIN_K_SHARED_VARIANTS_IN_SEALED_REFERENCE"))
    }
    rows <- univariate[match(phenotypes, univariate$phen)]
    if (any(rows$status != "TESTED")) {
      return(data.table(pair_id = pair_id, locus_id = cfg$locus_id,
        status = "UNIVARIATE_UNDERPOWERED", p = NA_real_, local_rg = NA_real_,
        reason = "ONE_OR_MORE_UNIQUE_TRAIT_LOCAL_TESTS_NOT_RUN"))
    }
    if (!all(rows$p < as.numeric(policy$univariate$gate_p_strictly_less_than))) {
      return(data.table(pair_id = pair_id, locus_id = cfg$locus_id,
        status = "UNIVARIATE_UNDERPOWERED", p = NA_real_, local_rg = NA_real_,
        reason = "ONE_OR_MORE_UNIQUE_TRAIT_LOCAL_TESTS_FAILED_FROZEN_GATE"))
    }
    need(file.exists(pair$input_info) && file.exists(pair$sample_overlap_file),
         paste("Pair inputs missing:", pair_id))
    input <- prepare_input(phenotypes, pair$sample_overlap_file)
    messages <- capture.output(loc <- LAVA::process.locus(
      locus_row, input, phenos = phenotypes,
      min.K = as.integer(policy$locus_processing$min_K),
      prune.thresh = as.numeric(policy$locus_processing$prune_thresh),
      max.prop.K = as.numeric(policy$locus_processing$max_prop_K),
      drop.failed = isTRUE(policy$locus_processing$drop_failed),
      max.block.size = as.integer(policy$locus_processing$max_block_size),
      cap.estimates = isTRUE(policy$locus_processing$cap_estimates)
    ))
    if (is.null(loc)) {
      no_overlap <- any(grepl("Fewer than", messages, fixed = TRUE))
      low_h2 <- any(grepl("Negative variance estimate for all phenotypes", messages, fixed = TRUE))
      status <- if (no_overlap) "NO_OVERLAP" else if (low_h2) "UNIVARIATE_UNDERPOWERED" else "FAILED"
      reason <- if (no_overlap) "LAVA_PAIR_HARMONIZATION_LEFT_FEWER_THAN_MIN_K" else
        if (low_h2) "LOW_PAIRWISE_LOCAL_H2_NOT_EXECUTION_FAILURE" else "LAVA_PROCESS_LOCUS_RETURNED_NULL"
      return(data.table(pair_id = pair_id, locus_id = cfg$locus_id,
        status = status, p = NA_real_, local_rg = NA_real_, reason = reason))
    }
    if (length(intersect(phenotypes, as.character(loc$phenos))) < 2L) {
      return(data.table(pair_id = pair_id, locus_id = cfg$locus_id,
        status = "UNIVARIATE_UNDERPOWERED", p = NA_real_, local_rg = NA_real_,
        reason = "PAIRWISE_LOCAL_H2_DROPPED_A_TRAIT"))
    }
    b <- as.data.table(LAVA::run.bivar(loc, phenos = phenotypes,
      adap.thresh = as.numeric(unlist(policy$bivariate$adap_thresh)),
      p.values = isTRUE(policy$bivariate$p_values),
      CIs = isTRUE(policy$bivariate$confidence_intervals),
      param.lim = as.numeric(policy$bivariate$parameter_limit),
      cap.estimates = isTRUE(policy$bivariate$cap_estimates)))
    need(nrow(b) == 1L && all(c("rho", "p") %in% names(b)) &&
           is.finite(b$rho[[1L]]) && abs(b$rho[[1L]]) <= 1 &&
           is.finite(b$p[[1L]]) && b$p[[1L]] >= 0 && b$p[[1L]] <= 1,
         paste("Invalid LAVA bivariate result:", pair_id))
    data.table(pair_id = pair_id, locus_id = cfg$locus_id, status = "TESTED",
      p = b$p[[1L]], local_rg = b$rho[[1L]], reason = "")
  }, error = function(e) {
    data.table(pair_id = pair_id, locus_id = cfg$locus_id, status = "FAILED",
      p = NA_real_, local_rg = NA_real_, reason = substr(conditionMessage(e), 1L, 500L))
  })
}

pair_results <- rbindlist(lapply(cfg$pairs, run_pair), use.names = TRUE, fill = TRUE)
need(nrow(pair_results) == 5L && setequal(pair_results$pair_id,
  vapply(cfg$pairs, `[[`, "", "pair_id")), "Worker did not preserve all five pair slots")
fwrite(pair_results, file.path(cfg$output_dir, "pair_results.tsv"), sep = "\t", na = "NA")
fwrite(univariate, file.path(cfg$output_dir, "univariate.tsv"), sep = "\t", na = "NA")
write_json(list(status = "PASS", analysis_id = cfg$analysis_id,
  family_lock_sha256 = cfg$family_lock_sha256,
  execution_lock_sha256 = cfg$execution_lock_sha256,
  bivariate_run_id = cfg$bivariate_run_id,
  canonical_run_id = cfg$canonical_run_id,
  canonical_cell_sha256 = cfg$canonical_cell_sha256,
  canonical_receipt_sha256 = cfg$canonical_receipt_sha256,
  locus_id = cfg$locus_id, n_pairs = nrow(pair_results),
  n_tested = sum(pair_results$status == "TESTED"),
  n_underpowered = sum(pair_results$status == "UNIVARIATE_UNDERPOWERED"),
  n_failed = sum(pair_results$status == "FAILED"),
  unique_univariate_tested = sum(univariate$status == "TESTED"),
  unique_univariate_not_run = sum(univariate$status != "TESTED"),
  participant_level_overlap_verified = FALSE,
  lava_version = as.character(packageVersion("LAVA"))),
  file.path(cfg$output_dir, "status.json"), auto_unbox = TRUE, pretty = TRUE,
  na = "null", digits = NA)
cat(sprintf("BRAIN6_LAVA_LOCUS_V3_COMPLETE locus=%s tested=%d underpowered=%d failed=%d unique_univariate_tested=%d/7\n",
  cfg$locus_id, sum(pair_results$status == "TESTED"),
  sum(pair_results$status == "UNIVARIATE_UNDERPOWERED"),
  sum(pair_results$status == "FAILED"), sum(univariate$status == "TESTED")))
