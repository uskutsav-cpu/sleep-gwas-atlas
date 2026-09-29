#!/usr/bin/env Rscript
# One process per frozen locus; evaluate all five locked pairs at that locus.
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(LAVA))

fail <- function(message) stop(message, call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) fail("Usage: Rscript run_lava_family_locus.R CONFIG.json")
cfg <- jsonlite::fromJSON(args[[1]], simplifyVector = FALSE)
need <- function(test, message) if (!isTRUE(test)) fail(message)
need(as.character(packageVersion("LAVA")) == "0.1.5", "LAVA version drift")
need(isTRUE(cfg$runtime_validator_passed), "Coordinator did not pass the pinned LAVA runtime check")
need(grepl("^[0-9a-f]{64}$", cfg$runtime_validation_output_sha256),
     "Pinned LAVA runtime validation receipt is missing")
need(cfg$schema_version == 1L, "Unsupported locus-worker config")
need(cfg$analysis_id == "brain6-lava-local-rg-v2", "Unexpected analysis identity")
need(length(cfg$pairs) == 5L && length(unique(vapply(cfg$pairs, `[[`, "", "pair_id"))) == 5L,
     "Worker must receive the exact five-pair family")
need(length(cfg$locus_id) == 1L && !is.na(cfg$locus_id), "Worker must receive exactly one locus")
need(file.exists(cfg$loci_file), "Pinned locus definition is missing")
need(dir.exists(cfg$output_dir), "Worker output directory is missing")
set.seed(as.integer(cfg$random_seed))
loci <- LAVA::read.loci(cfg$loci_file)
locus_index <- match(as.character(cfg$locus_id), as.character(loci$LOC))
need(!is.na(locus_index) && sum(as.character(loci$LOC) == as.character(cfg$locus_id)) == 1L,
     "Requested locus is missing or duplicated")
locus_row <- loci[locus_index, , drop = FALSE]
policy <- cfg$execution_policy
pair_results <- vector("list", length(cfg$pairs))
univariate_results <- vector("list", length(cfg$pairs))

run_pair <- function(pair) {
  pair_id <- as.character(pair$pair_id)
  phenotypes <- as.character(unlist(pair$phenotypes, use.names = FALSE))
  tryCatch({
    need(length(phenotypes) == 2L && length(unique(phenotypes)) == 2L,
         paste("Invalid phenotype pair:", pair_id))
    need(file.exists(pair$input_info) && file.exists(pair$sample_overlap_file),
         paste("Pair inputs missing:", pair_id))
    if (as.integer(pair$shared_reference_variants) < as.integer(policy$locus_processing$min_K)) {
      return(list(result = data.table(pair_id = pair_id, locus_id = cfg$locus_id,
        status = "NO_OVERLAP", p = NA_real_, local_rg = NA_real_,
        reason = "FEWER_THAN_MIN_K_SHARED_VARIANTS_IN_SEALED_REFERENCE"),
        univariate = data.table()))
    }
    input <- LAVA::process.input(
      input.info.file = pair$input_info,
      sample.overlap.file = pair$sample_overlap_file,
      ref.prefix = pair$reference_prefix,
      phenos = phenotypes,
      input.dir = dirname(pair$input_info)
    )
    locus_messages <- capture.output(loc <- LAVA::process.locus(
      locus_row, input, phenos = phenotypes,
      min.K = as.integer(policy$locus_processing$min_K),
      prune.thresh = as.numeric(policy$locus_processing$prune_thresh),
      max.prop.K = as.numeric(policy$locus_processing$max_prop_K),
      drop.failed = isTRUE(policy$locus_processing$drop_failed),
      max.block.size = as.integer(policy$locus_processing$max_block_size),
      cap.estimates = isTRUE(policy$locus_processing$cap_estimates)
    ))
    if (is.null(loc)) {
      no_overlap <- any(grepl("Fewer than", locus_messages, fixed = TRUE))
      low_h2 <- any(grepl("Negative variance estimate for all phenotypes", locus_messages, fixed = TRUE))
      status <- if (no_overlap) "NO_OVERLAP" else if (low_h2) "UNIVARIATE_UNDERPOWERED" else "FAILED"
      reason <- if (no_overlap) "LAVA_REFERENCE_HARMONIZATION_LEFT_FEWER_THAN_MIN_K_SNPs" else
        if (low_h2) "LOW_LOCAL_H2_NOT_EXECUTION_FAILURE" else "LAVA_PROCESS_LOCUS_RETURNED_NULL"
      return(list(result = data.table(pair_id = pair_id, locus_id = cfg$locus_id,
        status = status, p = NA_real_, local_rg = NA_real_, reason = reason),
        univariate = data.table()))
    }
    available_phenotypes <- intersect(phenotypes, as.character(loc$phenos))
    if (length(available_phenotypes) < 2L) {
      u <- as.data.table(LAVA::run.univ(loc, phenos = available_phenotypes,
        cap.estimates = isTRUE(policy$univariate$cap_estimates)))
      if (nrow(u) > 0L) u[, `:=`(pair_id = pair_id, locus_id = cfg$locus_id)]
      return(list(result = data.table(pair_id = pair_id, locus_id = cfg$locus_id,
        status = "UNIVARIATE_UNDERPOWERED", p = NA_real_, local_rg = NA_real_,
        reason = "ONE_TRAIT_DROPPED_FOR_LOW_LOCAL_H2"), univariate = u))
    }
    u <- as.data.table(LAVA::run.univ(loc, phenos = phenotypes,
      cap.estimates = isTRUE(policy$univariate$cap_estimates)))
    need(nrow(u) == 2L && all(phenotypes %in% u$phen) && all(is.finite(u$p)) &&
           all(u$p >= 0 & u$p <= 1),
         paste("Invalid local univariate tests:", pair_id))
    u[, `:=`(pair_id = pair_id, locus_id = cfg$locus_id)]
    if (!all(u$p < policy$univariate$gate_p_strictly_less_than)) {
      return(list(result = data.table(pair_id = pair_id, locus_id = cfg$locus_id,
        status = "UNIVARIATE_UNDERPOWERED", p = NA_real_, local_rg = NA_real_,
        reason = "ONE_OR_MORE_TRAITS_DID_NOT_PASS_FROZEN_LOCAL_UNIVARIATE_GATE"),
        univariate = u))
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
    list(result = data.table(pair_id = pair_id, locus_id = cfg$locus_id,
      status = "TESTED", p = b$p[[1L]], local_rg = b$rho[[1L]], reason = ""),
      univariate = u)
  }, error = function(e) {
    list(result = data.table(pair_id = pair_id, locus_id = cfg$locus_id,
      status = "FAILED", p = NA_real_, local_rg = NA_real_,
      reason = substr(conditionMessage(e), 1L, 500L)), univariate = data.table())
  })
}

for (index in seq_along(cfg$pairs)) {
  value <- run_pair(cfg$pairs[[index]])
  pair_results[[index]] <- value$result
  univariate_results[[index]] <- value$univariate
}
results <- rbindlist(pair_results, use.names = TRUE, fill = TRUE)
univariate <- rbindlist(univariate_results, use.names = TRUE, fill = TRUE)
need(nrow(results) == 5L && setequal(results$pair_id,
  vapply(cfg$pairs, `[[`, "", "pair_id")), "Worker did not preserve all five pair slots")
fwrite(results, file.path(cfg$output_dir, "pair_results.tsv"), sep = "\t", na = "NA")
fwrite(univariate, file.path(cfg$output_dir, "univariate.tsv"), sep = "\t", na = "NA")
write_json(list(status = "PASS", analysis_id = cfg$analysis_id,
  family_lock_sha256 = cfg$family_lock_sha256,
  execution_lock_sha256 = cfg$execution_lock_sha256,
  locus_id = cfg$locus_id, n_pairs = nrow(results),
  n_tested = sum(results$status == "TESTED"),
  n_underpowered = sum(results$status == "UNIVARIATE_UNDERPOWERED"),
  n_failed = sum(results$status == "FAILED"),
  participant_level_overlap_verified = FALSE,
  lava_version = as.character(packageVersion("LAVA"))),
  file.path(cfg$output_dir, "status.json"), auto_unbox = TRUE, pretty = TRUE,
  na = "null", digits = NA)
cat(sprintf("BRAIN6_LAVA_LOCUS_COMPLETE locus=%s tested=%d underpowered=%d failed=%d\n",
  cfg$locus_id, sum(results$status == "TESTED"),
  sum(results$status == "UNIVARIATE_UNDERPOWERED"), sum(results$status == "FAILED")))
