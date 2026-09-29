#!/usr/bin/env Rscript
# Diagnostic: calculate one LAVA local univariate test using a trait-only SNP set.
# This is deliberately separate from the frozen v2 family runner.
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(LAVA))

fail <- function(message) stop(message, call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !file.exists(args[[1L]])) {
  fail("Usage: Rscript diagnose_lava_single_trait_univariate.R CONFIG.json")
}
cfg <- jsonlite::fromJSON(args[[1L]], simplifyVector = FALSE)
need <- function(test, message) if (!isTRUE(test)) fail(message)
need(as.character(packageVersion("LAVA")) == "0.1.5", "LAVA version drift")
for (field in c("input_info", "input_dir", "reference_prefix", "loci_file",
                "locus_id", "trait_id", "output_file", "execution_policy")) {
  need(!is.null(cfg[[field]]), paste("Missing config field:", field))
}
need(file.exists(cfg$input_info), "LAVA input info file is missing")
need(dir.exists(cfg$input_dir), "LAVA input directory is missing")
need(file.exists(cfg$loci_file), "Pinned LAVA locus definition is missing")
need(length(cfg$locus_id) == 1L && !is.na(cfg$locus_id), "Exactly one locus is required")
need(length(cfg$trait_id) == 1L && nzchar(cfg$trait_id), "Exactly one trait is required")

loci <- LAVA::read.loci(cfg$loci_file)
locus_index <- match(as.character(cfg$locus_id), as.character(loci$LOC))
need(!is.na(locus_index) && sum(as.character(loci$LOC) == as.character(cfg$locus_id)) == 1L,
     "Requested locus is missing or duplicated")

policy <- cfg$execution_policy
result <- tryCatch({
  # No sample-overlap matrix is involved in a single-trait h2 test.
  input <- LAVA::process.input(
    input.info.file = cfg$input_info,
    sample.overlap.file = NULL,
    ref.prefix = cfg$reference_prefix,
    phenos = cfg$trait_id,
    input.dir = cfg$input_dir
  )
  captured <- capture.output(loc <- LAVA::process.locus(
    loci[locus_index, , drop = FALSE], input, phenos = cfg$trait_id,
    min.K = as.integer(policy$locus_processing$min_K),
    prune.thresh = as.numeric(policy$locus_processing$prune_thresh),
    max.prop.K = as.numeric(policy$locus_processing$max_prop_K),
    drop.failed = isTRUE(policy$locus_processing$drop_failed),
    max.block.size = as.integer(policy$locus_processing$max_block_size),
    cap.estimates = isTRUE(policy$locus_processing$cap_estimates)
  ))
  if (is.null(loc)) {
    reason <- if (any(grepl("Fewer than", captured, fixed = TRUE)))
      "FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS" else
        "LAVA_PROCESS_LOCUS_RETURNED_NULL"
    data.table(trait_id = cfg$trait_id, locus_id = as.character(cfg$locus_id),
      status = "NOT_RUN", n_snps = NA_integer_, n_components = NA_integer_,
      h2_obs = NA_real_, h2_latent = NA_real_, p = NA_real_, reason = reason)
  } else {
    u <- as.data.table(LAVA::run.univ(loc, phenos = cfg$trait_id,
      cap.estimates = isTRUE(policy$univariate$cap_estimates)))
    need(nrow(u) == 1L && identical(as.character(u$phen[[1L]]), cfg$trait_id),
         "LAVA did not return exactly one requested univariate row")
    need(is.finite(u$p[[1L]]) && u$p[[1L]] >= 0 && u$p[[1L]] <= 1,
         "LAVA returned an invalid univariate p-value")
    data.table(trait_id = cfg$trait_id, locus_id = as.character(cfg$locus_id),
      status = "TESTED", n_snps = as.integer(loc$n.snps),
      n_components = as.integer(loc$K), h2_obs = as.numeric(u$h2.obs[[1L]]),
      h2_latent = if ("h2.latent" %in% names(u)) as.numeric(u$h2.latent[[1L]]) else NA_real_,
      p = as.numeric(u$p[[1L]]), reason = "")
  }
}, error = function(e) {
  data.table(trait_id = cfg$trait_id, locus_id = as.character(cfg$locus_id),
    status = "FAILED", n_snps = NA_integer_, n_components = NA_integer_,
    h2_obs = NA_real_, h2_latent = NA_real_, p = NA_real_,
    reason = substr(conditionMessage(e), 1L, 500L))
})

need(nzchar(cfg$output_file), "Diagnostic output file path is empty")
if (file.exists(cfg$output_file)) fail("Refusing to overwrite an existing diagnostic output")
dir.create(dirname(cfg$output_file), recursive = TRUE, showWarnings = FALSE)
tmp <- paste0(cfg$output_file, ".tmp")
fwrite(result, tmp, sep = "\t", na = "NA")
if (!file.rename(tmp, cfg$output_file)) {
  unlink(tmp)
  fail("Could not atomically publish diagnostic output")
}
print(result)
cat(sprintf("BRAIN6_LAVA_SINGLE_TRAIT_DIAGNOSTIC status=%s trait=%s locus=%s\n",
  result$status[[1L]], cfg$trait_id, cfg$locus_id))
if (result$status[[1L]] == "FAILED") {
  fail(paste("Single-trait LAVA diagnostic failed:", result$reason[[1L]]))
}
