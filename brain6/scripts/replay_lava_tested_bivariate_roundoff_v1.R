#!/usr/bin/env Rscript
# Post-run detail replay: only loci/pairs already marked TESTED in the frozen run.
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(LAVA))
fail <- function(message) stop(message, call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) fail("Usage: Rscript replay_lava_tested_bivariate.R CONFIG.json")
cfg <- jsonlite::fromJSON(args[[1]], simplifyVector = FALSE)
need <- function(test, message) if (!isTRUE(test)) fail(message)
need(as.character(packageVersion("LAVA")) == "0.1.5", "LAVA version drift")
need(as.character(getRversion()) == "4.3.3", "R version drift")
need(cfg$analysis_id == "brain6-lava-local-rg-v2", "Unexpected LAVA analysis identity")
if (isTRUE(cfg$roundoff_patch)) {
# Numerical-stability patch for the pinned LAVA 0.1.5 blockwise LD path.
# The projected matrix is mathematically symmetric; symmetrize only floating-point
# roundoff before the recursive eigendecomposition, and fail on material asymmetry.
install_lava_roundoff_patch <- function() {
  ns <- asNamespace("LAVA")
  original <- get("decompose.ld", envir = ns, inherits = FALSE)
  body_lines <- deparse(body(original), width.cutoff = 500L)
  target <- "M = t(R.base) %*% ld$ld %*% R.base"
  hit <- which(trimws(body_lines) == target)
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
}
loci <- LAVA::read.loci(cfg$loci_file)
result_rows <- vector("list", length(cfg$pairs))
for (i in seq_along(cfg$pairs)) {
  pair <- cfg$pairs[[i]]
  set.seed(as.integer(cfg$random_seed))
  phenotypes <- as.character(unlist(pair$phenotypes, use.names = FALSE))
  need(length(phenotypes) == 2L, paste("Invalid phenotype pair", pair$pair_id))
  locus_index <- match(as.character(pair$locus_id), as.character(loci$LOC))
  need(!is.na(locus_index) && sum(as.character(loci$LOC) == pair$locus_id) == 1L,
       paste("Unknown or duplicate locus", pair$locus_id))
  input <- LAVA::process.input(input.info.file = pair$input_info,
    sample.overlap.file = pair$sample_overlap_file, ref.prefix = pair$reference_prefix,
    phenos = phenotypes, input.dir = dirname(pair$input_info))
  loc <- LAVA::process.locus(loci[locus_index, , drop = FALSE], input,
    phenos = phenotypes, min.K = as.integer(cfg$locus_processing$min_K),
    prune.thresh = as.numeric(cfg$locus_processing$prune_thresh),
    max.prop.K = as.numeric(cfg$locus_processing$max_prop_K),
    drop.failed = isTRUE(cfg$locus_processing$drop_failed),
    max.block.size = as.integer(cfg$locus_processing$max_block_size),
    cap.estimates = isTRUE(cfg$locus_processing$cap_estimates))
  need(!is.null(loc), paste("LAVA process.locus returned NULL at", pair$pair_id, pair$locus_id))
  univ <- as.data.table(LAVA::run.univ(loc, phenos = phenotypes,
    cap.estimates = isTRUE(cfg$univariate$cap_estimates)))
  need(nrow(univ) == 2L && all(univ$p < cfg$univariate$gate_p_strictly_less_than),
       paste("Frozen univariate gate no longer passes at", pair$pair_id, pair$locus_id))
  biv <- as.data.table(LAVA::run.bivar(loc, phenos = phenotypes,
    adap.thresh = as.numeric(unlist(cfg$bivariate$adap.thresh)),
    p.values = isTRUE(cfg$bivariate$p_values), CIs = isTRUE(cfg$bivariate$confidence_intervals),
    param.lim = as.numeric(cfg$bivariate$parameter_limit),
    cap.estimates = isTRUE(cfg$bivariate$cap_estimates)))
  need(nrow(biv) == 1L && all(c("rho", "rho.lower", "rho.upper", "p") %in% names(biv)),
       paste("Unexpected LAVA bivariate detail schema at", pair$pair_id, pair$locus_id))
  result_rows[[i]] <- data.table(pair_id = pair$pair_id, locus_id = pair$locus_id,
    rho = biv$rho[[1L]], rho.lower = biv$rho.lower[[1L]], rho.upper = biv$rho.upper[[1L]],
    r2 = biv$r2[[1L]], r2.lower = biv$r2.lower[[1L]], r2.upper = biv$r2.upper[[1L]],
    p = biv$p[[1L]])
}
result <- rbindlist(result_rows, use.names = TRUE)
need(!anyDuplicated(result[, .(pair_id, locus_id)]), "Duplicate replayed LAVA result")
fwrite(result, cfg$output_path, sep = "\t", na = "NA")
cat(sprintf("BRAIN6_LAVA_DETAILS_REPLAYED n=%d\n", nrow(result)))
