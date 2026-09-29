#!/usr/bin/env Rscript
# Persistent chromosome-scoped worker for canonical, pair-independent LAVA h2 cells.
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(LAVA))

fail <- function(message) stop(message, call. = FALSE)
need <- function(test, message) if (!isTRUE(test)) fail(message)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) fail("Usage: Rscript run_lava_canonical_batch_v3.R CONFIG.json")
cfg <- fromJSON(args[[1L]], simplifyVector = FALSE)
need(as.character(packageVersion("LAVA")) == "0.1.5", "LAVA version drift")
need(isTRUE(cfg$runtime_validator_passed) && grepl("^[0-9a-f]{64}$", cfg$runtime_validation_output_sha256),
     "Pinned LAVA runtime validation receipt is missing")
need(cfg$analysis_id == "brain6-lava-canonical-v3", "Unexpected analysis identity")
need(length(cfg$locus_ids) > 0L && length(unique(unlist(cfg$locus_ids))) == length(unlist(cfg$locus_ids)),
     "Batch locus ids must be nonempty and unique")
need(length(cfg$trait_ids) == 7L && length(unique(unlist(cfg$trait_ids))) == 7L,
     "Canonical family must contain exactly seven unique traits")
need(dir.exists(cfg$input_root) && file.exists(cfg$loci_file) && dir.exists(cfg$output_dir),
     "Worker input/output path is missing")

# Same roundoff-only fix as the separate v2 follow-up; exact source anchor protects against drift.
ns <- asNamespace("LAVA")
original <- get("decompose.ld", envir = ns, inherits = FALSE)
body_lines <- deparse(body(original), width.cutoff = 500L)
hit <- which(trimws(body_lines) == "M = t(R.base) %*% ld$ld %*% R.base")
need(length(hit) == 1L, "Pinned LAVA blockwise decomposition source changed")
insert <- c(
  "            asymmetry <- max(abs(M - t(M)))",
  "            tolerance <- 100 * .Machine$double.eps * max(1, nrow(M)) * max(1, max(abs(M)))",
  "            if (!is.finite(asymmetry) || asymmetry > tolerance) stop('Material asymmetry in block-reduced LD matrix; refusing numerical symmetrization', call. = FALSE)",
  "            if (asymmetry > 0) M <- (M + t(M)) / 2"
)
body_lines <- append(body_lines, insert, after = hit)
patched <- original
body(patched) <- parse(text = paste(body_lines, collapse = "\n"))[[1L]]
environment(patched) <- ns
unlockBinding("decompose.ld", ns); assign("decompose.ld", patched, envir = ns); lockBinding("decompose.ld", ns)
set.seed(as.integer(cfg$random_seed))

loci <- LAVA::read.loci(cfg$loci_file)
trait_ids <- as.character(unlist(cfg$trait_ids, use.names = FALSE))
input_root <- cfg$input_root
output_dir <- cfg$output_dir
policy <- cfg$execution_policy$locus_processing
univ_policy <- cfg$execution_policy$univariate
lava_ns <- asNamespace("LAVA")
check_reference <- get("check.reference", envir = lava_ns, inherits = FALSE)
load_reference <- get("load.reference", envir = lava_ns, inherits = FALSE)
reference_timer <- system.time(reference <- load_reference(check_reference(cfg$reference_prefix)))
read_sumstats <- get("read.sumstats.file", envir = lava_ns, inherits = FALSE)
harmonize_snps <- get("harmonize.snps", envir = lava_ns, inherits = FALSE)
align_sumstats <- get("align", envir = lava_ns, inherits = FALSE)
cat(sprintf("BRAIN6_WORKER_READY reference_load_seconds=%.6f\n", reference_timer[["elapsed"]]))
compute_seconds <- 0
result_write_seconds <- 0

atomic_fwrite <- function(x, path) {
  tmp <- paste0(path, ".", Sys.getpid(), ".partial")
  fwrite(x, tmp, sep = "\t", na = "NA")
  if (!file.rename(tmp, path)) { unlink(tmp); fail(paste("Cannot publish cell output:", path)) }
}

for (locus_id in as.character(unlist(cfg$locus_ids, use.names = FALSE))) {
  locus_timer <- proc.time()[["elapsed"]]
  locus_index <- match(locus_id, as.character(loci$LOC))
  need(!is.na(locus_index) && sum(as.character(loci$LOC) == locus_id) == 1L,
       paste("Missing or duplicated locus:", locus_id))
  locus_row <- loci[locus_index, , drop = FALSE]
  locus_dir <- file.path(input_root, paste0("locus_", locus_id))
  result_path <- file.path(output_dir, "cells", paste0("locus_", locus_id, ".tsv"))
  dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
  info <- fread(file.path(locus_dir, "input_info.tsv"))
  need(all(c("phenotype", "cases", "controls", "filename") %in% names(info)) &&
         setequal(info$phenotype, trait_ids), paste("Input-info trait coverage differs:", locus_id))
  info$filename <- file.path(locus_dir, info$filename)
  info$N <- info$cases + info$controls
  info$prop_cases <- info$cases / info$N
  info$binary <- !is.na(info$prop_cases) & info$prop_cases != 1
  rows <- lapply(trait_ids, function(trait) {
    tryCatch({
      trait_info <- info[match(trait, info$phenotype), , drop = FALSE]
      need(file.exists(trait_info$filename[[1L]]), paste("Missing sumstats:", trait, locus_id))
      input <- new.env(parent = globalenv())
      input$info <- data.table::copy(trait_info)
      input$P <- 1L
      input$reference <- reference
      input$sample.overlap <- NULL
      input$sum.stats <- setNames(list(read_sumstats(trait_info$filename[[1L]], trait)), trait)
      harmonize_snps(input); align_sumstats(input)
      messages <- capture.output(loc <- LAVA::process.locus(
        locus_row, input, phenos = trait,
        min.K = as.integer(policy$min_K), prune.thresh = as.numeric(policy$prune_thresh),
        max.prop.K = as.numeric(policy$max_prop_K), drop.failed = isTRUE(policy$drop_failed),
        max.block.size = as.integer(policy$max_block_size), cap.estimates = isTRUE(policy$cap_estimates)))
      if (is.null(loc)) {
        fewer_than_min_k <- any(grepl("Fewer than", messages, fixed = TRUE))
        negative_local_h2 <- any(grepl("Negative variance estimate", messages, fixed = TRUE))
        why <- if (fewer_than_min_k) "FEWER_THAN_MIN_K" else
          if (negative_local_h2) "LOW_LOCAL_H2_UNDERPOWERED" else "LAVA_PROCESS_LOCUS_RETURNED_NULL"
        return(data.table(phen = trait, locus_id = locus_id, status = "NOT_RUN", n_snps = NA_integer_,
          n_components = NA_integer_, h2.obs = NA_real_, h2.latent = NA_real_, p = NA_real_, reason = why))
      }
      u <- as.data.table(LAVA::run.univ(loc, phenos = trait, cap.estimates = isTRUE(univ_policy$cap_estimates)))
      need(nrow(u) == 1L && identical(as.character(u$phen[[1L]]), trait) &&
             is.finite(u$p[[1L]]) && u$p[[1L]] >= 0 && u$p[[1L]] <= 1,
           paste("Invalid canonical univariate output:", trait, locus_id))
      data.table(phen = trait, locus_id = locus_id, status = "TESTED", n_snps = as.integer(loc$n.snps),
        n_components = as.integer(loc$K), h2.obs = as.numeric(u$h2.obs[[1L]]),
        h2.latent = if ("h2.latent" %in% names(u)) as.numeric(u$h2.latent[[1L]]) else NA_real_,
        p = as.numeric(u$p[[1L]]), reason = "")
    }, error = function(e) {
      message <- conditionMessage(e)
      too_few <- grepl("Less than 3 SNPs shared across data sets", message, fixed = TRUE)
      data.table(phen = trait, locus_id = locus_id,
        status = if (too_few) "NOT_RUN" else "FAILED",
        n_snps = NA_integer_, n_components = NA_integer_, h2.obs = NA_real_, h2.latent = NA_real_, p = NA_real_,
        reason = if (too_few) "FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS" else substr(message, 1L, 500L))
    })
  })
  cell_rows <- rbindlist(rows, use.names = TRUE)
  need(nrow(cell_rows) == 7L && setequal(cell_rows$phen, trait_ids), paste("Incomplete locus cells:", locus_id))
  compute_seconds <- compute_seconds + proc.time()[["elapsed"]] - locus_timer
  write_timer <- system.time(atomic_fwrite(cell_rows, result_path))
  result_write_seconds <- result_write_seconds + write_timer[["elapsed"]]
  cat(sprintf("BRAIN6_CANONICAL_LOCUS_METRICS %s compute_wall_seconds=%.6f write_seconds=%.6f\n",
    locus_id, proc.time()[["elapsed"]] - locus_timer, write_timer[["elapsed"]]))
  cat(sprintf("BRAIN6_CANONICAL_LOCUS_COMPLETE %s\n", locus_id))
}
cat(sprintf("BRAIN6_CANONICAL_BATCH_COMPLETE reference_load_seconds=%.6f compute_seconds=%.6f result_write_seconds=%.6f\n",
  reference_timer[["elapsed"]], compute_seconds, result_write_seconds))
