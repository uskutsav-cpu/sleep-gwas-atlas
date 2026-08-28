#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(LAVA))

if (as.character(packageVersion("LAVA")) != "0.1.5") {
  stop("LAVA 0.1.5 is required")
}

defaults <- list(
  input_info = "results/tables/lava_input_info.tsv",
  sample_overlap = "results/tables/lava_sample_overlap.txt",
  pair_manifest = "results/tables/lava_pair_manifest.tsv",
  input_provenance = "results/tables/lava_input_provenance.tsv",
  runtime_policy = "results/tables/lava_runtime_policy.tsv",
  checkpoint_dir = "results/checkpoints/lava",
  out_dir = "results/tables",
  locus_start = 1L,
  locus_end = NA_integer_,
  finalize_only = FALSE
)

parse_args <- function(args, values) {
  i <- 1L
  while (i <= length(args)) {
    argument <- args[[i]]
    if (argument == "--finalize-only") {
      values$finalize_only <- TRUE
      i <- i + 1L
      next
    }
    mapping <- c(
      "--input-info" = "input_info", "--sample-overlap" = "sample_overlap",
      "--pair-manifest" = "pair_manifest", "--input-provenance" = "input_provenance",
      "--runtime-policy" = "runtime_policy", "--checkpoint-dir" = "checkpoint_dir",
      "--out-dir" = "out_dir", "--locus-start" = "locus_start", "--locus-end" = "locus_end"
    )
    if (!(argument %in% names(mapping)) || i == length(args)) {
      stop(paste("Unknown or incomplete argument:", argument))
    }
    key <- unname(mapping[[argument]])
    values[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  values$locus_start <- as.integer(values$locus_start)
  if (!is.na(values$locus_end)) values$locus_end <- as.integer(values$locus_end)
  values
}

values <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
required_files <- c(
  values$input_info, values$sample_overlap, values$pair_manifest,
  values$input_provenance, values$runtime_policy
)
if (!all(file.exists(required_files))) {
  stop(paste("Missing prepared LAVA input(s):", paste(required_files[!file.exists(required_files)], collapse = ", ")))
}

policy_table <- read.delim(values$runtime_policy, stringsAsFactors = FALSE, check.names = FALSE)
policy <- setNames(as.list(policy_table$value), policy_table$key)
numeric_policy <- function(key) as.numeric(policy[[key]])
integer_policy <- function(key) as.integer(policy[[key]])
logical_policy <- function(key) identical(tolower(policy[[key]]), "true")
if (policy$lava_version != "0.1.5") stop("Prepared policy does not pin LAVA 0.1.5")

input_info <- read.delim(values$input_info, stringsAsFactors = FALSE, check.names = FALSE)
pairs <- read.delim(values$pair_manifest, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(input_info) != integer_policy("expected_traits") || length(unique(input_info$phenotype)) != 45L) {
  stop("Input info is not the exact 45-trait panel")
}
if (nrow(pairs) != integer_policy("expected_sleep_non_sleep_pairs")) {
  stop("Pair manifest is not the exact 396-pair family")
}

contract_output <- suppressWarnings(system2(
  "python3", c("scripts/lava_contract.py", "--run-fingerprint"),
  stdout = TRUE, stderr = TRUE
))
contract_status <- attr(contract_output, "status")
if (!is.null(contract_status) && contract_status != 0L) {
  stop(paste("LAVA immutable execution contract failed:", paste(contract_output, collapse = " | ")))
}
fingerprint_lines <- contract_output[grepl("^[0-9a-f]{64}$", contract_output)]
if (length(fingerprint_lines) != 1L) stop("LAVA contract did not return one run fingerprint")
fingerprint <- fingerprint_lines[[1L]]

reference_prefix <- policy$reference_prefix
reference_files <- unlist(lapply(1:22, function(chromosome) {
  paste0(reference_prefix, "_chr", chromosome, c(".info", ".bcor"))
}))
if (!all(file.exists(reference_files))) {
  missing <- reference_files[!file.exists(reference_files)]
  stop(paste0(
    "Complete LAVA UKB v1.1 reference is required; missing ", length(missing),
    " chromosome file(s), beginning with: ", paste(head(missing, 4), collapse = ", ")
  ))
}

locus_file <- policy$locus_definition
if (!file.exists(locus_file)) stop(paste("Missing locus definition:", locus_file))
loci <- read.loci(locus_file)
expected_loci <- integer_policy("expected_loci")
if (nrow(loci) != expected_loci || anyDuplicated(loci$LOC)) stop("Locus definition is not the exact 2,495-locus file")
if (is.na(values$locus_end)) values$locus_end <- expected_loci
if (values$locus_start < 1L || values$locus_end > expected_loci || values$locus_start > values$locus_end) {
  stop("Invalid locus range")
}

dir.create(values$checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(values$out_dir, recursive = TRUE, showWarnings = FALSE)

empty_univ <- function(locus_row, status, message = NA_character_, n_snps = NA_integer_, K = NA_integer_) {
  data.frame(
    LOC = locus_row$LOC, CHR = locus_row$CHR, START = locus_row$START, STOP = locus_row$STOP,
    phen = input_info$phenotype, h2.obs = NA_real_, h2.latent = NA_real_,
    ascertained = NA, p = NA_real_, analysis_status = status, error = message,
    n_snps = n_snps, K = K, stringsAsFactors = FALSE
  )
}

empty_bivar <- function() {
  data.frame(
    LOC = numeric(0), CHR = numeric(0), START = numeric(0), STOP = numeric(0),
    sleep_trait = character(0), non_sleep_trait = character(0),
    analysis_tier = character(0), interpretation_status = character(0),
    global_rg = numeric(0), global_rg_se = numeric(0), global_rg_p = numeric(0), global_rg_fdr = numeric(0),
    rho = numeric(0), rho.lower = numeric(0), rho.upper = numeric(0),
    r2 = numeric(0), r2.lower = numeric(0), r2.upper = numeric(0), p = numeric(0),
    analysis_status = character(0), error = character(0), stringsAsFactors = FALSE
  )
}

process_one <- function(index, input) {
  locus_row <- loci[index, , drop = FALSE]
  started <- Sys.time()
  locus <- tryCatch(
    process.locus(
      locus_row, input,
      min.K = integer_policy("min_K"),
      prune.thresh = numeric_policy("prune_threshold"),
      max.prop.K = numeric_policy("max_proportion_K"),
      drop.failed = TRUE,
      max.block.size = integer_policy("max_block_size"),
      cap.estimates = logical_policy("cap_estimates")
    ),
    error = function(error) structure(list(message = conditionMessage(error)), class = "lava_process_error")
  )
  if (inherits(locus, "lava_process_error") || is.null(locus)) {
    message <- if (inherits(locus, "lava_process_error")) locus$message else "process.locus returned NULL"
    univ <- empty_univ(locus_row, "LOCUS_PROCESS_FAILED", message)
    status <- "PROCESS_FAILED"
    bivar <- empty_bivar()
    n_snps <- K <- NA_integer_
  } else {
    n_snps <- as.integer(locus$n.snps)
    K <- as.integer(locus$K)
    univ_result <- tryCatch(run.univ(locus), error = function(error) structure(list(message = conditionMessage(error)), class = "lava_univ_error"))
    if (inherits(univ_result, "lava_univ_error")) {
      univ <- empty_univ(locus_row, "UNIVARIATE_FAILED", univ_result$message, n_snps, K)
      status <- "UNIVARIATE_FAILED"
      bivar <- empty_bivar()
    } else {
      univ <- empty_univ(locus_row, "PHENOTYPE_DROPPED", "phenotype removed during locus processing", n_snps, K)
      matched <- match(univ_result$phen, univ$phen)
      univ$h2.obs[matched] <- univ_result$h2.obs
      univ$h2.latent[matched] <- univ_result$h2.latent
      univ$ascertained[matched] <- univ_result$ascertained
      univ$p[matched] <- univ_result$p
      univ$analysis_status[matched] <- "TESTED"
      univ$error[matched] <- NA_character_
      eligible_traits <- univ$phen[univ$analysis_status == "TESTED" & is.finite(univ$p) & univ$p <= numeric_policy("univariate_p_threshold")]
      eligible_pairs <- pairs[pairs$sleep_trait %in% eligible_traits & pairs$non_sleep_trait %in% eligible_traits, , drop = FALSE]
      bivar_rows <- vector("list", nrow(eligible_pairs))
      if (nrow(eligible_pairs) > 0L) {
        for (pair_index in seq_len(nrow(eligible_pairs))) {
          pair <- eligible_pairs[pair_index, , drop = FALSE]
          set.seed(integer_policy("random_seed") + index * 1000L + pair_index)
          result <- tryCatch(
            run.bivar(locus, phenos = c(pair$sleep_trait, pair$non_sleep_trait)),
            error = function(error) structure(list(message = conditionMessage(error)), class = "lava_bivar_error")
          )
          base <- data.frame(
            LOC = locus_row$LOC, CHR = locus_row$CHR, START = locus_row$START, STOP = locus_row$STOP,
            sleep_trait = pair$sleep_trait, non_sleep_trait = pair$non_sleep_trait,
            analysis_tier = pair$analysis_tier, interpretation_status = pair$interpretation_status,
            global_rg = as.numeric(pair$global_rg), global_rg_se = as.numeric(pair$global_rg_se),
            global_rg_p = as.numeric(pair$global_rg_p), global_rg_fdr = as.numeric(pair$global_rg_fdr),
            stringsAsFactors = FALSE
          )
          if (inherits(result, "lava_bivar_error") || is.null(result) || nrow(result) != 1L) {
            message <- if (inherits(result, "lava_bivar_error")) result$message else "run.bivar returned no unique result"
            bivar_rows[[pair_index]] <- cbind(
              base, rho = NA_real_, rho.lower = NA_real_, rho.upper = NA_real_,
              r2 = NA_real_, r2.lower = NA_real_, r2.upper = NA_real_, p = NA_real_,
              analysis_status = "BIVARIATE_FAILED", error = message
            )
          } else {
            bivar_rows[[pair_index]] <- cbind(
              base, rho = result$rho, rho.lower = result$rho.lower, rho.upper = result$rho.upper,
              r2 = result$r2, r2.lower = result$r2.lower, r2.upper = result$r2.upper, p = result$p,
              analysis_status = "TESTED", error = NA_character_
            )
          }
        }
      }
      bivar <- if (length(bivar_rows)) do.call(rbind, bivar_rows) else empty_bivar()
      status <- "PROCESSED"
    }
  }
  locus_status <- data.frame(
    LOC = locus_row$LOC, CHR = locus_row$CHR, START = locus_row$START, STOP = locus_row$STOP,
    status = status, n_snps = n_snps, K = K,
    univariate_tested = sum(univ$analysis_status == "TESTED"),
    eligible_bivariate_pairs = nrow(bivar),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    stringsAsFactors = FALSE
  )
  list(fingerprint = fingerprint, status = locus_status, univ = univ, bivar = bivar)
}

if (!values$finalize_only) {
  input <- process.input(
    input.info.file = values$input_info,
    sample.overlap.file = values$sample_overlap,
    ref.prefix = reference_prefix,
    phenos = input_info$phenotype
  )
  for (index in seq.int(values$locus_start, values$locus_end)) {
    checkpoint <- file.path(values$checkpoint_dir, sprintf("locus_%04d.rds", index))
    if (file.exists(checkpoint)) {
      existing <- readRDS(checkpoint)
      if (!identical(existing$fingerprint, fingerprint)) stop(paste("Stale checkpoint:", checkpoint))
      message(sprintf("[%d/%d] locus %s already complete", index, expected_loci, loci$LOC[index]))
      next
    }
    message(sprintf("[%d/%d] processing locus %s", index, expected_loci, loci$LOC[index]))
    result <- process_one(index, input)
    temporary <- paste0(checkpoint, ".tmp")
    saveRDS(result, temporary, version = 3)
    if (!file.rename(temporary, checkpoint)) stop(paste("Could not publish checkpoint:", checkpoint))
  }
}

checkpoint_files <- file.path(values$checkpoint_dir, sprintf("locus_%04d.rds", seq_len(expected_loci)))
if (!all(file.exists(checkpoint_files))) {
  message(sprintf("LAVA range complete; %d/%d locus checkpoints exist. Canonical outputs are not published.", sum(file.exists(checkpoint_files)), expected_loci))
  quit(save = "no", status = 3L)
}
results <- lapply(checkpoint_files, readRDS)
if (any(!vapply(results, function(result) identical(result$fingerprint, fingerprint), logical(1)))) {
  stop("One or more LAVA checkpoints have a stale input fingerprint")
}
locus_status <- do.call(rbind, lapply(results, `[[`, "status"))
univ <- do.call(rbind, lapply(results, `[[`, "univ"))
bivar <- do.call(rbind, lapply(results, `[[`, "bivar"))
if (nrow(locus_status) != expected_loci || nrow(univ) != integer_policy("planned_univariate_tests")) {
  stop("Checkpoint collation did not yield the complete planned test family")
}
univ$univariate_test_family_n <- integer_policy("planned_univariate_tests")
univ$univariate_p_threshold <- numeric_policy("univariate_p_threshold")
univ$p_bonferroni <- pmin(1, univ$p * integer_policy("planned_univariate_tests"))
tested_univ <- which(univ$analysis_status == "TESTED" & is.finite(univ$p))
univ$p_fdr <- NA_real_
univ$p_fdr[tested_univ] <- p.adjust(univ$p[tested_univ], method = "BH")
bivar$bivariate_test_family_n <- rep.int(nrow(bivar), nrow(bivar))
bivar$p_fdr <- rep(NA_real_, nrow(bivar))
tested_bivar <- which(bivar$analysis_status == "TESTED" & is.finite(bivar$p))
bivar$p_fdr[tested_bivar] <- p.adjust(bivar$p[tested_bivar], method = "BH")
bivar$fdr_significant <- !is.na(bivar$p_fdr) & bivar$p_fdr <= numeric_policy("bivariate_fdr_alpha")
locus_status$analysis_fingerprint <- fingerprint
locus_status$lava_version <- as.character(packageVersion("LAVA"))
locus_status$reference_prefix <- reference_prefix

write_atomic <- function(table, path) {
  temporary <- paste0(path, ".tmp")
  write.table(table, temporary, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  if (!file.rename(temporary, path)) stop(paste("Could not publish:", path))
}
write_atomic(locus_status, file.path(values$out_dir, "lava_locus_status.tsv"))
write_atomic(univ, file.path(values$out_dir, "lava_univariate.tsv"))
write_atomic(bivar, file.path(values$out_dir, "lava_bivariate.tsv"))
message(sprintf("Published complete LAVA analysis: %d loci, %d univariate rows, %d eligible bivariate rows", nrow(locus_status), nrow(univ), nrow(bivar)))
