#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(LAVA))
suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop("Usage: 44_run_lava_sensitivity_locus.R <trait> <pair_order> <locus_index> <input_info> <overlap_file> <locus_file> <reference_prefix> <input_directory> <receipt_json>")
}
trait <- args[[1L]]
pair_order <- as.integer(args[[2L]])
locus_index <- as.integer(args[[3L]])
input_info_file <- args[[4L]]
overlap_file <- args[[5L]]
locus_file <- args[[6L]]
reference_prefix <- args[[7L]]
input_directory <- args[[8L]]
receipt_file <- args[[9L]]

if (!identical(as.character(getRversion()), "4.3.3")) stop("R 4.3.3 is required")
if (!identical(as.character(packageVersion("LAVA")), "0.1.5")) stop("LAVA 0.1.5 is required")
if (!grepl("^(accel_sleep_duration|chronotype|insomnia|longsleep|napping|shortsleep|sleep_apnea|sleep_efficiency|sleep_timing|sleepdur|sleepiness|snoring)$", trait)) stop("Trait is outside the frozen sleep panel")
if (is.na(pair_order) || pair_order < 1L || pair_order > 12L) stop("Invalid pair order")
if (is.na(locus_index) || locus_index < 1L || locus_index > 2495L) stop("Invalid locus index")
if (file.exists(receipt_file)) stop(paste("Refusing to overwrite immutable receipt:", receipt_file))
for (path in c(input_info_file, overlap_file, locus_file, input_directory)) if (!file.exists(path)) stop(paste("Missing input:", path))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1L]
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), "../.."), mustWork = TRUE)
runtime_validator <- file.path(repo_root, "scripts/133_validate_track_b_lava_runtime.R")
validator_output <- system2(file.path(R.home("bin"), "Rscript"), runtime_validator, stdout = TRUE, stderr = TRUE)
validator_status <- attr(validator_output, "status")
if (!is.null(validator_status) && validator_status != 0L) stop(paste("Pinned LAVA runtime validation failed:", paste(validator_output, collapse = " | ")))

input_info <- read.delim(input_info_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(input_info$phenotype, c("frailty", trait))) stop("Pair input.info must contain FI then the selected sleep trait")
loci <- read.loci(locus_file)
if (nrow(loci) != 2495L || anyDuplicated(loci$LOC)) stop("Frozen locus definition drifted")
locus_row <- loci[locus_index, , drop = FALSE]
started <- Sys.time()
seed <- 20260924L + pair_order * 100000L + locus_index

result <- list(
  schema_version = 1L, trait = trait, pair_order = pair_order, locus_index = locus_index,
  LOC = as.character(locus_row$LOC), CHR = as.integer(locus_row$CHR),
  START = as.integer(locus_row$START), STOP = as.integer(locus_row$STOP), seed = seed,
  process_status = "NOT_RUN", univ_status = "NOT_RUN", fi_status = "NOT_RUN",
  sleep_status = "NOT_RUN", bivar_status = "NOT_ELIGIBLE",
  n_snps = NA_integer_, K = NA_integer_, fi_h2_obs = NA_real_, fi_p = NA_real_,
  sleep_h2_obs = NA_real_, sleep_p = NA_real_, rho = NA_real_, rho_lower = NA_real_,
  rho_upper = NA_real_, r2 = NA_real_, r2_lower = NA_real_, r2_upper = NA_real_,
  bivar_p = NA_real_, error = NA_character_, elapsed_seconds = NA_real_
)

fail <- function(field, value, message) {
  result[[field]] <<- value
  result$error <<- substr(as.character(message), 1L, 4000L)
}

chrom_prefix <- paste0(file.path(reference_prefix, basename(reference_prefix)), "_chr", as.integer(locus_row$CHR))
if (!all(file.exists(paste0(chrom_prefix, c(".info", ".bcor"))))) stop(paste("Missing LAVA reference chromosome:", chrom_prefix))
tryCatch({
  input <- process.input(input.info.file = input_info_file, sample.overlap.file = overlap_file,
                         ref.prefix = chrom_prefix, phenos = c("frailty", trait), input.dir = input_directory)
  locus <- tryCatch(process.locus(locus_row, input, min.K = 2L, prune.thresh = 99,
                                  max.prop.K = 0.75, drop.failed = TRUE,
                                  max.block.size = 3000L, cap.estimates = TRUE),
                    error = function(e) structure(list(message = conditionMessage(e)), class = "lava_process_error"))
  if (inherits(locus, "lava_process_error") || is.null(locus)) {
    message <- if (inherits(locus, "lava_process_error")) locus$message else "process.locus returned NULL"
    fail("process_status", "PROCESS_FAILED", message)
    result$fi_status <- "PROCESS_FAILED"
    result$sleep_status <- "PROCESS_FAILED"
  } else {
    result$n_snps <- as.integer(locus$n.snps)
    result$K <- as.integer(locus$K)
    result$process_status <- "PROCESSED"
    univ <- tryCatch(run.univ(locus), error = function(e) structure(list(message = conditionMessage(e)), class = "lava_univ_error"))
    if (inherits(univ, "lava_univ_error") || !is.data.frame(univ) || !all(c("phen", "h2.obs", "p") %in% names(univ))) {
      message <- if (inherits(univ, "lava_univ_error")) univ$message else "run.univ returned an invalid result"
      fail("univ_status", "UNIVARIATE_FAILED", message)
      result$fi_status <- "UNIVARIATE_FAILED"
      result$sleep_status <- "UNIVARIATE_FAILED"
    } else {
      fi_i <- match("frailty", univ$phen)
      sleep_i <- match(trait, univ$phen)
      valid_result <- function(i) {
        if (is.na(i)) return(FALSE)
        h <- as.numeric(univ$h2.obs[i])
        p <- as.numeric(univ$p[i])
        is.finite(h) && h >= 0 && is.finite(p) && p >= 0 && p <= 1
      }
      fi_valid <- valid_result(fi_i)
      sleep_valid <- valid_result(sleep_i)
      if (fi_valid) {
        result$fi_h2_obs <- as.numeric(univ$h2.obs[fi_i])
        result$fi_p <- as.numeric(univ$p[fi_i])
        result$fi_status <- "TESTED"
      } else {
        result$fi_status <- if (is.na(fi_i)) "PHENOTYPE_DROPPED" else "UNIVARIATE_INVALID"
      }
      if (sleep_valid) {
        result$sleep_h2_obs <- as.numeric(univ$h2.obs[sleep_i])
        result$sleep_p <- as.numeric(univ$p[sleep_i])
        result$sleep_status <- "TESTED"
      } else {
        result$sleep_status <- if (is.na(sleep_i)) "PHENOTYPE_DROPPED" else "UNIVARIATE_INVALID"
      }
      result$univ_status <- if (fi_valid && sleep_valid) "TESTED" else if (fi_valid || sleep_valid) "PARTIAL" else "PHENOTYPE_DROPPED"
      if (result$univ_status == "PHENOTYPE_DROPPED") {
        result$error <- "Neither phenotype yielded a valid local univariate result"
      } else if (result$univ_status == "PARTIAL") {
        dropped <- c(if (!fi_valid) "frailty" else NULL, if (!sleep_valid) trait else NULL)
        result$error <- paste("No valid local univariate result for:", paste(dropped, collapse = ", "))
      }
      threshold <- 4.45335114673792e-7
      if (fi_valid && sleep_valid && result$fi_p <= threshold && result$sleep_p <= threshold) {
        set.seed(seed)
        bivar <- tryCatch(run.bivar(locus, phenos = c("frailty", trait)), error = function(e) structure(list(message = conditionMessage(e)), class = "lava_bivar_error"))
        fields <- c("rho", "rho.lower", "rho.upper", "r2", "r2.lower", "r2.upper", "p")
        if (inherits(bivar, "lava_bivar_error") || !is.data.frame(bivar) || nrow(bivar) != 1L || !all(fields %in% names(bivar))) {
          message <- if (inherits(bivar, "lava_bivar_error")) bivar$message else "run.bivar returned an incomplete result"
          fail("bivar_status", "BIVARIATE_FAILED", message)
        } else {
          vals <- vapply(fields, function(field) as.numeric(bivar[[field]][1L]), numeric(1))
          names(vals) <- fields
          if (any(!is.finite(vals)) || vals[["p"]] < 0 || vals[["p"]] > 1) {
            fail("bivar_status", "BIVARIATE_INVALID", "run.bivar returned non-finite or out-of-range estimates")
          } else {
            for (field in fields[fields != "p"]) result[[gsub("\\.", "_", field)]] <- vals[[field]]
            result$bivar_p <- vals[["p"]]
            result$bivar_status <- "TESTED"
          }
        }
      }
    }
  }
}, error = function(e) fail("process_status", "PROCESS_FAILED", conditionMessage(e)))

result$elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
dir.create(dirname(receipt_file), recursive = TRUE, showWarnings = FALSE)
temporary <- paste0(receipt_file, ".", Sys.getpid(), ".partial")
write_json(result, temporary, auto_unbox = TRUE, na = "null", pretty = TRUE, digits = NA)
if (!file.rename(temporary, receipt_file)) stop(paste("Could not publish immutable receipt:", receipt_file))
cat(sprintf("LOCUS_RECEIPT trait=%s index=%d process=%s univ=%s bivar=%s seconds=%.2f\n",
            trait, locus_index, result$process_status, result$univ_status,
            result$bivar_status, result$elapsed_seconds))
