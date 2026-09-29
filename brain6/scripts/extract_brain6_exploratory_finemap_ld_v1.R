#!/usr/bin/env Rscript
# Frozen 100k-UKB LD extraction for the predeclared exploratory RSS inputs.
suppressPackageStartupMessages(library(LAVA))
if (as.character(packageVersion("LAVA")) != "0.1.5") stop("Expected LAVA 0.1.5")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: script <repository-output-dir> <reference-prefix>")
out <- normalizePath(args[[1]], mustWork = TRUE)
prefix <- args[[2]]
trait <- utils::read.delim(file.path(out, "fine_mapping_input_manifest.tsv"), check.names = FALSE)
pair <- utils::read.delim(file.path(out, "trait_coloc_input_manifest.tsv"), check.names = FALSE)
tasks <- data.frame(
  analysis_id = c(as.character(trait$analysis_id[trait$status == "READY"]),
                  paste0(as.character(pair$candidate_locus_id[pair$status == "READY"]), "__shared")),
  chromosome = c(as.integer(trait$chromosome[trait$status == "READY"]),
                 as.integer(sub(".*_chr([0-9]+)_[0-9]+_[0-9]+$", "\\1",
                                as.character(pair$candidate_locus_id[pair$status == "READY"])))),
  input_tsv = c(as.character(trait$input_tsv[trait$status == "READY"]),
                as.character(pair$trait1_input_tsv[pair$status == "READY"])),
  ld_rds = c(as.character(trait$ld_rds[trait$status == "READY"]),
             as.character(pair$ld_rds[pair$status == "READY"])),
  stringsAsFactors = FALSE
)
if (anyNA(tasks$chromosome) || anyDuplicated(tasks$analysis_id)) stop("Invalid LD task manifest")
if (!nrow(tasks)) stop("No admitted LD tasks")
results <- list()
for (chr in sort(unique(tasks$chromosome))) {
  reference <- LAVA:::load.reference(list(prefix = prefix, mode = "ld", chromosomes = chr))
  if (reference$sample.size != 100000) stop("Unexpected reference sample size")
  for (i in which(tasks$chromosome == chr)) {
    row <- tasks[i, ]
    result <- tryCatch({
      input <- utils::read.delim(row$input_tsv, check.names = FALSE)
      snps <- tolower(as.character(input$SNP))
      if (length(snps) < 500L || anyDuplicated(snps)) stop("Invalid admitted SNP set")
      ld <- LAVA:::read.ld(reference, snps)
      matrix <- as.matrix(ld$ld)
      if (!setequal(rownames(matrix), snps) || !setequal(colnames(matrix), snps))
        stop("LD SNP identity mismatch")
      matrix <- matrix[snps, snps, drop = FALSE]
      if (!all(is.finite(matrix)) || max(abs(matrix - t(matrix))) > 1e-6 ||
          max(abs(diag(matrix) - 1)) > 1e-6)
        stop("Malformed reference LD matrix")
      dir.create(dirname(row$ld_rds), recursive = TRUE, showWarnings = FALSE)
      saveRDS(list(snp = snps, ld = matrix, reference_n = reference$sample.size,
                   LAVA_version = "0.1.5", chromosome = chr), row$ld_rds,
              compress = FALSE)
      data.frame(analysis_id = row$analysis_id, chromosome = chr, status = "READY",
                 reason = "", n_snps = length(snps), raw_max_abs_r = max(abs(matrix)),
                 ld_rds = row$ld_rds,
                 stringsAsFactors = FALSE)
    }, error = function(e) {
      data.frame(analysis_id = row$analysis_id, chromosome = chr,
                 status = "RUNTIME_FAILURE", reason = conditionMessage(e),
                 n_snps = NA_integer_, raw_max_abs_r = NA_real_, ld_rds = row$ld_rds,
                 stringsAsFactors = FALSE)
    })
    results[[length(results) + 1L]] <- result
    cat(row$analysis_id, result$status, result$n_snps, "\n")
    flush.console()
  }
  rm(reference)
  gc()
}
result <- do.call(rbind, results)
utils::write.table(result, file.path(out, "ld_extraction_manifest.tsv"),
                   sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   na = "")
if (any(result$status != "READY")) stop("One or more LD extraction tasks failed; see manifest")
