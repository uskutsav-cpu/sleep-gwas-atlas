#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

fail <- function(...) stop(paste0(...), call. = FALSE)

parse_args <- function(args) {
  defaults <- list(
    panel = "config/analysis_panel.tsv",
    munged_dir = "data/munged",
    ld_dir = "ref/eur_w_ld_chr",
    out_dir = "results/tables",
    log_prefix = "results/logs/genomicsem/ldsc_45_trait"
  )
  if (length(args) == 0L) return(defaults)
  if (length(args) %% 2 != 0) {
    fail("Arguments must be supplied as --name value pairs")
  }
  for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[[i]])
    key <- gsub("-", "_", key)
    if (!key %in% names(defaults)) fail("Unknown argument: ", args[[i]])
    defaults[[key]] <- args[[i + 1]]
  }
  defaults
}

write_matrix <- function(value, row_ids, col_ids, path, row_name) {
  if (!identical(dim(value), c(length(row_ids), length(col_ids)))) {
    fail("Matrix dimension mismatch for ", path)
  }
  frame <- data.frame(row_ids, value, check.names = FALSE)
  names(frame) <- c(row_name, col_ids)
  write.table(frame, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

write_gzip_matrix <- function(value, ids, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  frame <- data.frame(ids, value, check.names = FALSE)
  names(frame) <- c("element_id", ids)
  write.table(frame, connection, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
panel <- read.delim(args$panel, check.names = FALSE, na.strings = c("NA", ""))
required <- c("panel_version", "trait_id", "domain", "type", "ncase", "ncontrol", "n_total", "pop_prev")
missing_columns <- setdiff(required, names(panel))
if (length(missing_columns)) fail("Panel is missing columns: ", paste(missing_columns, collapse = ", "))
if (nrow(panel) != 45L || length(unique(panel$trait_id)) != 45L) {
  fail("Expected exactly 45 unique ordered traits; found ", nrow(panel), " rows")
}
if (sum(panel$domain == "sleep") != 12L) fail("Expected exactly 12 sleep traits")
if (length(unique(panel$panel_version)) != 1L || unique(panel$panel_version) != "atlas-v1.0") {
  fail("Panel version must be exactly atlas-v1.0")
}
if (any(!panel$type %in% c("binary", "continuous"))) fail("Unexpected trait type")

is_binary <- panel$type == "binary"
if (any(is_binary & (is.na(panel$ncase) | is.na(panel$ncontrol) | is.na(panel$pop_prev)))) {
  fail("Every binary trait needs ncase, ncontrol, and pop_prev")
}
if (any(is_binary & (panel$pop_prev <= 0 | panel$pop_prev >= 1))) {
  fail("Binary population prevalences must be strictly between zero and one")
}

trait_ids <- panel$trait_id
trait_files <- file.path(args$munged_dir, paste0(trait_ids, ".sumstats.gz"))
missing_files <- trait_files[!file.exists(trait_files) | file.info(trait_files)$size <= 0]
if (length(missing_files)) fail("Missing munged input(s): ", paste(missing_files, collapse = ", "))
for (chromosome in 1:22) {
  for (suffix in c(".l2.ldscore.gz", ".l2.M_5_50")) {
    reference <- file.path(args$ld_dir, paste0(chromosome, suffix))
    if (!file.exists(reference) || file.info(reference)$size <= 0) fail("Missing LD reference: ", reference)
  }
}

sample_prev <- rep(NA_real_, nrow(panel))
sample_prev[is_binary] <- panel$ncase[is_binary] / (panel$ncase[is_binary] + panel$ncontrol[is_binary])
population_prev <- rep(NA_real_, nrow(panel))
population_prev[is_binary] <- panel$pop_prev[is_binary]

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args$log_prefix), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(GenomicSEM))
if (getRversion() != "4.3.3") fail("Pinned R 4.3.3 required; found ", getRversion())
if (as.character(packageVersion("GenomicSEM")) != "0.0.5") {
  fail("Pinned GenomicSEM 0.0.5 required; found ", packageVersion("GenomicSEM"))
}

started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
covstruc <- GenomicSEM::ldsc(
  traits = trait_files,
  sample.prev = sample_prev,
  population.prev = population_prev,
  ld = paste0(sub("/+$", "", args$ld_dir), "/"),
  wld = paste0(sub("/+$", "", args$ld_dir), "/"),
  trait.names = trait_ids,
  ldsc.log = args$log_prefix,
  stand = TRUE
)
finished_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)

expected_elements <- as.integer(length(trait_ids) * (length(trait_ids) + 1L) / 2L)
if (!identical(dim(covstruc$S), c(45L, 45L))) fail("Genetic covariance matrix is not 45x45")
if (!identical(dim(covstruc$I), c(45L, 45L))) fail("Intercept matrix is not 45x45")
if (!identical(dim(covstruc$V), c(expected_elements, expected_elements))) {
  fail("Sampling covariance matrix is not ", expected_elements, "x", expected_elements)
}
if (any(!is.finite(covstruc$S)) || any(!is.finite(covstruc$V)) || any(!is.finite(covstruc$I))) {
  fail("Covariance structure contains non-finite values")
}

rownames(covstruc$S) <- colnames(covstruc$S) <- trait_ids
rownames(covstruc$I) <- colnames(covstruc$I) <- trait_ids
rownames(covstruc$S_Stand) <- colnames(covstruc$S_Stand) <- trait_ids

pair_i <- integer(expected_elements)
pair_j <- integer(expected_elements)
cursor <- 1L
for (j in seq_along(trait_ids)) {
  for (k in j:length(trait_ids)) {
    pair_i[[cursor]] <- j
    pair_j[[cursor]] <- k
    cursor <- cursor + 1L
  }
}
element_ids <- paste(trait_ids[pair_i], trait_ids[pair_j], sep = "__")
rownames(covstruc$V) <- colnames(covstruc$V) <- element_ids
names(covstruc$N) <- element_ids

se_elements <- sqrt(diag(covstruc$V))
pair_table <- data.frame(
  element_id = element_ids,
  trait_1 = trait_ids[pair_i],
  trait_2 = trait_ids[pair_j],
  genetic_covariance = covstruc$S[cbind(pair_i, pair_j)],
  se = se_elements,
  z = covstruc$S[cbind(pair_i, pair_j)] / se_elements,
  p = 2 * pnorm(abs(covstruc$S[cbind(pair_i, pair_j)] / se_elements), lower.tail = FALSE),
  genetic_correlation = covstruc$S_Stand[cbind(pair_i, pair_j)],
  cross_trait_intercept = covstruc$I[cbind(pair_i, pair_j)],
  effective_n = as.numeric(covstruc$N),
  scale_1 = ifelse(is_binary[pair_i], "liability", "observed"),
  scale_2 = ifelse(is_binary[pair_j], "liability", "observed")
)
pair_table$fdr_off_diagonal <- NA_real_
off_diagonal <- pair_table$trait_1 != pair_table$trait_2
pair_table$fdr_off_diagonal[off_diagonal] <- p.adjust(pair_table$p[off_diagonal], method = "BH")

trait_scales <- data.frame(
  order = seq_along(trait_ids),
  trait_id = trait_ids,
  type = panel$type,
  sample_prevalence = sample_prev,
  population_prevalence = population_prev,
  covariance_scale = ifelse(is_binary, "liability", "observed")
)

rds_tmp <- file.path(args$out_dir, "ldsc_covariance_structure.rds.tmp")
rds_path <- file.path(args$out_dir, "ldsc_covariance_structure.rds")
saveRDS(covstruc, rds_tmp, compress = "xz")
if (!file.rename(rds_tmp, rds_path)) fail("Could not atomically publish ", rds_path)

write_matrix(covstruc$S, trait_ids, trait_ids,
             file.path(args$out_dir, "ldsc_covariance_45x45.tsv"), "trait_id")
write_matrix(covstruc$S_Stand, trait_ids, trait_ids,
             file.path(args$out_dir, "ldsc_genetic_correlation_45x45.tsv"), "trait_id")
write_matrix(covstruc$I, trait_ids, trait_ids,
             file.path(args$out_dir, "ldsc_intercept_45x45.tsv"), "trait_id")
write_gzip_matrix(covstruc$V, element_ids,
                  file.path(args$out_dir, "ldsc_sampling_covariance_1035x1035.tsv.gz"))
write.table(pair_table, file.path(args$out_dir, "ldsc_covariance_pairs.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write.table(trait_scales, file.path(args$out_dir, "ldsc_covariance_trait_scales.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

metadata <- data.frame(
  field = c("panel_version", "trait_count", "covariance_elements", "jackknife_blocks",
            "ld_reference", "r_version", "genomicsem_version", "genomicsem_commit",
            "started_at_utc", "finished_at_utc", "m_total"),
  value = c("atlas-v1.0", "45", as.character(expected_elements), as.character(expected_elements + 1L),
            args$ld_dir, as.character(getRversion()), as.character(packageVersion("GenomicSEM")),
            "6b65ca5db39fdade08b0d811477be1cdd57b5039", started_at, finished_at,
            as.character(covstruc$m))
)
write.table(metadata, file.path(args$out_dir, "ldsc_covariance_metadata.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Published complete 45-trait GenomicSEM covariance structure\n")
cat("S: 45x45; V: ", expected_elements, "x", expected_elements,
    "; pair rows: ", nrow(pair_table), "\n", sep = "")
