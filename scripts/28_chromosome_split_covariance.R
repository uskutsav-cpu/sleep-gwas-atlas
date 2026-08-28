#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
fail <- function(...) stop(paste0(...), call. = FALSE)

parse_args <- function(args) {
  defaults <- list(
    panel = "config/analysis_panel.tsv",
    inclusion = "results/tables/genomicsem_trait_inclusion.tsv",
    munged_dir = "data/munged_chromosome_split/odd",
    ld_dir = "ref/eur_w_ld_chr_odd",
    split = "odd",
    out = "results/tables/genomicsem_discovery_odd.rds",
    metadata = "results/tables/genomicsem_discovery_odd_metadata.tsv",
    log_prefix = "results/logs/genomicsem/discovery_odd"
  )
  if (!length(args)) return(defaults)
  if (length(args) %% 2L) fail("Arguments must be --name value pairs")
  for (i in seq(1L, length(args), by = 2L)) {
    key <- gsub("-", "_", sub("^--", "", args[[i]]))
    if (!key %in% names(defaults)) fail("Unknown argument: ", args[[i]])
    defaults[[key]] <- args[[i + 1L]]
  }
  defaults
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!args$split %in% c("odd", "even")) fail("split must be odd or even")
panel <- read.delim(args$panel, check.names = FALSE, na.strings = c("NA", ""))
if (nrow(panel) != 45L || length(unique(panel$trait_id)) != 45L ||
    sum(panel$domain == "sleep") != 12L || unique(panel$panel_version) != "atlas-v1.0") {
  fail("Expected the exact locked atlas-v1.0 panel")
}
inclusion <- read.delim(args$inclusion, check.names = FALSE)
if (!identical(inclusion$trait_id, panel$trait_id)) fail("SEM inclusion ledger differs from panel order")
selected <- inclusion$include_genomic_sem == TRUE
if (sum(selected) != 42L) fail("Expected exactly 42 QC-passing GenomicSEM traits")
model_panel <- panel[selected, ]
trait_ids <- model_panel$trait_id
trait_files <- file.path(args$munged_dir, paste0(trait_ids, ".sumstats.gz"))
if (any(!file.exists(trait_files) | file.info(trait_files)$size <= 0)) fail("Split munged inputs are incomplete")
for (chromosome in 1:22) {
  for (suffix in c(".l2.ldscore.gz", ".l2.M_5_50")) {
    path <- file.path(args$ld_dir, paste0(chromosome, suffix))
    if (!file.exists(path) || file.info(path)$size <= 0) fail("Missing split reference: ", path)
  }
}

is_binary <- model_panel$type == "binary"
sample_prev <- rep(NA_real_, nrow(model_panel))
sample_prev[is_binary] <- model_panel$ncase[is_binary] /
  (model_panel$ncase[is_binary] + model_panel$ncontrol[is_binary])
population_prev <- rep(NA_real_, nrow(model_panel))
population_prev[is_binary] <- model_panel$pop_prev[is_binary]

suppressPackageStartupMessages(library(GenomicSEM))
if (getRversion() != "4.3.3" || as.character(packageVersion("GenomicSEM")) != "0.0.5") {
  fail("Pinned R 4.3.3 and GenomicSEM 0.0.5 are required")
}
dir.create(dirname(args$out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args$log_prefix), recursive = TRUE, showWarnings = FALSE)
started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
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
finished <- format(Sys.time(), tz = "UTC", usetz = TRUE)

n_traits <- length(trait_ids)
n_elements <- as.integer(n_traits * (n_traits + 1L) / 2L)
if (!identical(dim(covstruc$S), c(n_traits, n_traits)) ||
    !identical(dim(covstruc$V), c(n_elements, n_elements)) ||
    !identical(dim(covstruc$I), c(n_traits, n_traits)) ||
    any(!is.finite(covstruc$S)) || any(!is.finite(covstruc$V)) ||
    any(!is.finite(covstruc$I))) fail("Invalid split covariance structure")
rownames(covstruc$S) <- colnames(covstruc$S) <- trait_ids
rownames(covstruc$S_Stand) <- colnames(covstruc$S_Stand) <- trait_ids
rownames(covstruc$I) <- colnames(covstruc$I) <- trait_ids

tmp <- paste0(args$out, ".tmp")
saveRDS(covstruc, tmp, compress = "xz")
if (!file.rename(tmp, args$out)) fail("Could not atomically publish ", args$out)
metadata <- data.frame(
  field = c("panel_version", "chromosome_split", "trait_count", "covariance_elements",
            "jackknife_blocks", "r_version", "genomicsem_version", "genomicsem_commit",
            "started_at_utc", "finished_at_utc", "traits"),
  value = c("atlas-v1.0", args$split, n_traits, n_elements,
            ((n_traits + 1L) * (n_traits + 2L) / 2L) + 1L,
            as.character(getRversion()), as.character(packageVersion("GenomicSEM")),
            "6b65ca5db39fdade08b0d811477be1cdd57b5039", started, finished,
            paste(trait_ids, collapse = ","))
)
write.table(metadata, args$metadata, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Published ", args$split, "-chromosome covariance: ", n_traits, " traits and ",
    n_elements, " elements\n", sep = "")
