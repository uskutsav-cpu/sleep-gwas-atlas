#!/usr/bin/env Rscript

# Additive post-discovery continuation for the immutable Track B V1 discovery
# family.  This runner deliberately has no discovery mode and never publishes a
# canonical checkpoint link; a supervisor must validate and commit the attempt.

suppressPackageStartupMessages(library(LAVA))

schema_version <- "track-b-lava-continuation-checkpoint.1"
allowed_phases <- c("aggregate-discovery", "conditional", "finalize")

values <- list(
  input_info = "results/track_b/lava_input_info.tsv",
  sample_overlap = "results/track_b/lava_sample_overlap.txt",
  pair_manifest = "results/track_b/lava_pair_manifest.tsv",
  conditional_manifest = "results/track_b/local_conditional_manifest.tsv",
  runtime_policy = "results/track_b/lava_runtime_policy.tsv",
  phase = NA_character_,
  locus_index = NA_integer_,
  source_fingerprint = NA_character_,
  continuation_fingerprint = NA_character_,
  source_checkpoint_root = NA_character_,
  continuation_checkpoint_root = NA_character_,
  worker_output = NA_character_
)

args <- commandArgs(trailingOnly = TRUE)
mapping <- c(
  "--phase" = "phase",
  "--locus-index" = "locus_index",
  "--source-fingerprint" = "source_fingerprint",
  "--continuation-fingerprint" = "continuation_fingerprint",
  "--source-checkpoint-root" = "source_checkpoint_root",
  "--continuation-checkpoint-root" = "continuation_checkpoint_root",
  "--worker-output" = "worker_output"
)
seen <- character(0)
i <- 1L
while (i <= length(args)) {
  option <- args[[i]]
  if (!(option %in% names(mapping)) || i == length(args)) {
    stop(paste("Unknown/incomplete argument:", option))
  }
  if (option %in% seen) stop(paste("Duplicate argument:", option))
  seen <- c(seen, option)
  values[[unname(mapping[[option]])]] <- args[[i + 1L]]
  i <- i + 2L
}
values$locus_index <- suppressWarnings(as.integer(values$locus_index))

scalar_string <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}
valid_fingerprint <- function(value) {
  scalar_string(value) && grepl("^[0-9a-f]{64}$", value)
}
if (!scalar_string(values$phase) || !(values$phase %in% allowed_phases)) {
  stop("--phase must be aggregate-discovery, conditional, or finalize; discovery is forbidden in V2")
}
if (!valid_fingerprint(values$source_fingerprint)) stop("Invalid --source-fingerprint")
if (!valid_fingerprint(values$continuation_fingerprint)) stop("Invalid --continuation-fingerprint")
if (identical(values$source_fingerprint, values$continuation_fingerprint)) {
  stop("Source and continuation fingerprints must differ")
}
if (values$phase == "conditional") {
  if (is.na(values$locus_index) || values$locus_index < 1L || values$locus_index > 2495L) {
    stop("Conditional phase requires --locus-index in 1..2495")
  }
} else if (!is.na(values$locus_index)) {
  stop("Aggregate/finalize phases do not accept --locus-index")
}

repository_root <- normalizePath(".", mustWork = TRUE)
path_is_symlink <- function(path) {
  target <- Sys.readlink(path)
  length(target) == 1L && !is.na(target) && nzchar(target)
}
direct_directory <- function(path, label) {
  if (!scalar_string(path) || !dir.exists(path)) stop(paste(label, "does not exist"))
  expanded <- sub("/+$", "", path.expand(path))
  direct <- file.path(
    normalizePath(dirname(expanded), mustWork = TRUE), basename(expanded)
  )
  if (path_is_symlink(direct)) stop(paste(label, "must not be a symbolic link"))
  normalized <- normalizePath(direct, mustWork = TRUE)
  if (!identical(normalized, direct)) {
    stop(paste(label, "must be a direct directory path without symlink traversal"))
  }
  if (grepl("[\r\n\t]", normalized)) {
    stop(paste(label, "contains a forbidden control character"))
  }
  normalized
}
normalize_namespace_root <- function(path, fingerprint, expected_parent, label) {
  normalized <- direct_directory(path, label)
  expected <- file.path(repository_root, expected_parent, fingerprint)
  if (!identical(normalized, expected) || !identical(basename(normalized), fingerprint)) {
    stop(paste(label, "must be the exact direct fingerprint namespace root"))
  }
  normalized
}
source_root <- normalize_namespace_root(
  values$source_checkpoint_root, values$source_fingerprint,
  "results/track_b/checkpoints/lava", "Source checkpoint root"
)
continuation_root <- normalize_namespace_root(
  values$continuation_checkpoint_root, values$continuation_fingerprint,
  "results/track_b/checkpoints/lava_continuations/runs",
  "Continuation checkpoint root"
)
if (identical(source_root, continuation_root)) stop("Source and continuation roots must differ")

if (!scalar_string(values$worker_output)) stop("--worker-output is required")
if (!identical(basename(values$worker_output), "result.rds")) {
  stop("Worker output basename must be result.rds")
}
attempt_root <- direct_directory(
  file.path(continuation_root, ".attempts"), "Continuation attempt root"
)
if (!identical(attempt_root, file.path(continuation_root, ".attempts"))) {
  stop("Continuation attempt root must be directly inside its namespace")
}
worker_parent <- direct_directory(
  dirname(values$worker_output), "Worker attempt directory"
)
if (!identical(dirname(worker_parent), attempt_root)) {
  stop("Worker output must be inside one supervised continuation attempt directory")
}
values$worker_output <- file.path(worker_parent, "result.rds")
if (path_is_symlink(values$worker_output)) {
  stop("Worker attempt result must not be a symbolic link")
}
if (file.exists(values$worker_output)) {
  stop("Worker attempt result already exists; overwrite is forbidden")
}

canonical_checkpoint_path <- function(phase, index = NA_integer_) {
  unit <- if (phase %in% c("aggregate-discovery", "finalize")) {
    "unit_all"
  } else {
    sprintf("locus_%04d", index)
  }
  file.path(continuation_root, phase, unit, "result.rds")
}
if (file.exists(canonical_checkpoint_path(values$phase, values$locus_index))) {
  stop("Supervisor invoked a continuation worker whose canonical ready bundle already exists")
}

if (as.character(packageVersion("LAVA")) != "0.1.5") stop("LAVA 0.1.5 is required")
runtime_check <- suppressWarnings(system2(
  ".r-env/bin/Rscript", "scripts/133_validate_track_b_lava_runtime.R",
  stdout = TRUE, stderr = TRUE
))
if (!is.null(attr(runtime_check, "status")) && attr(runtime_check, "status") != 0L) {
  stop(paste("Pinned LAVA implementation failed verification:", paste(runtime_check, collapse = " | ")))
}
if (sum(grepl("^TRACK_B_LAVA_RUNTIME\\tstatus=PASS\\t", runtime_check)) != 1L) {
  stop("Pinned LAVA runtime did not emit one semantic PASS marker")
}

required <- c(
  values$input_info, values$sample_overlap, values$pair_manifest,
  values$conditional_manifest, values$runtime_policy
)
if (!all(file.exists(required))) {
  stop(paste("Missing Track B LAVA input:", paste(required[!file.exists(required)], collapse = ", ")))
}

policy_table <- read.delim(values$runtime_policy, stringsAsFactors = FALSE, check.names = FALSE)
policy <- setNames(as.list(policy_table$value), policy_table$key)
policy_num <- function(key) as.numeric(policy[[key]])
policy_int <- function(key) as.integer(policy[[key]])
policy_bool <- function(key) identical(tolower(policy[[key]]), "true")
if (policy$analysis_id != "track-b-v1.0-local" || policy$lava_version != "0.1.5") {
  stop("Track B runtime policy drifted")
}
if (policy$ram_execution_unit != "WHOLE_PREDECLARED_LAVA_LOCUS" ||
    policy_int("maximum_loci_per_process") != 1L) {
  stop("Track B RAM-aware execution unit drifted")
}

input_info <- read.delim(values$input_info, stringsAsFactors = FALSE, check.names = FALSE)
pairs <- read.delim(values$pair_manifest, stringsAsFactors = FALSE, check.names = FALSE)
conditioners <- read.delim(values$conditional_manifest, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(input_info) != policy_int("expected_traits") ||
    length(unique(input_info$phenotype)) != policy_int("expected_traits")) {
  stop("Input is not the exact eight-trait family")
}
if (nrow(pairs) != policy_int("expected_pairs") ||
    !identical(as.character(pairs$pair_id), c("A", "B", "CONTROL"))) {
  stop("Pair manifest is not A/B/CONTROL")
}
if (any(conditioners$selection_timing != "BEFORE_LOCAL_RESULT_ACCESS")) {
  stop("Conditional covariates were not frozen before results")
}
if (!identical(
  as.character(conditioners$conditional_model_id),
  c("A_BMI_ONLY", "A_SLEEP_APNEA_ONLY", "B_MDD_ONLY", "NONE")
)) {
  stop("Conditional model family drifted")
}

locus_file <- "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
loci <- read.loci(locus_file)
if (nrow(loci) != policy_int("expected_loci") || anyDuplicated(loci$LOC)) {
  stop("Locus file is not the exact 2,495 family")
}

path_within <- function(path, root) {
  resolved <- normalizePath(path, mustWork = TRUE)
  identical(resolved, root) || startsWith(paste0(resolved, "/"), paste0(root, "/"))
}

sha256_files <- function(paths) {
  if (!length(paths)) return(character(0))
  if (any(!file.exists(paths))) stop("Cannot hash a missing dependency")
  hashes <- character(length(paths))
  groups <- split(seq_along(paths), ceiling(seq_along(paths) / 64L))
  for (indices in groups) {
    output <- suppressWarnings(system2(
      "shasum",
      c("-a", "256", vapply(paths[indices], shQuote, character(1L))),
      stdout = TRUE, stderr = TRUE
    ))
    status <- attr(output, "status")
    if ((!is.null(status) && status != 0L) || length(output) != length(indices) ||
        any(!grepl("^[0-9a-f]{64}[[:space:]]", output))) {
      stop(paste("SHA-256 dependency hashing failed:", paste(output, collapse = " | ")))
    }
    hashes[indices] <- substr(output, 1L, 64L)
  }
  hashes
}

read_rds_with_sha256 <- function(path, label) {
  if (!file.exists(path) || file.info(path)$size <= 0L) {
    stop(paste(label, "is missing or empty"))
  }
  resolved_before <- normalizePath(path, mustWork = TRUE)
  before <- file.info(path)[1L, c("size", "mode", "mtime", "ctime")]
  connection <- file(path, open = "rb")
  raw <- tryCatch(
    readBin(connection, what = "raw", n = as.numeric(before$size) + 1),
    finally = close(connection)
  )
  after <- file.info(path)[1L, c("size", "mode", "mtime", "ctime")]
  resolved_after <- normalizePath(path, mustWork = TRUE)
  if (!identical(resolved_before, resolved_after) || !identical(before, after) ||
      length(raw) != as.numeric(before$size)) {
    stop(paste(label, "changed while being read"))
  }
  temporary <- tempfile(pattern = ".track-b-lava-rds-bytes-")
  on.exit(unlink(temporary), add = TRUE)
  writeBin(raw, temporary)
  digest <- sha256_files(temporary)[[1L]]
  payload <- suppressWarnings(memDecompress(raw, type = "unknown"))
  value <- tryCatch(
    unserialize(payload),
    error = function(error) stop(paste(label, "is not a valid RDS:", conditionMessage(error)))
  )
  list(value = value, sha256 = digest)
}

source_relative_paths <- file.path(
  "discovery", sprintf("locus_%04d", seq_len(nrow(loci))), "result.rds"
)
source_paths <- file.path(source_root, source_relative_paths)

status_pre_fields <- c(
  "LOC", "CHR", "START", "STOP", "status", "n_snps", "K",
  "univariate_tested", "eligible_bivariate_pairs", "elapsed_seconds"
)
univ_pre_fields <- c(
  "LOC", "CHR", "START", "STOP", "phen", "h2.obs", "h2.latent",
  "ascertained", "p", "analysis_status", "error", "n_snps", "K"
)
bivar_pre_fields <- c(
  "LOC", "CHR", "START", "STOP", "pair_id", "trait1", "trait2",
  "discovery_rg", "discovery_SE", "discovery_P", "discovery_FDR",
  "local_covariance", "rho", "rho.lower", "rho.upper", "r2", "r2.lower",
  "r2.upper", "p", "analysis_status", "error"
)
conditional_pre_fields <- c(
  "LOC", "CHR", "START", "STOP", "pair_id", "trait1", "trait2",
  "conditional_model_id", "covariates", "pcor", "ci.lower", "ci.upper",
  "p", "r2.trait1_z", "r2.trait2_z", "analysis_status", "error"
)

exact_frame_names <- function(value, expected, label) {
  if (!is.data.frame(value) || !identical(names(value), expected)) {
    stop(paste(label, "schema drifted"))
  }
}
same_coordinates <- function(table, row) {
  if (!nrow(table)) return(TRUE)
  all(
    as.character(table$LOC) == as.character(row$LOC),
    as.character(table$CHR) == as.character(row$CHR),
    as.character(table$START) == as.character(row$START),
    as.character(table$STOP) == as.character(row$STOP)
  )
}
finite_or_na_vector <- function(value) all(is.na(value) | is.finite(value))
finite_or_na_scalar <- function(value) {
  if (length(value) == 1L && is.finite(value)) as.numeric(value) else NA_real_
}

validate_source_discovery <- function(checkpoint, index) {
  expected_names <- c("fingerprint", "phase", "locus_index", "LOC", "CHR", "status", "univ", "bivar")
  official <- loci[index, , drop = FALSE]
  if (!is.list(checkpoint) || !identical(names(checkpoint), expected_names) ||
      !identical(checkpoint$fingerprint, values$source_fingerprint) ||
      !identical(checkpoint$phase, "discovery") ||
      !identical(as.integer(checkpoint$locus_index), as.integer(index)) ||
      !identical(as.character(checkpoint$LOC), as.character(official$LOC)) ||
      !identical(as.character(checkpoint$CHR), as.character(official$CHR))) {
    stop(paste("Stale or malformed V1 discovery checkpoint for locus index", index))
  }
  exact_frame_names(checkpoint$status, status_pre_fields, "Source discovery status")
  exact_frame_names(checkpoint$univ, univ_pre_fields, "Source discovery univariate")
  exact_frame_names(checkpoint$bivar, bivar_pre_fields, "Source discovery bivariate")
  if (nrow(checkpoint$status) != 1L || !same_coordinates(checkpoint$status, official)) {
    stop(paste("Source discovery status coordinates drifted at locus index", index))
  }
  if (nrow(checkpoint$univ) != policy_int("expected_traits") ||
      !identical(as.character(checkpoint$univ$phen), as.character(input_info$phenotype)) ||
      !same_coordinates(checkpoint$univ, official)) {
    stop(paste("Source discovery univariate family drifted at locus index", index))
  }
  if (!(checkpoint$status$status %in% c("PROCESSED", "PROCESS_FAILED", "UNIVARIATE_FAILED")) ||
      any(!(checkpoint$univ$analysis_status %in%
        c("TESTED", "PHENOTYPE_DROPPED", "LOCUS_PROCESS_FAILED", "UNIVARIATE_FAILED")))) {
    stop(paste("Source discovery terminal status drifted at locus index", index))
  }
  if (any(!(checkpoint$bivar$pair_id %in% pairs$pair_id)) ||
      anyDuplicated(checkpoint$bivar$pair_id) || !same_coordinates(checkpoint$bivar, official)) {
    stop(paste("Source discovery bivariate family drifted at locus index", index))
  }
  pair_index <- match(checkpoint$bivar$pair_id, pairs$pair_id)
  if (nrow(checkpoint$bivar) &&
      any(checkpoint$bivar$trait1 != pairs$trait1[pair_index] |
          checkpoint$bivar$trait2 != pairs$trait2[pair_index])) {
    stop(paste("Source discovery bivariate traits drifted at locus index", index))
  }
  if (any(!(checkpoint$bivar$analysis_status %in% c("TESTED", "BIVARIATE_FAILED"))) ||
      !finite_or_na_vector(checkpoint$univ$p) ||
      !finite_or_na_vector(checkpoint$bivar$p)) {
    stop(paste("Source discovery numerical/status fields drifted at locus index", index))
  }
  tested_univ <- checkpoint$univ$analysis_status == "TESTED"
  eligible_traits <- checkpoint$univ$phen[
    tested_univ & is.finite(checkpoint$univ$p) &
      checkpoint$univ$p <= policy_num("univariate_p_threshold")
  ]
  expected_pairs <- pairs$pair_id[
    pairs$trait1 %in% eligible_traits & pairs$trait2 %in% eligible_traits
  ]
  if (!identical(as.character(checkpoint$bivar$pair_id), as.character(expected_pairs))) {
    stop(paste("Source discovery bivariate rows differ from the local-h2 gate at locus index", index))
  }
  checkpoint
}

read_source_family <- function() {
  if (!all(file.exists(source_paths))) {
    stop(sprintf(
      "Incomplete V1 discovery checkpoint family: %d/%d",
      sum(file.exists(source_paths)), length(source_paths)
    ))
  }
  if (any(!vapply(source_paths, path_within, logical(1L), root = source_root))) {
    stop("A V1 discovery result resolves outside the source checkpoint root")
  }
  sealed <- lapply(seq_along(source_paths), function(index) {
    observed <- read_rds_with_sha256(
      source_paths[[index]], paste("Source discovery checkpoint", index)
    )
    observed$value <- validate_source_discovery(observed$value, index)
    observed
  })
  manifest <- data.frame(
    locus_index = seq_len(nrow(loci)),
    LOC = as.character(loci$LOC),
    relative_path = source_relative_paths,
    sha256 = vapply(sealed, function(item) item$sha256, character(1L)),
    stringsAsFactors = FALSE
  )
  results <- lapply(sealed, function(item) item$value)
  list(results = results, manifest = manifest)
}

empty_conditional <- function() data.frame(
  LOC = numeric(0), CHR = numeric(0), START = numeric(0), STOP = numeric(0),
  pair_id = character(0), trait1 = character(0), trait2 = character(0),
  conditional_model_id = character(0), covariates = character(0),
  pcor = numeric(0), ci.lower = numeric(0), ci.upper = numeric(0),
  p = numeric(0), r2.trait1_z = numeric(0), r2.trait2_z = numeric(0),
  analysis_status = character(0), error = character(0),
  stringsAsFactors = FALSE
)

build_aggregate <- function(source_family) {
  status <- do.call(rbind, lapply(source_family$results, function(value) value$status))
  univ <- do.call(rbind, lapply(source_family$results, function(value) value$univ))
  bivar <- do.call(rbind, lapply(source_family$results, function(value) value$bivar))
  if (nrow(status) != policy_int("expected_loci") ||
      nrow(univ) != policy_int("planned_univariate_tests")) {
    stop("Incomplete V1 discovery checkpoint contents")
  }

  planned_univ <- policy_int("planned_univariate_tests")
  univ$univariate_test_family_n <- planned_univ
  univ$univariate_p_threshold <- policy_num("univariate_p_threshold")
  univ$p_bonferroni <- pmin(1, univ$p * planned_univ)
  univ$p_fdr <- NA_real_
  tested_univ <- which(univ$analysis_status == "TESTED" & is.finite(univ$p))
  univ_family_p <- rep.int(1, nrow(univ))
  univ_family_p[tested_univ] <- univ$p[tested_univ]
  univ_adjusted <- p.adjust(univ_family_p, method = "BH")
  univ$p_fdr[tested_univ] <- univ_adjusted[tested_univ]

  bivar$bivariate_test_family_n <- rep.int(nrow(bivar), nrow(bivar))
  bivar$p_fdr <- rep(NA_real_, nrow(bivar))
  tested_bivar <- which(bivar$analysis_status == "TESTED" & is.finite(bivar$p))
  bivar_family_p <- rep.int(1, nrow(bivar))
  bivar_family_p[tested_bivar] <- bivar$p[tested_bivar]
  bivar_adjusted <- p.adjust(bivar_family_p, method = "BH")
  bivar$p_fdr[tested_bivar] <- bivar_adjusted[tested_bivar]
  bivar$fdr_significant <- !is.na(bivar$p_fdr) &
    bivar$p_fdr <= policy_num("bivariate_fdr_alpha")

  status$analysis_fingerprint <- values$source_fingerprint
  status$lava_version <- as.character(packageVersion("LAVA"))
  status$reference_prefix <- policy$reference_prefix

  list(
    schema_version = schema_version,
    source_discovery_fingerprint = values$source_fingerprint,
    continuation_execution_fingerprint = values$continuation_fingerprint,
    phase = "aggregate-discovery",
    source_checkpoint_manifest = source_family$manifest,
    status = status,
    univ = univ,
    bivar = bivar
  )
}

build_candidates <- function(family) {
  candidate_rows <- list()
  candidate_index <- 0L
  for (index in seq_len(nrow(loci))) {
    locus_row <- loci[index, , drop = FALSE]
    targets <- family$bivar[
      as.character(family$bivar$LOC) == as.character(locus_row$LOC) &
        family$bivar$fdr_significant &
        family$bivar$pair_id %in% c("A", "B"),
      , drop = FALSE
    ]
    eligible_models <- 0L
    eligible_pairs <- character(0)
    if (nrow(targets)) for (target_index in seq_len(nrow(targets))) {
      target <- targets[target_index, , drop = FALSE]
      models <- conditioners[
        conditioners$pair_id == target$pair_id &
          conditioners$conditional_model_id != "NONE",
        , drop = FALSE
      ]
      if (nrow(models)) for (model_index in seq_len(nrow(models))) {
        covars <- strsplit(models$covariates[[model_index]], ";", fixed = TRUE)[[1L]]
        local_univ <- family$univ[
          as.character(family$univ$LOC) == as.character(locus_row$LOC) &
            family$univ$phen %in% covars,
          , drop = FALSE
        ]
        eligible <- length(covars) > 0L && nrow(local_univ) == length(covars) &&
          all(
            local_univ$analysis_status == "TESTED" &
              is.finite(local_univ$p) &
              local_univ$p <= policy_num("univariate_p_threshold")
          )
        if (eligible) {
          eligible_models <- eligible_models + 1L
          eligible_pairs <- unique(c(eligible_pairs, as.character(target$pair_id)))
        }
      }
    }
    if (eligible_models > 0L) {
      candidate_index <- candidate_index + 1L
      candidate_rows[[candidate_index]] <- data.frame(
        locus_index = index,
        locus = locus_row$LOC,
        chromosome = locus_row$CHR,
        start = locus_row$START,
        stop = locus_row$STOP,
        n_snps = family$status$n_snps[
          match(as.character(locus_row$LOC), as.character(family$status$LOC))
        ],
        eligible_models = eligible_models,
        eligible_pairs = paste(eligible_pairs, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(candidate_rows)) {
    do.call(rbind, candidate_rows)
  } else {
    data.frame(
      locus_index = integer(0), locus = integer(0), chromosome = integer(0),
      start = integer(0), stop = integer(0), n_snps = integer(0),
      eligible_models = integer(0), eligible_pairs = character(0),
      stringsAsFactors = FALSE
    )
  }
}

candidate_fields <- c(
  "locus_index", "locus", "chromosome", "start", "stop", "n_snps",
  "eligible_models", "eligible_pairs"
)

canonical_candidates <- function(table, label) {
  if (!is.data.frame(table) || !identical(names(table), candidate_fields)) {
    stop(paste(label, "schema drifted"))
  }
  integer_fields <- candidate_fields[candidate_fields != "eligible_pairs"]
  if (nrow(table) && any(vapply(
    table[integer_fields],
    function(column) any(!is.finite(column) | column != as.integer(column)),
    logical(1L)
  ))) {
    stop(paste(label, "contains a non-integer or non-finite identity"))
  }
  canonical <- data.frame(
    locus_index = as.integer(table$locus_index),
    locus = as.integer(table$locus),
    chromosome = as.integer(table$chromosome),
    start = as.integer(table$start),
    stop = as.integer(table$stop),
    n_snps = as.integer(table$n_snps),
    eligible_models = as.integer(table$eligible_models),
    eligible_pairs = as.character(table$eligible_pairs),
    stringsAsFactors = FALSE
  )
  if (nrow(canonical) && (
    any(canonical$locus_index < 1L | canonical$locus_index > nrow(loci)) ||
      any(canonical$n_snps < 1L) || any(canonical$eligible_models < 1L) ||
      any(!(canonical$eligible_pairs %in% c("A", "B", "A;B"))) ||
      is.unsorted(canonical$locus_index, strictly = TRUE)
  )) {
    stop(paste(label, "contains an invalid or duplicated candidate identity"))
  }
  canonical
}

read_conditional_candidates <- function(path, family) {
  if (!file.exists(path) || file.info(path)$size <= 0L) {
    stop("Aggregate checkpoint lacks its conditional candidate attestation")
  }
  observed <- read.delim(
    path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA",
    colClasses = c(rep("integer", 7L), "character")
  )
  observed <- canonical_candidates(observed, "Conditional candidate attestation")
  expected <- canonical_candidates(
    build_candidates(family), "Recomputed conditional candidates"
  )
  if (!identical(observed, expected)) {
    stop("Conditional candidate attestation differs from the sealed aggregate gates")
  }
  observed
}

publish_rds_no_replace <- function(value, path) {
  if (file.exists(path)) stop(paste("Checkpoint already exists; overwrite is forbidden:", path))
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(value, temporary, version = 3)
  if (!file.link(temporary, path)) {
    stop(paste("Could not atomically publish checkpoint without replacement:", path))
  }
  invisible(path)
}

write_table_no_replace <- function(table, path) {
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  write.table(table, temporary, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  if (file.exists(path)) {
    if (!identical(unname(tools::md5sum(temporary)), unname(tools::md5sum(path)))) {
      stop(paste("Existing output differs; overwrite is forbidden:", path))
    }
    return(invisible(path))
  }
  if (!file.link(temporary, path)) {
    stop(paste("Could not atomically publish output without replacement:", path))
  }
  invisible(path)
}

write_lines_no_replace <- function(lines, path) {
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  writeLines(lines, temporary, useBytes = TRUE)
  if (file.exists(path)) {
    if (!identical(unname(tools::md5sum(temporary)), unname(tools::md5sum(path)))) {
      stop(paste("Existing attestation differs; overwrite is forbidden:", path))
    }
    return(invisible(path))
  }
  if (!file.link(temporary, path)) {
    stop(paste("Could not atomically publish attestation without replacement:", path))
  }
  invisible(path)
}

worker_marker <- function(
  phase, index, locus, chromosome, n_snps, pair, qc,
  construction_complete, output
) {
  marker_n_snps <- if (length(n_snps) == 1L && is.finite(n_snps)) {
    as.character(as.integer(n_snps))
  } else {
    "NA"
  }
  cat(paste(
    "TRACK_B_LAVA_CONTINUATION_WORKER",
    paste0("phase=", phase),
    paste0("index=", index),
    paste0("locus=", locus),
    paste0("chromosome=", chromosome),
    paste0("n_snps=", marker_n_snps),
    paste0("pair=", pair),
    paste0("qc=", qc),
    paste0("construction_complete=", construction_complete),
    paste0("output=", output),
    paste0("source_fingerprint=", values$source_fingerprint),
    paste0("continuation_fingerprint=", values$continuation_fingerprint),
    sep = "\t"
  ), "\n", sep = "")
}

validate_aggregate_identity <- function(family) {
  expected_names <- c(
    "schema_version", "source_discovery_fingerprint",
    "continuation_execution_fingerprint", "phase",
    "source_checkpoint_manifest", "status", "univ", "bivar"
  )
  if (!is.list(family) || !identical(names(family), expected_names) ||
      !identical(family$schema_version, schema_version) ||
      !identical(family$source_discovery_fingerprint, values$source_fingerprint) ||
      !identical(family$continuation_execution_fingerprint, values$continuation_fingerprint) ||
      !identical(family$phase, "aggregate-discovery")) {
    stop("Continuation aggregate checkpoint identity drifted")
  }
  manifest_fields <- c("locus_index", "LOC", "relative_path", "sha256")
  if (!is.data.frame(family$source_checkpoint_manifest) ||
      !identical(names(family$source_checkpoint_manifest), manifest_fields) ||
      nrow(family$source_checkpoint_manifest) != nrow(loci) ||
      !identical(as.integer(family$source_checkpoint_manifest$locus_index), seq_len(nrow(loci))) ||
      !identical(as.character(family$source_checkpoint_manifest$LOC), as.character(loci$LOC)) ||
      !identical(as.character(family$source_checkpoint_manifest$relative_path), source_relative_paths) ||
      any(!grepl("^[0-9a-f]{64}$", family$source_checkpoint_manifest$sha256))) {
    stop("Continuation aggregate source manifest drifted")
  }
  family
}

if (values$phase == "aggregate-discovery") {
  source_family <- read_source_family()
  aggregate <- build_aggregate(source_family)
  candidates <- build_candidates(aggregate)
  write_table_no_replace(
    candidates, file.path(dirname(values$worker_output), "conditional_candidates.tsv")
  )
  publish_rds_no_replace(aggregate, values$worker_output)
  total_snps <- sum(aggregate$status$n_snps, na.rm = TRUE)
  rm(source_family, aggregate, candidates)
  invisible(gc(full = TRUE))
  worker_marker(
    "aggregate-discovery", 0L, "ALL", "ALL", total_snps,
    "A;B;CONTROL", "FULL_FAMILY_BH_COMPLETE", "TRUE", values$worker_output
  )
  quit(save = "no", status = 0L)
}

aggregate_path <- canonical_checkpoint_path("aggregate-discovery")
if (!file.exists(aggregate_path) || file.info(aggregate_path)$size <= 0L ||
    !path_within(aggregate_path, continuation_root)) {
  stop("Canonical continuation aggregate checkpoint is missing or escapes its namespace")
}
aggregate_read <- read_rds_with_sha256(aggregate_path, "Continuation aggregate checkpoint")
sealed_aggregate <- validate_aggregate_identity(aggregate_read$value)
aggregate_sha256 <- aggregate_read$sha256
rm(aggregate_read)

scoped_input <- function(locus_row) {
  chromosome <- as.integer(locus_row$CHR)
  chromosome_prefix <- paste0(policy$reference_prefix, "_chr", chromosome)
  reference_files <- paste0(chromosome_prefix, c(".info", ".bcor"))
  if (!all(file.exists(reference_files))) {
    stop(paste(
      "Eligible conditional locus has an incomplete LAVA reference for chromosome",
      chromosome, "missing",
      paste(reference_files[!file.exists(reference_files)], collapse = ", ")
    ))
  }
  chromosome_input_info <- file.path(
    "results/track_b/lava_chromosome_inputs",
    sprintf("chr%02d", chromosome),
    "lava_input_info.tsv"
  )
  if (!file.exists(chromosome_input_info)) {
    stop(paste("Verified chromosome input-info is missing:", chromosome_input_info))
  }
  process.input(
    input.info.file = chromosome_input_info,
    sample.overlap.file = values$sample_overlap,
    ref.prefix = chromosome_prefix,
    phenos = input_info$phenotype
  )
}

process_conditional_locus <- function(index, family, candidate_attested) {
  locus_row <- loci[index, , drop = FALSE]
  targets <- family$bivar[
    as.character(family$bivar$LOC) == as.character(locus_row$LOC) &
      family$bivar$fdr_significant &
      family$bivar$pair_id %in% c("A", "B"),
    , drop = FALSE
  ]
  conditional_rows <- list()
  conditional_index <- 0L
  locus <- NULL
  input <- NULL
  locus_process_attempted <- FALSE
  locus_process_ok <- FALSE

  if (nrow(targets)) for (target_index in seq_len(nrow(targets))) {
    target <- targets[target_index, , drop = FALSE]
    models <- conditioners[
      conditioners$pair_id == target$pair_id &
        conditioners$conditional_model_id != "NONE",
      , drop = FALSE
    ]
    if (!nrow(models)) stop(paste("No frozen conditional model for pair", target$pair_id))
    for (model_index in seq_len(nrow(models))) {
      conditional_index <- conditional_index + 1L
      model <- models[model_index, , drop = FALSE]
      covars <- strsplit(model$covariates, ";", fixed = TRUE)[[1L]]
      base <- data.frame(
        LOC = target$LOC, CHR = target$CHR, START = target$START, STOP = target$STOP,
        pair_id = target$pair_id, trait1 = target$trait1, trait2 = target$trait2,
        conditional_model_id = model$conditional_model_id,
        covariates = model$covariates,
        stringsAsFactors = FALSE
      )
      local_univ <- family$univ[
        as.character(family$univ$LOC) == as.character(target$LOC) &
          family$univ$phen %in% covars,
        , drop = FALSE
      ]
      eligible <- length(covars) > 0L && nrow(local_univ) == length(covars) &&
        all(
          local_univ$analysis_status == "TESTED" &
            is.finite(local_univ$p) &
            local_univ$p <= policy_num("univariate_p_threshold")
        )
      if (!eligible) {
        conditional_rows[[conditional_index]] <- cbind(
          base,
          pcor = NA_real_, ci.lower = NA_real_, ci.upper = NA_real_, p = NA_real_,
          r2.trait1_z = NA_real_, r2.trait2_z = NA_real_,
          analysis_status = "CONDITIONER_LOCAL_H2_INELIGIBLE",
          error = "At least one predeclared conditioner did not pass the same-locus univariate gate"
        )
        next
      }

      if (!isTRUE(candidate_attested)) {
        stop("An eligible conditional model is absent from the sealed candidate attestation")
      }

      if (is.null(input)) input <- scoped_input(locus_row)
      if (is.null(locus)) {
        locus_process_attempted <- TRUE
        locus <- tryCatch(
          process.locus(
            locus_row, input,
            min.K = policy_int("min_K"),
            prune.thresh = policy_num("prune_threshold"),
            max.prop.K = policy_num("max_proportion_K"),
            drop.failed = TRUE,
            max.block.size = policy_int("max_block_size"),
            cap.estimates = policy_bool("cap_estimates")
          ),
          error = function(error) {
            structure(list(message = conditionMessage(error)), class = "conditional_error")
          }
        )
        locus_process_ok <- !inherits(locus, "conditional_error") && !is.null(locus)
      }
      model_seed <- match(
        model$conditional_model_id,
        c("A_BMI_ONLY", "A_SLEEP_APNEA_ONLY", "B_MDD_ONLY")
      )
      set.seed(policy_int("random_seed") + index * 1000L + 500L + model_seed)
      partial <- if (inherits(locus, "conditional_error") || is.null(locus)) {
        locus
      } else {
        tryCatch(
          run.pcor(
            locus,
            target = c(target$trait1, target$trait2),
            phenos = c(target$trait1, target$trait2, covars),
            max.r2 = policy_num("conditional_max_r2")
          ),
          error = function(error) {
            structure(list(message = conditionMessage(error)), class = "conditional_error")
          }
        )
      }
      unique_partial <- !inherits(partial, "conditional_error") &&
        !is.null(partial) && nrow(partial) == 1L
      partial_values <- if (unique_partial) {
        unlist(partial[1L, c(
          "pcor", "ci.lower", "ci.upper", "p", "r2.phen1_z", "r2.phen2_z"
        )])
      } else {
        rep(NA_real_, 6L)
      }
      pcor_valid <- unique_partial && is.finite(partial[["pcor"]]) &&
        partial[["pcor"]] >= -1 && partial[["pcor"]] <= 1
      ci_valid <- unique_partial &&
        is.finite(partial[["ci.lower"]]) && is.finite(partial[["ci.upper"]]) &&
        partial[["ci.lower"]] <= partial[["pcor"]] &&
        partial[["pcor"]] <= partial[["ci.upper"]]
      p_valid <- unique_partial && is.finite(partial[["p"]]) &&
        partial[["p"]] >= 0 && partial[["p"]] <= 1
      r2_in_domain <- unique_partial &&
        is.finite(partial[["r2.phen1_z"]]) && is.finite(partial[["r2.phen2_z"]]) &&
        partial[["r2.phen1_z"]] >= 0 && partial[["r2.phen1_z"]] <= 1 &&
        partial[["r2.phen2_z"]] >= 0 && partial[["r2.phen2_z"]] <= 1
      r2_stable <- r2_in_domain &&
        partial[["r2.phen1_z"]] < policy_num("conditional_max_r2") &&
        partial[["r2.phen2_z"]] < policy_num("conditional_max_r2")
      valid_partial <- unique_partial && all(is.finite(partial_values)) &&
        pcor_valid && ci_valid && p_valid && r2_stable
      unstable_max_r2 <- unique_partial &&
        pcor_valid && ci_valid && r2_in_domain &&
        (
          partial[["r2.phen1_z"]] >= policy_num("conditional_max_r2") ||
            partial[["r2.phen2_z"]] >= policy_num("conditional_max_r2")
        )
      if (!valid_partial) {
        message_text <- if (inherits(partial, "conditional_error")) {
          partial$message
        } else if (!unique_partial) {
          "run.pcor returned no unique result"
        } else if (unstable_max_r2) {
          "run.pcor exceeded the predeclared conditional max-r2 stability threshold"
        } else {
          paste(
            "run.pcor returned an incomplete or invalid partial estimate,",
            "confidence interval, r2, or p-value; raw diagnostics:",
            paste(
              paste0(
                c("pcor", "ci.lower", "ci.upper", "p", "r2.phen1_z", "r2.phen2_z"),
                "=", vapply(partial_values, as.character, character(1L))
              ),
              collapse = ","
            )
          )
        }
        terminal_status <- if (unstable_max_r2) {
          "CONDITIONAL_UNSTABLE_MAX_R2"
        } else {
          "CONDITIONAL_FAILED"
        }
        conditional_rows[[conditional_index]] <- cbind(
          base,
          pcor = if (pcor_valid && ci_valid) as.numeric(partial$pcor) else NA_real_,
          ci.lower = if (pcor_valid && ci_valid) as.numeric(partial$ci.lower) else NA_real_,
          ci.upper = if (pcor_valid && ci_valid) as.numeric(partial$ci.upper) else NA_real_,
          p = if (p_valid && !unstable_max_r2) {
            finite_or_na_scalar(partial$p)
          } else {
            NA_real_
          },
          r2.trait1_z = if (r2_in_domain) {
            finite_or_na_scalar(partial[["r2.phen1_z"]])
          } else {
            NA_real_
          },
          r2.trait2_z = if (r2_in_domain) {
            finite_or_na_scalar(partial[["r2.phen2_z"]])
          } else {
            NA_real_
          },
          analysis_status = terminal_status,
          error = message_text
        )
      } else {
        conditional_rows[[conditional_index]] <- cbind(
          base,
          pcor = partial$pcor,
          ci.lower = partial$ci.lower,
          ci.upper = partial$ci.upper,
          p = partial$p,
          r2.trait1_z = partial[["r2.phen1_z"]],
          r2.trait2_z = partial[["r2.phen2_z"]],
          analysis_status = "TESTED",
          error = NA_character_
        )
      }
    }
  }

  conditional <- if (length(conditional_rows)) {
    do.call(rbind, conditional_rows)
  } else {
    empty_conditional()
  }
  state <- if (!nrow(targets)) {
    "NOT_BIVARIATE_FDR_ELIGIBLE"
  } else if (!nrow(conditional)) {
    "NO_FROZEN_MODEL"
  } else if (!locus_process_attempted) {
    "CONDITIONER_LOCAL_H2_INELIGIBLE"
  } else if (locus_process_ok) {
    "CONDITIONAL_LOCUS_PROCESSED"
  } else {
    "CONDITIONAL_LOCUS_PROCESS_FAILED"
  }
  if (!identical(isTRUE(candidate_attested), locus_process_attempted)) {
    stop("Conditional process-locus use differs from the sealed candidate attestation")
  }
  n_snps <- family$status$n_snps[
    match(as.character(locus_row$LOC), as.character(family$status$LOC))
  ]
  if (!is.null(locus)) rm(locus)
  if (!is.null(input)) rm(input)
  invisible(gc(full = TRUE))
  list(
    schema_version = schema_version,
    source_discovery_fingerprint = values$source_fingerprint,
    continuation_execution_fingerprint = values$continuation_fingerprint,
    phase = "conditional",
    locus_index = as.integer(index),
    LOC = as.character(locus_row$LOC),
    CHR = as.integer(locus_row$CHR),
    aggregate_checkpoint_sha256 = aggregate_sha256,
    n_snps = as.integer(n_snps),
    pair = if (nrow(targets)) paste(targets$pair_id, collapse = ";") else "NONE",
    state = state,
    conditional = conditional
  )
}

if (values$phase == "conditional") {
  candidates <- read_conditional_candidates(
    file.path(dirname(aggregate_path), "conditional_candidates.tsv"),
    sealed_aggregate
  )
  candidate_attested <- values$locus_index %in% candidates$locus_index
  output <- process_conditional_locus(
    values$locus_index, sealed_aggregate, candidate_attested
  )
  publish_rds_no_replace(output, values$worker_output)
  marker <- list(
    locus = output$LOC,
    chromosome = output$CHR,
    n_snps = output$n_snps,
    pair = output$pair,
    qc = output$state
  )
  rm(output, sealed_aggregate, candidates)
  invisible(gc(full = TRUE))
  worker_marker(
    "conditional", values$locus_index, marker$locus, marker$chromosome,
    marker$n_snps, marker$pair, marker$qc,
    if (marker$qc == "CONDITIONAL_LOCUS_PROCESSED") "TRUE" else "FALSE",
    values$worker_output
  )
  quit(save = "no", status = 0L)
}

conditional_relative_paths <- file.path(
  "conditional", sprintf("locus_%04d", seq_len(nrow(loci))), "result.rds"
)
conditional_paths <- file.path(continuation_root, conditional_relative_paths)
if (!all(file.exists(conditional_paths))) {
  stop(sprintf(
    "Incomplete continuation conditional checkpoint family: %d/%d",
    sum(file.exists(conditional_paths)), length(conditional_paths)
  ))
}
if (any(!vapply(conditional_paths, path_within, logical(1L), root = continuation_root))) {
  stop("A continuation conditional result resolves outside its checkpoint root")
}

validate_conditional_identity <- function(checkpoint, index) {
  expected_names <- c(
    "schema_version", "source_discovery_fingerprint",
    "continuation_execution_fingerprint", "phase", "locus_index", "LOC", "CHR",
    "aggregate_checkpoint_sha256", "n_snps", "pair", "state", "conditional"
  )
  official <- loci[index, , drop = FALSE]
  if (!is.list(checkpoint) || !identical(names(checkpoint), expected_names) ||
      !identical(checkpoint$schema_version, schema_version) ||
      !identical(checkpoint$source_discovery_fingerprint, values$source_fingerprint) ||
      !identical(checkpoint$continuation_execution_fingerprint, values$continuation_fingerprint) ||
      !identical(checkpoint$phase, "conditional") ||
      !identical(as.integer(checkpoint$locus_index), as.integer(index)) ||
      !identical(as.character(checkpoint$LOC), as.character(official$LOC)) ||
      !identical(as.character(checkpoint$CHR), as.character(official$CHR)) ||
      !identical(checkpoint$aggregate_checkpoint_sha256, aggregate_sha256)) {
    stop(paste("Stale or malformed continuation conditional checkpoint for locus index", index))
  }
  exact_frame_names(checkpoint$conditional, conditional_pre_fields, "Conditional result")
  checkpoint
}

conditional_reads <- lapply(seq_along(conditional_paths), function(index) {
  observed <- read_rds_with_sha256(
    conditional_paths[[index]], paste("Conditional checkpoint", index)
  )
  observed$value <- validate_conditional_identity(observed$value, index)
  observed
})
conditional_results <- lapply(conditional_reads, function(item) item$value)
conditional_manifest <- data.frame(
  locus_index = seq_len(nrow(loci)),
  LOC = as.character(loci$LOC),
  relative_path = conditional_relative_paths,
  sha256 = vapply(conditional_reads, function(item) item$sha256, character(1L)),
  stringsAsFactors = FALSE
)
rm(conditional_reads)
conditional_parts <- lapply(conditional_results, function(value) value$conditional)
conditional <- if (any(vapply(conditional_parts, nrow, integer(1L)) > 0L)) {
  do.call(rbind, conditional_parts)
} else {
  empty_conditional()
}
eligible_conditional <- which(
  conditional$analysis_status != "CONDITIONER_LOCAL_H2_INELIGIBLE"
)
conditional$conditional_test_family_n <- rep.int(
  length(eligible_conditional), nrow(conditional)
)
conditional$p_fdr <- rep(NA_real_, nrow(conditional))
tested_conditional <- which(
  conditional$analysis_status == "TESTED" & is.finite(conditional$p)
)
conditional_family_p <- rep.int(1, length(eligible_conditional))
names(conditional_family_p) <- eligible_conditional
conditional_family_p[as.character(tested_conditional)] <- conditional$p[tested_conditional]
conditional_adjusted <- p.adjust(conditional_family_p, method = "BH")
conditional$p_fdr[tested_conditional] <- conditional_adjusted[
  as.character(tested_conditional)
]
conditional$fdr_significant <- !is.na(conditional$p_fdr) &
  conditional$p_fdr <= policy_num("conditional_fdr_alpha")

attempt_dir <- dirname(values$worker_output)
write_table_no_replace(sealed_aggregate$status, file.path(attempt_dir, "lava_locus_status.tsv"))
write_table_no_replace(sealed_aggregate$univ, file.path(attempt_dir, "lava_univariate.tsv"))
write_table_no_replace(sealed_aggregate$bivar, file.path(attempt_dir, "lava_bivariate.tsv"))
write_table_no_replace(conditional, file.path(attempt_dir, "lava_conditional.tsv"))

staged_validator <- "scripts/138_validate_track_b_lava_results_v2.py"
if (!file.exists(staged_validator)) {
  stop(paste("Required continuation staged validator is missing:", staged_validator))
}
terminal_qc_path <- file.path(attempt_dir, "terminal_qc.json")
validation <- suppressWarnings(system2(
  "python3",
  c(
    staged_validator,
    "--validate-staged-results",
    "--staging-dir", shQuote(attempt_dir),
    "--source-fingerprint", values$source_fingerprint,
    "--continuation-fingerprint", values$continuation_fingerprint,
    "--terminal-qc-output", shQuote(terminal_qc_path)
  ),
  stdout = TRUE, stderr = TRUE
))
validation_status <- attr(validation, "status")
if (is.null(validation_status)) validation_status <- 0L
total_snps <- sum(sealed_aggregate$status$n_snps, na.rm = TRUE)
if (identical(as.integer(validation_status), 78L)) {
  if (!file.exists(terminal_qc_path) || file.info(terminal_qc_path)$size <= 0L) {
    stop("V2 staged validator exited 78 without a nonempty terminal_qc.json")
  }
  terminal_payload <- list(
    schema_version = schema_version,
    source_discovery_fingerprint = values$source_fingerprint,
    continuation_execution_fingerprint = values$continuation_fingerprint,
    phase = "terminal-qc",
    locus_index = 0L,
    LOC = "ALL",
    CHR = "ALL",
    aggregate_checkpoint_sha256 = aggregate_sha256,
    conditional_checkpoint_manifest = conditional_manifest,
    terminal_qc_sha256 = sha256_files(terminal_qc_path)[[1L]],
    total_snps = total_snps,
    state = "TERMINAL_FAILED_QC"
  )
  publish_rds_no_replace(terminal_payload, values$worker_output)
  worker_marker(
    "terminal-qc", 0L, "ALL", "ALL", total_snps, "A;B;CONTROL",
    "TERMINAL_FAILED_QC", "TRUE", values$worker_output
  )
  quit(save = "no", status = 78L)
}
if (validation_status != 0L) {
  stop(paste("Track B LAVA V2 staged validation failed:", paste(validation, collapse = " | ")))
}
if (sum(grepl(
  "^TRACK_B_LAVA_V2_STAGED_RESULTS_SEMANTICALLY_VALIDATED",
  validation
)) != 1L) {
  stop("V2 staged result validator did not emit one PASS marker")
}
staging_validation_path <- file.path(attempt_dir, "staging_validation.txt")
write_lines_no_replace(validation, staging_validation_path)

final_payload <- list(
  schema_version = schema_version,
  source_discovery_fingerprint = values$source_fingerprint,
  continuation_execution_fingerprint = values$continuation_fingerprint,
  phase = "finalize",
  locus_index = 0L,
  LOC = "ALL",
  CHR = "ALL",
  aggregate_checkpoint_sha256 = aggregate_sha256,
  conditional_checkpoint_manifest = conditional_manifest,
  staging_validation_sha256 = sha256_files(staging_validation_path)[[1L]],
  total_snps = total_snps,
  state = "STAGED_SEMANTIC_VALIDATION_COMPLETE"
)
publish_rds_no_replace(final_payload, values$worker_output)
rm(
  sealed_aggregate, conditional_results, conditional_parts, conditional,
  conditional_manifest, final_payload
)
invisible(gc(full = TRUE))
worker_marker(
  "finalize", 0L, "ALL", "ALL", total_snps, "A;B;CONTROL",
  "STAGED_SEMANTIC_VALIDATION_COMPLETE", "TRUE", values$worker_output
)
