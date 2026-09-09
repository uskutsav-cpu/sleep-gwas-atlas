#!/usr/bin/env Rscript

# Semantic validator for additive Track B post-discovery continuation
# checkpoints.  V1 discovery artifacts remain identified by their source
# fingerprint; V2 artifacts additionally carry a distinct execution identity.

schema_version <- "track-b-lava-continuation-checkpoint.1"
allowed_phases <- c("aggregate-discovery", "conditional", "finalize", "terminal-qc")

values <- list(
  phase = NA_character_,
  locus_index = NA_integer_,
  source_fingerprint = NA_character_,
  continuation_fingerprint = NA_character_,
  source_checkpoint_root = NA_character_,
  continuation_checkpoint_root = NA_character_,
  rds = NA_character_
)
args <- commandArgs(trailingOnly = TRUE)
mapping <- c(
  "--phase" = "phase",
  "--locus-index" = "locus_index",
  "--source-fingerprint" = "source_fingerprint",
  "--continuation-fingerprint" = "continuation_fingerprint",
  "--source-checkpoint-root" = "source_checkpoint_root",
  "--continuation-checkpoint-root" = "continuation_checkpoint_root",
  "--rds" = "rds"
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
  stop("Invalid V2 checkpoint phase")
}
if (!valid_fingerprint(values$source_fingerprint)) stop("Invalid source fingerprint")
if (!valid_fingerprint(values$continuation_fingerprint)) stop("Invalid continuation fingerprint")
if (identical(values$source_fingerprint, values$continuation_fingerprint)) {
  stop("Source and continuation fingerprints must differ")
}
if (values$phase == "conditional") {
  if (is.na(values$locus_index) || values$locus_index < 1L || values$locus_index > 2495L) {
    stop("Conditional checkpoint requires a locus index in 1..2495")
  }
} else if (!is.na(values$locus_index)) {
  stop("Aggregate/final/terminal checkpoint cannot have a locus index")
}
if (!scalar_string(values$rds) || !file.exists(values$rds) ||
    file.info(values$rds)$size <= 0L) {
  stop("Missing checkpoint RDS")
}

normalize_namespace_root <- function(path, fingerprint, label) {
  if (!scalar_string(path) || !dir.exists(path)) stop(paste(label, "does not exist"))
  normalized <- normalizePath(path, mustWork = TRUE)
  if (!identical(basename(normalized), fingerprint)) {
    stop(paste(label, "must be the exact fingerprint namespace root"))
  }
  if (grepl("[\r\n\t]", normalized)) stop(paste(label, "contains a forbidden control character"))
  normalized
}
source_root <- normalize_namespace_root(
  values$source_checkpoint_root, values$source_fingerprint, "Source checkpoint root"
)
continuation_root <- normalize_namespace_root(
  values$continuation_checkpoint_root, values$continuation_fingerprint,
  "Continuation checkpoint root"
)
if (identical(source_root, continuation_root)) stop("Source and continuation roots must differ")

path_within <- function(path, root) {
  resolved <- normalizePath(path, mustWork = TRUE)
  identical(resolved, root) || startsWith(paste0(resolved, "/"), paste0(root, "/"))
}
if (!path_within(values$rds, continuation_root)) {
  stop("Checkpoint RDS resolves outside the continuation namespace")
}

locus_lines <- strsplit(
  readLines("ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile", warn = FALSE),
  "[[:space:]]+"
)
locus_header <- locus_lines[[1L]]
if (!identical(locus_header, c("LOC", "CHR", "START", "STOP")) ||
    length(locus_lines) != 2496L) {
  stop("Official locus family drifted")
}
loci <- as.data.frame(do.call(rbind, locus_lines[-1L]), stringsAsFactors = FALSE)
names(loci) <- locus_header
if (anyDuplicated(loci$LOC)) stop("Official locus identifiers are not unique")

input_info <- read.delim(
  "results/track_b/lava_input_info.tsv", stringsAsFactors = FALSE, check.names = FALSE
)
pairs <- read.delim(
  "results/track_b/lava_pair_manifest.tsv", stringsAsFactors = FALSE, check.names = FALSE
)
conditioners <- read.delim(
  "results/track_b/local_conditional_manifest.tsv",
  stringsAsFactors = FALSE, check.names = FALSE
)
policy_table <- read.delim(
  "results/track_b/lava_runtime_policy.tsv", stringsAsFactors = FALSE, check.names = FALSE
)
policy <- setNames(as.list(policy_table$value), policy_table$key)
policy_num <- function(key) as.numeric(policy[[key]])
policy_int <- function(key) as.integer(policy[[key]])
if (policy$analysis_id != "track-b-v1.0-local" || policy$lava_version != "0.1.5") {
  stop("Track B runtime policy drifted")
}
if (nrow(input_info) != policy_int("expected_traits") ||
    anyDuplicated(input_info$phenotype)) {
  stop("Eight-trait input family drifted")
}
if (!identical(as.character(pairs$pair_id), c("A", "B", "CONTROL"))) {
  stop("Three-pair family drifted")
}
if (!identical(
  as.character(conditioners$conditional_model_id),
  c("A_BMI_ONLY", "A_SLEEP_APNEA_ONLY", "B_MDD_ONLY", "NONE")
)) {
  stop("Conditional model family drifted")
}

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
aggregate_names <- c(
  "schema_version", "source_discovery_fingerprint",
  "continuation_execution_fingerprint", "phase",
  "source_checkpoint_manifest", "status", "univ", "bivar"
)
conditional_names <- c(
  "schema_version", "source_discovery_fingerprint",
  "continuation_execution_fingerprint", "phase", "locus_index", "LOC", "CHR",
  "aggregate_checkpoint_sha256", "n_snps", "pair", "state", "conditional"
)
manifest_fields <- c("locus_index", "LOC", "relative_path", "sha256")

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
finite_or_na <- function(value) all(is.na(value) | is.finite(value))
has_explanation <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    nzchar(trimws(value))
}
marker_number <- function(value) {
  if (length(value) == 1L && is.finite(value) && value >= 0) {
    as.character(as.integer(value))
  } else {
    "NA"
  }
}
same_scalar_integer_or_na <- function(left, right) {
  (length(left) == 1L && length(right) == 1L) &&
    ((is.na(left) && is.na(right)) ||
      (!is.na(left) && !is.na(right) && as.integer(left) == as.integer(right)))
}

frames_equivalent <- function(observed, expected, tolerance = 1e-12) {
  if (!is.data.frame(observed) || !is.data.frame(expected) ||
      !identical(names(observed), names(expected)) ||
      nrow(observed) != nrow(expected)) {
    return(FALSE)
  }
  for (name in names(expected)) {
    left <- observed[[name]]
    right <- expected[[name]]
    if (is.factor(left)) left <- as.character(left)
    if (is.factor(right)) right <- as.character(right)
    if (is.numeric(left) || is.numeric(right)) {
      left <- as.numeric(left)
      right <- as.numeric(right)
      if (!identical(is.na(left), is.na(right))) return(FALSE)
      keep <- !is.na(left)
      if (any(keep) && !isTRUE(all.equal(
        left[keep], right[keep], tolerance = tolerance, check.attributes = FALSE
      ))) return(FALSE)
    } else if (!identical(as.character(left), as.character(right))) {
      return(FALSE)
    }
  }
  TRUE
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
conditional_relative_paths <- file.path(
  "conditional", sprintf("locus_%04d", seq_len(nrow(loci))), "result.rds"
)
conditional_paths <- file.path(continuation_root, conditional_relative_paths)
aggregate_path <- file.path(
  continuation_root, "aggregate-discovery", "unit_all", "result.rds"
)

validate_source_discovery <- function(checkpoint, index) {
  expected_names <- c("fingerprint", "phase", "locus_index", "LOC", "CHR", "status", "univ", "bivar")
  official <- loci[index, , drop = FALSE]
  if (!is.list(checkpoint) || !identical(names(checkpoint), expected_names) ||
      !identical(checkpoint$fingerprint, values$source_fingerprint) ||
      !identical(checkpoint$phase, "discovery") ||
      !identical(as.integer(checkpoint$locus_index), as.integer(index)) ||
      !identical(as.character(checkpoint$LOC), official$LOC) ||
      !identical(as.character(checkpoint$CHR), official$CHR)) {
    stop(paste("Source discovery identity drifted at locus index", index))
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
  status_value <- as.character(checkpoint$status$status[[1L]])
  if (!(status_value %in% c("PROCESSED", "PROCESS_FAILED", "UNIVARIATE_FAILED"))) {
    stop(paste("Invalid source discovery status at locus index", index))
  }
  allowed_univ <- c(
    "TESTED", "PHENOTYPE_DROPPED", "LOCUS_PROCESS_FAILED", "UNIVARIATE_FAILED"
  )
  if (any(!(checkpoint$univ$analysis_status %in% allowed_univ))) {
    stop(paste("Invalid source univariate status at locus index", index))
  }
  if (any(!(checkpoint$bivar$pair_id %in% pairs$pair_id)) ||
      anyDuplicated(checkpoint$bivar$pair_id) ||
      !same_coordinates(checkpoint$bivar, official)) {
    stop(paste("Source discovery bivariate family drifted at locus index", index))
  }
  pair_index <- match(checkpoint$bivar$pair_id, pairs$pair_id)
  if (nrow(checkpoint$bivar) &&
      any(checkpoint$bivar$trait1 != pairs$trait1[pair_index] |
          checkpoint$bivar$trait2 != pairs$trait2[pair_index])) {
    stop(paste("Source discovery bivariate traits drifted at locus index", index))
  }
  if (any(!(checkpoint$bivar$analysis_status %in% c("TESTED", "BIVARIATE_FAILED"))) ||
      !finite_or_na(checkpoint$univ$p) || !finite_or_na(checkpoint$bivar$p)) {
    stop(paste("Source discovery result status/numerics drifted at locus index", index))
  }
  tested_univ <- checkpoint$univ$analysis_status == "TESTED"
  if (any(tested_univ & (
    !is.finite(checkpoint$univ$p) |
      checkpoint$univ$p < 0 | checkpoint$univ$p > 1 |
      !is.finite(checkpoint$univ$h2.obs) | checkpoint$univ$h2.obs < 0
  ))) {
    stop(paste("Invalid tested source univariate row at locus index", index))
  }
  binary_traits <- as.character(input_info$phenotype[!is.na(input_info$cases)])
  tested_binary <- tested_univ & checkpoint$univ$phen %in% binary_traits
  tested_continuous <- tested_univ & !(checkpoint$univ$phen %in% binary_traits)
  if (any(tested_binary & (
    !is.finite(checkpoint$univ$h2.latent) | checkpoint$univ$h2.latent < 0
  )) || any(tested_continuous & !is.na(checkpoint$univ$h2.latent))) {
    stop(paste("Source latent h2 semantics drifted at locus index", index))
  }
  if (status_value == "PROCESS_FAILED" &&
      (!all(checkpoint$univ$analysis_status == "LOCUS_PROCESS_FAILED") ||
        nrow(checkpoint$bivar) != 0L)) {
    stop(paste("PROCESS_FAILED source rows drifted at locus index", index))
  }
  if (status_value == "UNIVARIATE_FAILED" &&
      (!all(checkpoint$univ$analysis_status == "UNIVARIATE_FAILED") ||
        nrow(checkpoint$bivar) != 0L)) {
    stop(paste("UNIVARIATE_FAILED source rows drifted at locus index", index))
  }
  if (status_value == "PROCESSED" &&
      any(checkpoint$univ$analysis_status %in%
        c("LOCUS_PROCESS_FAILED", "UNIVARIATE_FAILED"))) {
    stop(paste("PROCESSED source rows drifted at locus index", index))
  }
  eligible_traits <- checkpoint$univ$phen[
    tested_univ & is.finite(checkpoint$univ$p) &
      checkpoint$univ$p <= policy_num("univariate_p_threshold")
  ]
  expected_pairs <- pairs$pair_id[
    pairs$trait1 %in% eligible_traits & pairs$trait2 %in% eligible_traits
  ]
  if (!identical(as.character(checkpoint$bivar$pair_id), as.character(expected_pairs))) {
    stop(paste("Source bivariate rows differ from local-h2 gate at locus index", index))
  }
  tested_bivar <- checkpoint$bivar$analysis_status == "TESTED"
  numeric_bivar <- c(
    "local_covariance", "rho", "rho.lower", "rho.upper",
    "r2", "r2.lower", "r2.upper", "p"
  )
  if (nrow(checkpoint$bivar) && any(vapply(
    checkpoint$bivar[numeric_bivar],
    function(column) any(tested_bivar & !is.finite(column)),
    logical(1L)
  ))) {
    stop(paste("Invalid tested source bivariate row at locus index", index))
  }
  if (any(tested_bivar & (checkpoint$bivar$p < 0 | checkpoint$bivar$p > 1)) ||
      any(!tested_bivar & (!is.na(checkpoint$bivar$p) | is.na(checkpoint$bivar$error)))) {
    stop(paste("Source bivariate failure semantics drifted at locus index", index))
  }
  checkpoint
}

validate_manifest_shape <- function(manifest, relative_paths, label) {
  if (!is.data.frame(manifest) || !identical(names(manifest), manifest_fields) ||
      nrow(manifest) != nrow(loci) ||
      !identical(as.integer(manifest$locus_index), seq_len(nrow(loci))) ||
      !identical(as.character(manifest$LOC), loci$LOC) ||
      !identical(as.character(manifest$relative_path), relative_paths) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$sha256)) ||
      anyDuplicated(manifest$sha256)) {
    stop(paste(label, "manifest drifted"))
  }
  invisible(TRUE)
}

validate_aggregate_identity <- function(checkpoint) {
  if (!is.list(checkpoint) || !identical(names(checkpoint), aggregate_names) ||
      !identical(checkpoint$schema_version, schema_version) ||
      !identical(checkpoint$source_discovery_fingerprint, values$source_fingerprint) ||
      !identical(
        checkpoint$continuation_execution_fingerprint,
        values$continuation_fingerprint
      ) ||
      !identical(checkpoint$phase, "aggregate-discovery")) {
    stop("Continuation aggregate identity drifted")
  }
  validate_manifest_shape(
    checkpoint$source_checkpoint_manifest, source_relative_paths, "Source checkpoint"
  )
  exact_frame_names(
    checkpoint$status,
    c(status_pre_fields, "analysis_fingerprint", "lava_version", "reference_prefix"),
    "Aggregated status"
  )
  exact_frame_names(
    checkpoint$univ,
    c(
      univ_pre_fields, "univariate_test_family_n", "univariate_p_threshold",
      "p_bonferroni", "p_fdr"
    ),
    "Aggregated univariate"
  )
  exact_frame_names(
    checkpoint$bivar,
    c(bivar_pre_fields, "bivariate_test_family_n", "p_fdr", "fdr_significant"),
    "Aggregated bivariate"
  )
  if (nrow(checkpoint$status) != nrow(loci) ||
      !identical(as.character(checkpoint$status$LOC), loci$LOC) ||
      any(checkpoint$status$analysis_fingerprint != values$source_fingerprint) ||
      any(checkpoint$status$lava_version != policy$lava_version) ||
      any(checkpoint$status$reference_prefix != policy$reference_prefix)) {
    stop("Aggregated status lineage drifted")
  }
  expected_univ_loci <- rep(loci$LOC, each = policy_int("expected_traits"))
  expected_univ_traits <- rep(
    as.character(input_info$phenotype), times = policy_int("expected_loci")
  )
  if (nrow(checkpoint$univ) != policy_int("planned_univariate_tests") ||
      !identical(as.character(checkpoint$univ$LOC), expected_univ_loci) ||
      !identical(as.character(checkpoint$univ$phen), expected_univ_traits) ||
      any(checkpoint$univ$univariate_test_family_n !=
        policy_int("planned_univariate_tests"))) {
    stop("Aggregated univariate family drifted")
  }
  if (any(!(checkpoint$bivar$pair_id %in% pairs$pair_id)) ||
      anyDuplicated(paste(checkpoint$bivar$LOC, checkpoint$bivar$pair_id)) ||
      (nrow(checkpoint$bivar) &&
        any(checkpoint$bivar$bivariate_test_family_n != nrow(checkpoint$bivar)))) {
    stop("Aggregated bivariate family drifted")
  }
  checkpoint
}

build_expected_aggregate_tables <- function(source_results) {
  status <- do.call(rbind, lapply(source_results, function(value) value$status))
  univ <- do.call(rbind, lapply(source_results, function(value) value$univ))
  bivar <- do.call(rbind, lapply(source_results, function(value) value$bivar))
  if (nrow(status) != policy_int("expected_loci") ||
      nrow(univ) != policy_int("planned_univariate_tests")) {
    stop("Source discovery contents are incomplete")
  }
  planned_univ <- policy_int("planned_univariate_tests")
  univ$univariate_test_family_n <- planned_univ
  univ$univariate_p_threshold <- policy_num("univariate_p_threshold")
  univ$p_bonferroni <- pmin(1, univ$p * planned_univ)
  univ$p_fdr <- NA_real_
  tested_univ <- which(univ$analysis_status == "TESTED" & is.finite(univ$p))
  univ_family_p <- rep.int(1, nrow(univ))
  univ_family_p[tested_univ] <- univ$p[tested_univ]
  expected_univ_fdr <- p.adjust(univ_family_p, method = "BH")
  univ$p_fdr[tested_univ] <- expected_univ_fdr[tested_univ]

  bivar$bivariate_test_family_n <- rep.int(nrow(bivar), nrow(bivar))
  bivar$p_fdr <- rep(NA_real_, nrow(bivar))
  tested_bivar <- which(bivar$analysis_status == "TESTED" & is.finite(bivar$p))
  bivar_family_p <- rep.int(1, nrow(bivar))
  bivar_family_p[tested_bivar] <- bivar$p[tested_bivar]
  expected_bivar_fdr <- p.adjust(bivar_family_p, method = "BH")
  bivar$p_fdr[tested_bivar] <- expected_bivar_fdr[tested_bivar]
  bivar$fdr_significant <- !is.na(bivar$p_fdr) &
    bivar$p_fdr <= policy_num("bivariate_fdr_alpha")

  status$analysis_fingerprint <- values$source_fingerprint
  status$lava_version <- policy$lava_version
  status$reference_prefix <- policy$reference_prefix
  list(status = status, univ = univ, bivar = bivar)
}

build_expected_candidates <- function(family) {
  rows <- list()
  output_index <- 0L
  for (index in seq_len(nrow(loci))) {
    locus_id <- loci$LOC[[index]]
    targets <- family$bivar[
      as.character(family$bivar$LOC) == locus_id &
        family$bivar$fdr_significant &
        family$bivar$pair_id %in% c("A", "B"),
      , drop = FALSE
    ]
    eligible_models <- 0L
    eligible_pairs <- character(0)
    if (nrow(targets)) for (target_index in seq_len(nrow(targets))) {
      models <- conditioners[
        conditioners$pair_id == targets$pair_id[[target_index]] &
          conditioners$conditional_model_id != "NONE",
        , drop = FALSE
      ]
      if (nrow(models)) for (model_index in seq_len(nrow(models))) {
        covars <- strsplit(models$covariates[[model_index]], ";", fixed = TRUE)[[1L]]
        local_univ <- family$univ[
          as.character(family$univ$LOC) == locus_id &
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
          eligible_pairs <- unique(c(
            eligible_pairs, as.character(targets$pair_id[[target_index]])
          ))
        }
      }
    }
    if (eligible_models > 0L) {
      output_index <- output_index + 1L
      rows[[output_index]] <- data.frame(
        locus_index = index,
        locus = locus_id,
        chromosome = loci$CHR[[index]],
        start = loci$START[[index]],
        stop = loci$STOP[[index]],
        n_snps = family$status$n_snps[
          match(locus_id, as.character(family$status$LOC))
        ],
        eligible_models = eligible_models,
        eligible_pairs = paste(eligible_pairs, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame(
      locus_index = integer(0), locus = integer(0), chromosome = integer(0),
      start = integer(0), stop = integer(0), n_snps = integer(0),
      eligible_models = integer(0), eligible_pairs = character(0),
      stringsAsFactors = FALSE
    )
  }
}

validate_full_aggregate <- function(checkpoint, candidate_path) {
  if (!all(file.exists(source_paths))) stop("Complete source family is unavailable")
  if (any(!vapply(source_paths, path_within, logical(1L), root = source_root))) {
    stop("A source discovery result resolves outside its checkpoint root")
  }
  source_reads <- lapply(seq_along(source_paths), function(index) {
    observed <- read_rds_with_sha256(
      source_paths[[index]], paste("Source discovery checkpoint", index)
    )
    observed$value <- validate_source_discovery(observed$value, index)
    observed
  })
  observed_manifest <- data.frame(
    locus_index = seq_len(nrow(loci)),
    LOC = loci$LOC,
    relative_path = source_relative_paths,
    sha256 = vapply(source_reads, function(item) item$sha256, character(1L)),
    stringsAsFactors = FALSE
  )
  if (!identical(
    as.character(checkpoint$source_checkpoint_manifest$sha256),
    as.character(observed_manifest$sha256)
  )) {
    stop("Aggregate source-result hashes drifted")
  }
  source_results <- lapply(source_reads, function(item) item$value)
  rm(source_reads)
  expected <- build_expected_aggregate_tables(source_results)
  if (!frames_equivalent(checkpoint$status, expected$status) ||
      !frames_equivalent(checkpoint$univ, expected$univ) ||
      !frames_equivalent(checkpoint$bivar, expected$bivar)) {
    stop("Aggregate tables differ from the immutable V1 discovery family")
  }
  tested_univ <- checkpoint$univ$analysis_status == "TESTED" &
    is.finite(checkpoint$univ$p)
  univ_family_p <- rep.int(1, nrow(checkpoint$univ))
  univ_family_p[tested_univ] <- checkpoint$univ$p[tested_univ]
  recomputed_univ <- p.adjust(univ_family_p, method = "BH")
  if (!identical(
    as.numeric(checkpoint$univ$p_fdr[tested_univ]),
    as.numeric(recomputed_univ[tested_univ])
  ) || any(!tested_univ & !is.na(checkpoint$univ$p_fdr))) {
    stop("Full-family univariate BH values drifted")
  }
  tested_bivar <- checkpoint$bivar$analysis_status == "TESTED" &
    is.finite(checkpoint$bivar$p)
  bivar_family_p <- rep.int(1, nrow(checkpoint$bivar))
  bivar_family_p[tested_bivar] <- checkpoint$bivar$p[tested_bivar]
  recomputed_bivar <- p.adjust(bivar_family_p, method = "BH")
  if (!identical(
    as.numeric(checkpoint$bivar$p_fdr[tested_bivar]),
    as.numeric(recomputed_bivar[tested_bivar])
  ) || any(!tested_bivar & !is.na(checkpoint$bivar$p_fdr)) ||
      !identical(
        as.logical(checkpoint$bivar$fdr_significant),
        !is.na(checkpoint$bivar$p_fdr) &
          checkpoint$bivar$p_fdr <= policy_num("bivariate_fdr_alpha")
      )) {
    stop("Full-family bivariate BH values drifted")
  }

  expected_pair_keys <- character(0)
  for (locus_id in loci$LOC) {
    local_univ <- checkpoint$univ[
      as.character(checkpoint$univ$LOC) == locus_id,
      , drop = FALSE
    ]
    eligible_traits <- local_univ$phen[
      local_univ$analysis_status == "TESTED" &
        is.finite(local_univ$p) &
        local_univ$p <= policy_num("univariate_p_threshold")
    ]
    local_pairs <- pairs[
      pairs$trait1 %in% eligible_traits & pairs$trait2 %in% eligible_traits,
      , drop = FALSE
    ]
    if (nrow(local_pairs)) {
      expected_pair_keys <- c(
        expected_pair_keys,
        paste(locus_id, as.character(local_pairs$pair_id), sep = "\t")
      )
    }
  }
  observed_pair_keys <- paste(
    as.character(checkpoint$bivar$LOC), checkpoint$bivar$pair_id, sep = "\t"
  )
  if (!identical(observed_pair_keys, expected_pair_keys)) {
    stop("Aggregate bivariate family omits or adds locally eligible tests")
  }

  if (!file.exists(candidate_path) || file.info(candidate_path)$size <= 0L) {
    stop("Aggregate checkpoint lacks conditional candidate attestation")
  }
  candidates <- read.delim(
    candidate_path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA"
  )
  expected_candidates <- build_expected_candidates(checkpoint)
  if (!frames_equivalent(candidates, expected_candidates)) {
    stop("Conditional candidates differ from the frozen full-family gates")
  }
  checkpoint
}

load_aggregate <- function(full_source_validation = FALSE) {
  if (!file.exists(aggregate_path) || file.info(aggregate_path)$size <= 0L ||
      !path_within(aggregate_path, continuation_root)) {
    stop("Canonical continuation aggregate checkpoint is missing or escapes its namespace")
  }
  observed <- read_rds_with_sha256(aggregate_path, "Continuation aggregate checkpoint")
  checkpoint <- validate_aggregate_identity(observed$value)
  if (full_source_validation) {
    checkpoint <- validate_full_aggregate(
      checkpoint, file.path(dirname(aggregate_path), "conditional_candidates.tsv")
    )
  }
  list(value = checkpoint, sha256 = observed$sha256)
}

expected_conditional_rows <- function(index, family) {
  official <- loci[index, , drop = FALSE]
  targets <- family$bivar[
    as.character(family$bivar$LOC) == official$LOC &
      family$bivar$fdr_significant &
      family$bivar$pair_id %in% c("A", "B"),
    , drop = FALSE
  ]
  rows <- list()
  row_index <- 0L
  eligibility <- logical(0)
  if (nrow(targets)) for (target_index in seq_len(nrow(targets))) {
    target <- targets[target_index, , drop = FALSE]
    models <- conditioners[
      conditioners$pair_id == target$pair_id &
        conditioners$conditional_model_id != "NONE",
      , drop = FALSE
    ]
    if (!nrow(models)) stop(paste("No frozen conditional model for pair", target$pair_id))
    for (model_index in seq_len(nrow(models))) {
      row_index <- row_index + 1L
      model <- models[model_index, , drop = FALSE]
      covars <- strsplit(model$covariates, ";", fixed = TRUE)[[1L]]
      local_univ <- family$univ[
        as.character(family$univ$LOC) == official$LOC &
          family$univ$phen %in% covars,
        , drop = FALSE
      ]
      eligibility[[row_index]] <- length(covars) > 0L &&
        nrow(local_univ) == length(covars) &&
        all(
          local_univ$analysis_status == "TESTED" &
            is.finite(local_univ$p) &
            local_univ$p <= policy_num("univariate_p_threshold")
        )
      rows[[row_index]] <- data.frame(
        LOC = target$LOC, CHR = target$CHR, START = target$START, STOP = target$STOP,
        pair_id = target$pair_id, trait1 = target$trait1, trait2 = target$trait2,
        conditional_model_id = model$conditional_model_id,
        covariates = model$covariates,
        stringsAsFactors = FALSE
      )
    }
  }
  base <- if (length(rows)) do.call(rbind, rows) else data.frame(
    LOC = numeric(0), CHR = numeric(0), START = numeric(0), STOP = numeric(0),
    pair_id = character(0), trait1 = character(0), trait2 = character(0),
    conditional_model_id = character(0), covariates = character(0),
    stringsAsFactors = FALSE
  )
  list(targets = targets, base = base, eligible = eligibility)
}

validate_conditional <- function(
  checkpoint, index, aggregate, aggregate_sha256, check_local_source = TRUE
) {
  official <- loci[index, , drop = FALSE]
  if (!is.list(checkpoint) || !identical(names(checkpoint), conditional_names) ||
      !identical(checkpoint$schema_version, schema_version) ||
      !identical(checkpoint$source_discovery_fingerprint, values$source_fingerprint) ||
      !identical(
        checkpoint$continuation_execution_fingerprint,
        values$continuation_fingerprint
      ) ||
      !identical(checkpoint$phase, "conditional") ||
      !identical(as.integer(checkpoint$locus_index), as.integer(index)) ||
      !identical(as.character(checkpoint$LOC), official$LOC) ||
      !identical(as.character(checkpoint$CHR), official$CHR) ||
      !identical(checkpoint$aggregate_checkpoint_sha256, aggregate_sha256)) {
    stop(paste("Continuation conditional identity drifted at locus index", index))
  }
  exact_frame_names(checkpoint$conditional, conditional_pre_fields, "Conditional result")
  if (!same_coordinates(checkpoint$conditional, official) ||
      any(!(checkpoint$conditional$pair_id %in% c("A", "B"))) ||
      anyDuplicated(paste(
        checkpoint$conditional$pair_id,
        checkpoint$conditional$conditional_model_id
      ))) {
    stop(paste("Conditional row identity drifted at locus index", index))
  }
  expected <- expected_conditional_rows(index, aggregate)
  if (nrow(checkpoint$conditional) != nrow(expected$base) ||
      !frames_equivalent(
        checkpoint$conditional[names(expected$base)], expected$base
      )) {
    stop(paste("Conditional checkpoint omits or adds a frozen model at locus index", index))
  }
  expected_pair <- if (nrow(expected$targets)) {
    paste(expected$targets$pair_id, collapse = ";")
  } else {
    "NONE"
  }
  expected_n_snps <- aggregate$status$n_snps[
    match(official$LOC, as.character(aggregate$status$LOC))
  ]
  if (!identical(as.character(checkpoint$pair), expected_pair) ||
      !same_scalar_integer_or_na(checkpoint$n_snps, expected_n_snps)) {
    stop(paste("Conditional summary identity drifted at locus index", index))
  }
  allowed_status <- c(
    "TESTED", "CONDITIONER_LOCAL_H2_INELIGIBLE",
    "CONDITIONAL_FAILED", "CONDITIONAL_UNSTABLE_MAX_R2"
  )
  numeric_fields <- c(
    "pcor", "ci.lower", "ci.upper", "p", "r2.trait1_z", "r2.trait2_z"
  )
  if (!is.character(checkpoint$conditional$analysis_status) ||
      !is.character(checkpoint$conditional$error) ||
      any(!vapply(checkpoint$conditional[numeric_fields], is.numeric, logical(1L))) ||
      any(!(checkpoint$conditional$analysis_status %in% allowed_status))) {
    stop(paste("Invalid conditional terminal status at locus index", index))
  }
  if (nrow(checkpoint$conditional)) for (row_index in seq_len(nrow(checkpoint$conditional))) {
    status_value <- checkpoint$conditional$analysis_status[[row_index]]
    error_value <- checkpoint$conditional$error[[row_index]]
    diagnostics <- vapply(
      numeric_fields,
      function(field) checkpoint$conditional[[field]][[row_index]],
      numeric(1L)
    )
    eligible <- expected$eligible[[row_index]]
    if (!eligible) {
      if (!identical(status_value, "CONDITIONER_LOCAL_H2_INELIGIBLE") ||
          !has_explanation(error_value) || any(!is.na(diagnostics))) {
        stop(paste("Conditioner-ineligible row semantics drifted at locus index", index))
      }
      next
    }
    if (identical(status_value, "TESTED")) {
      if (!is.na(error_value) || any(!is.finite(diagnostics)) ||
          diagnostics[["pcor"]] < -1 || diagnostics[["pcor"]] > 1 ||
          diagnostics[["ci.lower"]] > diagnostics[["pcor"]] ||
          diagnostics[["ci.upper"]] < diagnostics[["pcor"]] ||
          diagnostics[["p"]] < 0 || diagnostics[["p"]] > 1 ||
          diagnostics[["r2.trait1_z"]] < 0 ||
          diagnostics[["r2.trait1_z"]] >= policy_num("conditional_max_r2") ||
          diagnostics[["r2.trait2_z"]] < 0 ||
          diagnostics[["r2.trait2_z"]] >= policy_num("conditional_max_r2")) {
        stop(paste("Tested conditional row semantics drifted at locus index", index))
      }
      next
    }
    if (!(status_value %in% c("CONDITIONAL_FAILED", "CONDITIONAL_UNSTABLE_MAX_R2")) ||
        !has_explanation(error_value) ||
        any(!is.na(diagnostics) & !is.finite(diagnostics)) ||
        (!is.na(diagnostics[["pcor"]]) &&
          (diagnostics[["pcor"]] < -1 || diagnostics[["pcor"]] > 1)) ||
        (!is.na(diagnostics[["p"]]) &&
          (diagnostics[["p"]] < 0 || diagnostics[["p"]] > 1)) ||
        any(!is.na(diagnostics[c("r2.trait1_z", "r2.trait2_z")]) &
          (diagnostics[c("r2.trait1_z", "r2.trait2_z")] < 0 |
            diagnostics[c("r2.trait1_z", "r2.trait2_z")] > 1))) {
      stop(paste("Failed conditional row semantics drifted at locus index", index))
    }
    if (identical(status_value, "CONDITIONAL_UNSTABLE_MAX_R2")) {
      retained <- diagnostics[c(
        "pcor", "ci.lower", "ci.upper", "r2.trait1_z", "r2.trait2_z"
      )]
      if (any(is.na(retained)) || !is.na(diagnostics[["p"]]) ||
          max(diagnostics[c("r2.trait1_z", "r2.trait2_z")]) <
            policy_num("conditional_max_r2")) {
        stop(paste("Conditional max-r2 failure semantics drifted at locus index", index))
      }
    }
  }
  valid_states <- c(
    "NOT_BIVARIATE_FDR_ELIGIBLE", "NO_FROZEN_MODEL",
    "CONDITIONER_LOCAL_H2_INELIGIBLE", "CONDITIONAL_LOCUS_PROCESSED",
    "CONDITIONAL_LOCUS_PROCESS_FAILED"
  )
  if (!(checkpoint$state %in% valid_states)) {
    stop(paste("Invalid conditional process state at locus index", index))
  }
  if (!nrow(expected$targets) &&
      !identical(checkpoint$state, "NOT_BIVARIATE_FDR_ELIGIBLE")) {
    stop(paste("Ineligible conditional locus state drifted at locus index", index))
  }
  if (nrow(expected$targets) && length(expected$eligible) &&
      !any(expected$eligible) &&
      !identical(checkpoint$state, "CONDITIONER_LOCAL_H2_INELIGIBLE")) {
    stop(paste("Conditioner-ineligible locus state drifted at locus index", index))
  }
  if (any(expected$eligible) &&
      !(checkpoint$state %in%
        c("CONDITIONAL_LOCUS_PROCESSED", "CONDITIONAL_LOCUS_PROCESS_FAILED"))) {
    stop(paste("Conditionally processed locus state drifted at locus index", index))
  }
  if (identical(checkpoint$state, "CONDITIONAL_LOCUS_PROCESS_FAILED") &&
      any(
        expected$eligible &
          checkpoint$conditional$analysis_status != "CONDITIONAL_FAILED"
      )) {
    stop(paste("Conditional locus-process failure rows drifted at locus index", index))
  }

  if (check_local_source) {
    source_path <- source_paths[[index]]
    if (!file.exists(source_path) || !path_within(source_path, source_root)) {
      stop(paste("Local source checkpoint is missing at locus index", index))
    }
    source_read <- read_rds_with_sha256(
      source_path, paste("Local source checkpoint", index)
    )
    source_hash <- source_read$sha256
    if (!identical(
      source_hash,
      as.character(aggregate$source_checkpoint_manifest$sha256[[index]])
    )) {
      stop(paste("Local source lineage hash drifted at locus index", index))
    }
    validate_source_discovery(source_read$value, index)
  }
  checkpoint
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

build_final_conditional <- function(checkpoints) {
  parts <- lapply(checkpoints, function(value) value$conditional)
  conditional <- if (any(vapply(parts, nrow, integer(1L)) > 0L)) {
    do.call(rbind, parts)
  } else {
    empty_conditional()
  }
  eligible <- which(
    conditional$analysis_status != "CONDITIONER_LOCAL_H2_INELIGIBLE"
  )
  conditional$conditional_test_family_n <- rep.int(length(eligible), nrow(conditional))
  conditional$p_fdr <- rep(NA_real_, nrow(conditional))
  tested <- which(conditional$analysis_status == "TESTED" & is.finite(conditional$p))
  family_p <- rep.int(1, length(eligible))
  names(family_p) <- eligible
  family_p[as.character(tested)] <- conditional$p[tested]
  adjusted <- p.adjust(family_p, method = "BH")
  conditional$p_fdr[tested] <- adjusted[as.character(tested)]
  conditional$fdr_significant <- !is.na(conditional$p_fdr) &
    conditional$p_fdr <= policy_num("conditional_fdr_alpha")
  conditional
}

validate_staged_tables <- function(directory, aggregate, conditional, require_derived) {
  required_names <- c(
    "lava_locus_status.tsv", "lava_univariate.tsv",
    "lava_bivariate.tsv", "lava_conditional.tsv"
  )
  if (require_derived) {
    required_names <- c(
      required_names, "04_lava_local_results.tsv", "05_local_conditional_results.tsv"
    )
  }
  required_paths <- file.path(directory, required_names)
  if (!all(file.exists(required_paths)) || any(file.info(required_paths)$size <= 0L)) {
    stop("Terminal checkpoint lacks its staged result family")
  }
  observed_status <- read.delim(
    file.path(directory, "lava_locus_status.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA"
  )
  observed_univ <- read.delim(
    file.path(directory, "lava_univariate.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA"
  )
  observed_bivar <- read.delim(
    file.path(directory, "lava_bivariate.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA"
  )
  observed_conditional <- read.delim(
    file.path(directory, "lava_conditional.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA"
  )
  if (!frames_equivalent(observed_status, aggregate$status) ||
      !frames_equivalent(observed_univ, aggregate$univ) ||
      !frames_equivalent(observed_bivar, aggregate$bivar) ||
      !frames_equivalent(observed_conditional, conditional)) {
    stop("Staged result tables differ from their checkpoint dependencies")
  }
  invisible(TRUE)
}

validate_terminal_qc_json <- function(path, aggregate, conditional) {
  counts <- c(
    sum(aggregate$status$status != "PROCESSED"), nrow(aggregate$status),
    sum(aggregate$univ$analysis_status != "TESTED"), nrow(aggregate$univ),
    sum(aggregate$bivar$analysis_status == "BIVARIATE_FAILED"), nrow(aggregate$bivar),
    sum(conditional$analysis_status %in%
      c("CONDITIONAL_FAILED", "CONDITIONAL_UNSTABLE_MAX_R2")),
    sum(conditional$analysis_status != "CONDITIONER_LOCAL_H2_INELIGIBLE"),
    sum(conditional$analysis_status == "CONDITIONER_LOCAL_H2_INELIGIBLE"),
    nrow(conditional)
  )
  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != as.integer(counts))) {
    stop("Cannot derive exact terminal-QC counts from checkpoint dependencies")
  }
  thresholds <- c(
    policy_num("maximum_locus_failure_fraction"),
    policy_num("maximum_univariate_untested_fraction"),
    policy_num("maximum_bivariate_failure_fraction"),
    policy_num("maximum_conditional_failure_fraction")
  )
  if (any(!is.finite(thresholds)) || any(thresholds < 0 | thresholds > 1)) {
    stop("Cannot derive exact terminal-QC thresholds from runtime policy")
  }
  python <- paste(c(
    "import hashlib, json, math, re, sys",
    "from pathlib import Path",
    "def fail(message): raise SystemExit(message)",
    "path = Path(sys.argv[1])",
    "source, continuation = sys.argv[2], sys.argv[3]",
    "numbers = [int(value) for value in sys.argv[4:14]]",
    "limits = [float(value) for value in sys.argv[14:18]]",
    "root = Path(sys.argv[18]).resolve(strict=True)",
    "try: payload = json.loads(path.read_text(encoding='utf-8'))",
    "except (OSError, json.JSONDecodeError) as error: fail(f'unreadable terminal JSON: {error}')",
    "top = {'schema_version','state','analysis_id','source_discovery_fingerprint','continuation_execution_fingerprint','full_family_complete','scientific_validation_passed','canonical_publication_allowed','terminal','validation_basis','qc','staged_low_level_results','upstream_continuation_family','source_family_lock','source_checkpoint_family_sha256','lineage','validators'}",
    "if not isinstance(payload, dict) or set(payload) != top: fail('terminal JSON schema drifted')",
    "if payload['schema_version'] != 'track-b-lava-terminal-qc.1' or payload['state'] != 'TERMINAL_FAILED_QC' or payload['analysis_id'] != 'track-b-v1.0-local': fail('terminal JSON identity drifted')",
    "if payload['source_discovery_fingerprint'] != source or payload['continuation_execution_fingerprint'] != continuation: fail('terminal JSON fingerprint lineage drifted')",
    "if payload['full_family_complete'] is not True or payload['scientific_validation_passed'] is not False or payload['canonical_publication_allowed'] is not False or payload['terminal'] is not True: fail('terminal JSON disposition drifted')",
    "basis = 'UNMODIFIED_LEGACY_SCHEMA_TOPOLOGY_NUMERICAL_AND_FULL_FAMILY_BH_CHECKS_WITH_QC_THRESHOLDS_DIAGNOSTICALLY_DEFERRED_THEN_EXACT_FROZEN_GATE_EVALUATION'",
    "if payload['validation_basis'] != basis: fail('terminal JSON validation basis drifted')",
    "locus_failed,locus_n,univ_failed,univ_n,bivar_failed,bivar_n,cond_failed,cond_eligible,cond_ineligible,cond_n = numbers",
    "expected_counts = {'locus': {'failed': locus_failed, 'family': locus_n}, 'univariate': {'failed_or_untested': univ_failed, 'family': univ_n}, 'bivariate': {'failed': bivar_failed, 'eligible_family': bivar_n}, 'conditional': {'failed': cond_failed, 'eligible_family': cond_eligible, 'ineligible': cond_ineligible, 'reported_family': cond_n}}",
    "expected_thresholds = dict(zip(('maximum_locus_failure_fraction','maximum_univariate_untested_fraction','maximum_bivariate_failure_fraction','maximum_conditional_failure_fraction'), limits))",
    "qc = payload['qc']",
    "count_values = [value for family in expected_counts for value in qc.get('counts', {}).get(family, {}).values()] if isinstance(qc, dict) and isinstance(qc.get('counts'), dict) else []",
    "threshold_values = list(qc.get('thresholds', {}).values()) if isinstance(qc, dict) and isinstance(qc.get('thresholds'), dict) else []",
    "if any(type(value) is not int for value in count_values) or any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) for value in threshold_values): fail('terminal JSON QC numeric types drifted')",
    "if not isinstance(qc, dict) or set(qc) != {'thresholds','counts','reasons'} or qc['thresholds'] != expected_thresholds or qc['counts'] != expected_counts: fail('terminal JSON QC counts or thresholds drifted')",
    "specs = [('LOCUS_FAILURE_FRACTION_EXCEEDED','locus_failure_fraction',locus_failed,locus_n,limits[0]),('UNIVARIATE_UNTESTED_FRACTION_EXCEEDED','univariate_untested_fraction',univ_failed,univ_n,limits[1]),('BIVARIATE_FAILURE_FRACTION_EXCEEDED','bivariate_failure_fraction',bivar_failed,bivar_n,limits[2]),('CONDITIONAL_FAILURE_FRACTION_EXCEEDED','conditional_failure_fraction',cond_failed,cond_eligible,limits[3])]",
    "def reason(code, metric, numerator, denominator, threshold):",
    "    observed = numerator / denominator if denominator else 0.0",
    "    return {'code': code, 'metric': metric, 'numerator': numerator, 'denominator': denominator, 'observed_fraction': observed, 'observed_fraction_decimal': format(observed, '.17g'), 'comparison': '>', 'maximum_allowed_fraction': threshold, 'maximum_allowed_fraction_decimal': format(threshold, '.17g')}",
    "expected_reasons = [reason(*spec) for spec in specs if spec[3] and spec[2] / spec[3] > spec[4]]",
    "if not expected_reasons or qc['reasons'] != expected_reasons: fail('terminal JSON QC reasons drifted')",
    "if any(type(item['numerator']) is not int or type(item['denominator']) is not int or isinstance(item['observed_fraction'], bool) or not isinstance(item['observed_fraction'], (int, float)) or isinstance(item['maximum_allowed_fraction'], bool) or not isinstance(item['maximum_allowed_fraction'], (int, float)) for item in qc['reasons']): fail('terminal JSON QC reason types drifted')",
    "sha = re.compile(r'^[0-9a-f]{64}$')",
    "def identity(target):",
    "    resolved = target.resolve(strict=True)",
    "    try: relative = resolved.relative_to(root)",
    "    except ValueError: fail(f'identity path escapes repository: {target}')",
    "    if not resolved.is_file(): fail(f'identity path is not a file: {target}')",
    "    data = resolved.read_bytes()",
    "    return {'path': str(relative), 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}",
    "def validate_record(record, label):",
    "    if not isinstance(record, dict) or set(record) != {'path','bytes','sha256'} or isinstance(record['bytes'], bool) or not isinstance(record['bytes'], int) or record['bytes'] <= 0 or not isinstance(record['sha256'], str) or not sha.fullmatch(record['sha256']): fail(f'malformed {label} identity')",
    "    relative = Path(record['path'])",
    "    if relative.is_absolute() or '..' in relative.parts or record != identity(root / relative): fail(f'{label} identity drifted')",
    "validate_record(payload['source_family_lock'], 'source-family lock')",
    "expected_source_lock = root / 'results' / 'track_b' / 'checkpoints' / 'lava_continuations' / 'source_locks' / f'{source}.json'",
    "if payload['source_family_lock'] != identity(expected_source_lock): fail('source-family lock path drifted')",
    "source_lock = json.loads((root / payload['source_family_lock']['path']).read_text(encoding='utf-8'))",
    "if payload['source_checkpoint_family_sha256'] != source_lock.get('checkpoint_family_sha256') or not isinstance(payload['source_checkpoint_family_sha256'], str) or not sha.fullmatch(payload['source_checkpoint_family_sha256']): fail('source checkpoint-family identity drifted')",
    "validate_record(payload['lineage'], 'continuation lineage')",
    "expected_lineage = root / 'results' / 'track_b' / 'checkpoints' / 'lava_continuations' / 'runs' / continuation / 'lineage.json'",
    "if payload['lineage'] != identity(expected_lineage): fail('continuation lineage path drifted')",
    "validators = payload['validators']",
    "if not isinstance(validators, dict) or set(validators) != {'legacy_scientific_validator','continuation_result_validator'}: fail('terminal validator family drifted')",
    "if validators['legacy_scientific_validator'] != identity(root / 'scripts' / '121_validate_track_b_lava.py') or validators['continuation_result_validator'] != identity(root / 'scripts' / '138_validate_track_b_lava_results_v2.py'): fail('terminal validator identity drifted')",
    "staged = payload['staged_low_level_results']",
    "names = ('lava_locus_status.tsv','lava_univariate.tsv','lava_bivariate.tsv','lava_conditional.tsv')",
    "if not isinstance(staged, list) or staged != [identity(path.parent / name) for name in names]: fail('terminal staged-result identities drifted')",
    "upstream = payload['upstream_continuation_family']",
    "if not isinstance(upstream, dict) or set(upstream) != {'aggregate','conditional_checkpoint_count','conditional_checkpoint_family_sha256','complete'} or upstream['conditional_checkpoint_count'] != 2495 or upstream['complete'] is not True or not isinstance(upstream['conditional_checkpoint_family_sha256'], str) or not sha.fullmatch(upstream['conditional_checkpoint_family_sha256']): fail('terminal upstream-family summary drifted')",
    "aggregate_record = upstream['aggregate']",
    "if not isinstance(aggregate_record, dict) or set(aggregate_record) != {'phase','locus_index','receipt','result','qc'} or aggregate_record['phase'] != 'aggregate-discovery' or aggregate_record['locus_index'] is not None or aggregate_record['qc'] != 'FULL_FAMILY_BH_COMPLETE': fail('terminal aggregate summary drifted')",
    "validate_record(aggregate_record['receipt'], 'aggregate receipt')",
    "validate_record(aggregate_record['result'], 'aggregate result')",
    "print('TRACK_B_LAVA_TERMINAL_QC_JSON_VALID')"
  ), collapse = "\n")
  arguments <- c(
    "-c", shQuote(python), shQuote(path), values$source_fingerprint,
    values$continuation_fingerprint, as.character(as.integer(counts)),
    format(thresholds, digits = 17L, scientific = TRUE, trim = TRUE),
    shQuote(normalizePath(".", mustWork = TRUE))
  )
  check <- suppressWarnings(system2(
    "python3", arguments, stdout = TRUE, stderr = TRUE
  ))
  status <- attr(check, "status")
  if (is.null(status)) status <- 0L
  if (status != 0L || sum(check == "TRACK_B_LAVA_TERMINAL_QC_JSON_VALID") != 1L) {
    stop(paste("Terminal-QC report semantics drifted:", paste(check, collapse = " | ")))
  }
  invisible(TRUE)
}

checkpoint_read <- read_rds_with_sha256(values$rds, "Checkpoint RDS")
checkpoint <- checkpoint_read$value
rm(checkpoint_read)

if (values$phase == "aggregate-discovery") {
  checkpoint <- validate_aggregate_identity(checkpoint)
  checkpoint <- validate_full_aggregate(
    checkpoint, file.path(dirname(values$rds), "conditional_candidates.tsv")
  )
  semantic_locus <- "ALL"
  semantic_chromosome <- "ALL"
  semantic_n_snps <- marker_number(sum(checkpoint$status$n_snps, na.rm = TRUE))
  semantic_pair <- "A;B;CONTROL"
  semantic_qc <- "FULL_FAMILY_BH_COMPLETE"
  semantic_construction <- "TRUE"
} else if (values$phase == "conditional") {
  aggregate_read <- load_aggregate(full_source_validation = FALSE)
  aggregate <- aggregate_read$value
  aggregate_sha256 <- aggregate_read$sha256
  rm(aggregate_read)
  checkpoint <- validate_conditional(
    checkpoint, values$locus_index, aggregate, aggregate_sha256,
    check_local_source = TRUE
  )
  semantic_locus <- as.character(checkpoint$LOC)
  semantic_chromosome <- as.character(checkpoint$CHR)
  semantic_n_snps <- marker_number(checkpoint$n_snps)
  semantic_pair <- as.character(checkpoint$pair)
  semantic_qc <- as.character(checkpoint$state)
  semantic_construction <- if (identical(
    checkpoint$state, "CONDITIONAL_LOCUS_PROCESSED"
  )) "TRUE" else "FALSE"
} else {
  aggregate_read <- load_aggregate(full_source_validation = TRUE)
  aggregate <- aggregate_read$value
  aggregate_sha256 <- aggregate_read$sha256
  rm(aggregate_read)
  if (!all(file.exists(conditional_paths)) ||
      any(!vapply(
        conditional_paths, path_within, logical(1L), root = continuation_root
      ))) {
    stop("Complete canonical conditional family is unavailable")
  }
  conditional_reads <- lapply(seq_along(conditional_paths), function(index) {
    observed <- read_rds_with_sha256(
      conditional_paths[[index]], paste("Conditional checkpoint", index)
    )
    observed$value <- validate_conditional(
      observed$value, index, aggregate, aggregate_sha256, check_local_source = FALSE
    )
    observed
  })
  conditional_hashes <- vapply(
    conditional_reads, function(item) item$sha256, character(1L)
  )
  conditional_checkpoints <- lapply(conditional_reads, function(item) item$value)
  rm(conditional_reads)
  expected_manifest <- data.frame(
    locus_index = seq_len(nrow(loci)),
    LOC = loci$LOC,
    relative_path = conditional_relative_paths,
    sha256 = conditional_hashes,
    stringsAsFactors = FALSE
  )
  final_conditional <- build_final_conditional(conditional_checkpoints)
  attempt_dir <- dirname(values$rds)

  if (values$phase == "finalize") {
    expected_names <- c(
      "schema_version", "source_discovery_fingerprint",
      "continuation_execution_fingerprint", "phase", "locus_index", "LOC", "CHR",
      "aggregate_checkpoint_sha256", "conditional_checkpoint_manifest",
      "staging_validation_sha256", "total_snps", "state"
    )
    if (!is.list(checkpoint) || !identical(names(checkpoint), expected_names) ||
        !identical(checkpoint$schema_version, schema_version) ||
        !identical(checkpoint$source_discovery_fingerprint, values$source_fingerprint) ||
        !identical(
          checkpoint$continuation_execution_fingerprint,
          values$continuation_fingerprint
        ) ||
        !identical(checkpoint$phase, "finalize") ||
        !identical(as.integer(checkpoint$locus_index), 0L) ||
        !identical(checkpoint$LOC, "ALL") ||
        !identical(checkpoint$CHR, "ALL") ||
        !identical(checkpoint$aggregate_checkpoint_sha256, aggregate_sha256) ||
        !identical(checkpoint$state, "STAGED_SEMANTIC_VALIDATION_COMPLETE")) {
      stop("Finalization checkpoint identity drifted")
    }
    validate_manifest_shape(
      checkpoint$conditional_checkpoint_manifest,
      conditional_relative_paths, "Conditional checkpoint"
    )
    if (!identical(
      as.character(checkpoint$conditional_checkpoint_manifest$sha256),
      as.character(expected_manifest$sha256)
    )) {
      stop("Finalization conditional dependency hashes drifted")
    }
    validation_path <- file.path(attempt_dir, "staging_validation.txt")
    if (!file.exists(validation_path) || file.info(validation_path)$size <= 0L ||
        !identical(
          checkpoint$staging_validation_sha256,
          sha256_files(validation_path)[[1L]]
        )) {
      stop("Finalization staged-validation attestation drifted")
    }
    validation_lines <- readLines(validation_path, warn = FALSE)
    if (sum(grepl(
      "^TRACK_B_LAVA_V2_STAGED_RESULTS_SEMANTICALLY_VALIDATED",
      validation_lines
    )) != 1L) {
      stop("Finalization attestation lacks one V2 semantic PASS marker")
    }
    validate_staged_tables(
      attempt_dir, aggregate, final_conditional, require_derived = TRUE
    )
    expected_total_snps <- sum(aggregate$status$n_snps, na.rm = TRUE)
    if (!isTRUE(all.equal(
      as.numeric(checkpoint$total_snps), as.numeric(expected_total_snps),
      tolerance = 0, check.attributes = FALSE
    ))) {
      stop("Finalization total SNP count drifted")
    }
    semantic_qc <- "STAGED_SEMANTIC_VALIDATION_COMPLETE"
  } else {
    expected_names <- c(
      "schema_version", "source_discovery_fingerprint",
      "continuation_execution_fingerprint", "phase", "locus_index", "LOC", "CHR",
      "aggregate_checkpoint_sha256", "conditional_checkpoint_manifest",
      "terminal_qc_sha256", "total_snps", "state"
    )
    if (!is.list(checkpoint) || !identical(names(checkpoint), expected_names) ||
        !identical(checkpoint$schema_version, schema_version) ||
        !identical(checkpoint$source_discovery_fingerprint, values$source_fingerprint) ||
        !identical(
          checkpoint$continuation_execution_fingerprint,
          values$continuation_fingerprint
        ) ||
        !identical(checkpoint$phase, "terminal-qc") ||
        !identical(as.integer(checkpoint$locus_index), 0L) ||
        !identical(checkpoint$LOC, "ALL") ||
        !identical(checkpoint$CHR, "ALL") ||
        !identical(checkpoint$aggregate_checkpoint_sha256, aggregate_sha256) ||
        !identical(checkpoint$state, "TERMINAL_FAILED_QC")) {
      stop("Terminal-QC checkpoint identity drifted")
    }
    validate_manifest_shape(
      checkpoint$conditional_checkpoint_manifest,
      conditional_relative_paths, "Conditional checkpoint"
    )
    if (!identical(
      as.character(checkpoint$conditional_checkpoint_manifest$sha256),
      as.character(expected_manifest$sha256)
    )) {
      stop("Terminal-QC conditional dependency hashes drifted")
    }
    terminal_qc_path <- file.path(attempt_dir, "terminal_qc.json")
    if (!file.exists(terminal_qc_path) || file.info(terminal_qc_path)$size <= 0L ||
        !identical(
          checkpoint$terminal_qc_sha256,
          sha256_files(terminal_qc_path)[[1L]]
        )) {
      stop("Terminal-QC report drifted")
    }
    validate_terminal_qc_json(terminal_qc_path, aggregate, final_conditional)
    validate_staged_tables(
      attempt_dir, aggregate, final_conditional, require_derived = FALSE
    )
    expected_total_snps <- sum(aggregate$status$n_snps, na.rm = TRUE)
    if (!isTRUE(all.equal(
      as.numeric(checkpoint$total_snps), as.numeric(expected_total_snps),
      tolerance = 0, check.attributes = FALSE
    ))) {
      stop("Terminal-QC total SNP count drifted")
    }
    semantic_qc <- "TERMINAL_FAILED_QC"
  }
  semantic_locus <- "ALL"
  semantic_chromosome <- "ALL"
  semantic_n_snps <- marker_number(checkpoint$total_snps)
  semantic_pair <- "A;B;CONTROL"
  semantic_construction <- "TRUE"
}

cat(paste(
  "TRACK_B_LAVA_CONTINUATION_CHECKPOINT",
  paste0("phase=", values$phase),
  paste0("index=", if (is.na(values$locus_index)) 0L else values$locus_index),
  paste0("locus=", semantic_locus),
  paste0("chromosome=", semantic_chromosome),
  paste0("n_snps=", semantic_n_snps),
  paste0("pair=", semantic_pair),
  paste0("qc=", semantic_qc),
  paste0("construction_complete=", semantic_construction),
  paste0("source_fingerprint=", values$source_fingerprint),
  paste0("continuation_fingerprint=", values$continuation_fingerprint),
  "status=PASS",
  sep = "\t"
), "\n", sep = "")
