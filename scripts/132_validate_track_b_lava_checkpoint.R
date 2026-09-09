#!/usr/bin/env Rscript

values <- list(phase = NA_character_, locus_index = NA_integer_, fingerprint = NA_character_, rds = NA_character_)
args <- commandArgs(trailingOnly = TRUE)
mapping <- c("--phase" = "phase", "--locus-index" = "locus_index", "--fingerprint" = "fingerprint", "--rds" = "rds")
i <- 1L
while (i <= length(args)) {
  if (!(args[[i]] %in% names(mapping)) || i == length(args)) stop(paste("Unknown/incomplete argument:", args[[i]]))
  values[[unname(mapping[[args[[i]]]])]] <- args[[i + 1L]]
  i <- i + 2L
}
values$locus_index <- as.integer(values$locus_index)
if (!(values$phase %in% c("discovery", "aggregate-discovery", "conditional", "finalize"))) stop("Invalid phase")
if (!grepl("^[0-9a-f]{64}$", values$fingerprint)) stop("Invalid execution fingerprint")
if (is.na(values$rds) || !file.exists(values$rds) || file.info(values$rds)$size <= 0) stop("Missing checkpoint RDS")
if (values$phase %in% c("discovery", "conditional") && (is.na(values$locus_index) || values$locus_index < 1L || values$locus_index > 2495L)) stop("Invalid locus index")
if (!(values$phase %in% c("discovery", "conditional")) && !is.na(values$locus_index)) stop("Aggregate/final phase cannot have a locus index")

locus_lines <- strsplit(readLines("ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile", warn = FALSE), "[[:space:]]+")
locus_header <- locus_lines[[1L]]
if (!identical(locus_header, c("LOC", "CHR", "START", "STOP")) || length(locus_lines) != 2496L) stop("Official locus family drifted")
loci <- as.data.frame(do.call(rbind, locus_lines[-1L]), stringsAsFactors = FALSE)
names(loci) <- locus_header
if (anyDuplicated(loci$LOC)) stop("Official locus identifiers are not unique")

input_info <- read.delim("results/track_b/lava_input_info.tsv", stringsAsFactors = FALSE, check.names = FALSE)
pairs <- read.delim("results/track_b/lava_pair_manifest.tsv", stringsAsFactors = FALSE, check.names = FALSE)
conditioners <- read.delim("results/track_b/local_conditional_manifest.tsv", stringsAsFactors = FALSE, check.names = FALSE)
runtime_table <- read.delim("results/track_b/lava_runtime_policy.tsv", stringsAsFactors = FALSE, check.names = FALSE)
runtime <- setNames(as.list(runtime_table$value), runtime_table$key)
univariate_threshold <- as.numeric(runtime$univariate_p_threshold)
bivariate_fdr_alpha <- as.numeric(runtime$bivariate_fdr_alpha)
if (nrow(input_info) != 8L || anyDuplicated(input_info$phenotype)) stop("Eight-trait input family drifted")
if (!identical(pairs$pair_id, c("A", "B", "CONTROL"))) stop("Three-pair family drifted")
if (!identical(conditioners$conditional_model_id, c("A_BMI_ONLY", "A_SLEEP_APNEA_ONLY", "B_MDD_ONLY", "NONE"))) stop("Conditional model family drifted")
if (!is.finite(univariate_threshold) || !is.finite(bivariate_fdr_alpha)) stop("Runtime thresholds are invalid")

checkpoint <- readRDS(values$rds)
if (!is.list(checkpoint) || !identical(checkpoint$fingerprint, values$fingerprint) || !identical(checkpoint$phase, values$phase)) stop("Checkpoint identity drifted")

exact_names <- function(value, expected, label) {
  if (!is.data.frame(value) || !identical(names(value), expected)) stop(paste(label, "schema drifted"))
}
finite_or_na <- function(value) all(is.na(value) | is.finite(value))
marker_number <- function(value) {
  if (length(value) == 1L && is.finite(value) && value >= 0) as.character(as.integer(value)) else "NA"
}
same_coordinates <- function(table, row) {
  if (!nrow(table)) return(TRUE)
  all(as.character(table$LOC) == row$LOC & as.character(table$CHR) == row$CHR &
    as.character(table$START) == row$START & as.character(table$STOP) == row$STOP)
}

status_pre_fields <- c("LOC", "CHR", "START", "STOP", "status", "n_snps", "K", "univariate_tested", "eligible_bivariate_pairs", "elapsed_seconds")
univ_pre_fields <- c("LOC", "CHR", "START", "STOP", "phen", "h2.obs", "h2.latent", "ascertained", "p", "analysis_status", "error", "n_snps", "K")
bivar_pre_fields <- c("LOC", "CHR", "START", "STOP", "pair_id", "trait1", "trait2", "discovery_rg", "discovery_SE", "discovery_P", "discovery_FDR", "local_covariance", "rho", "rho.lower", "rho.upper", "r2", "r2.lower", "r2.upper", "p", "analysis_status", "error")
conditional_pre_fields <- c("LOC", "CHR", "START", "STOP", "pair_id", "trait1", "trait2", "conditional_model_id", "covariates", "pcor", "ci.lower", "ci.upper", "p", "r2.trait1_z", "r2.trait2_z", "analysis_status", "error")

if (values$phase == "discovery") {
  expected_names <- c("fingerprint", "phase", "locus_index", "LOC", "CHR", "status", "univ", "bivar")
  if (!identical(names(checkpoint), expected_names) || !identical(as.integer(checkpoint$locus_index), values$locus_index)) stop("Discovery checkpoint list drifted")
  official <- loci[values$locus_index, , drop = FALSE]
  if (!identical(as.character(checkpoint$LOC), official$LOC) || !identical(as.character(checkpoint$CHR), official$CHR)) stop("Discovery locus identity drifted")
  exact_names(checkpoint$status, status_pre_fields, "Discovery status")
  exact_names(checkpoint$univ, univ_pre_fields, "Discovery univariate")
  exact_names(checkpoint$bivar, bivar_pre_fields, "Discovery bivariate")
  if (nrow(checkpoint$status) != 1L || !same_coordinates(checkpoint$status, official)) stop("Discovery status coordinates drifted")
  if (nrow(checkpoint$univ) != 8L || !identical(as.character(checkpoint$univ$phen), as.character(input_info$phenotype)) || !same_coordinates(checkpoint$univ, official)) stop("Discovery univariate row family drifted")
  if (!(checkpoint$status$status %in% c("PROCESSED", "PROCESS_FAILED", "UNIVARIATE_FAILED"))) stop("Invalid discovery terminal status")
  if (any(!(checkpoint$univ$analysis_status %in% c("TESTED", "PHENOTYPE_DROPPED", "LOCUS_PROCESS_FAILED", "UNIVARIATE_FAILED")))) stop("Invalid discovery univariate status")
  if (any(!(checkpoint$bivar$pair_id %in% pairs$pair_id)) || anyDuplicated(checkpoint$bivar$pair_id) || !same_coordinates(checkpoint$bivar, official)) stop("Discovery bivariate pair family drifted")
  pair_index <- match(checkpoint$bivar$pair_id, pairs$pair_id)
  if (nrow(checkpoint$bivar) && any(checkpoint$bivar$trait1 != pairs$trait1[pair_index] | checkpoint$bivar$trait2 != pairs$trait2[pair_index])) stop("Discovery bivariate trait identity drifted")
  if (any(!(checkpoint$bivar$analysis_status %in% c("TESTED", "BIVARIATE_FAILED")))) stop("Invalid discovery bivariate status")
  if (!finite_or_na(checkpoint$univ$p) || !finite_or_na(checkpoint$bivar$p)) stop("Non-finite discovery P-value diagnostics")
  tested <- checkpoint$univ$analysis_status == "TESTED"
  if (any(tested & (!is.finite(checkpoint$univ$p) | checkpoint$univ$p < 0 | checkpoint$univ$p > 1 |
      !is.finite(checkpoint$univ$h2.obs) | checkpoint$univ$h2.obs < 0))) stop("A TESTED discovery univariate row is numerically invalid")
  binary_traits <- as.character(input_info$phenotype[!is.na(input_info$cases)])
  tested_binary <- tested & checkpoint$univ$phen %in% binary_traits
  tested_continuous <- tested & !(checkpoint$univ$phen %in% binary_traits)
  if (any(tested_binary & (!is.finite(checkpoint$univ$h2.latent) | checkpoint$univ$h2.latent < 0)) ||
      any(tested_continuous & !is.na(checkpoint$univ$h2.latent))) stop("Discovery latent-scale h2 semantics drifted")
  status_value <- as.character(checkpoint$status$status[[1L]])
  if (status_value == "PROCESS_FAILED" &&
      (!all(checkpoint$univ$analysis_status == "LOCUS_PROCESS_FAILED") || nrow(checkpoint$bivar) != 0L)) stop("PROCESS_FAILED discovery checkpoint has incompatible result rows")
  if (status_value == "UNIVARIATE_FAILED" &&
      (!all(checkpoint$univ$analysis_status == "UNIVARIATE_FAILED") || nrow(checkpoint$bivar) != 0L)) stop("UNIVARIATE_FAILED discovery checkpoint has incompatible result rows")
  if (status_value == "PROCESSED" && any(checkpoint$univ$analysis_status %in% c("LOCUS_PROCESS_FAILED", "UNIVARIATE_FAILED"))) stop("PROCESSED discovery checkpoint contains a locus-level failure row")
  eligible_traits <- checkpoint$univ$phen[tested & checkpoint$univ$p <= univariate_threshold]
  expected_pairs <- pairs$pair_id[pairs$trait1 %in% eligible_traits & pairs$trait2 %in% eligible_traits]
  if (!identical(as.character(checkpoint$bivar$pair_id), as.character(expected_pairs))) stop("Discovery bivariate row family differs from the local-h2 gate")
  tested_bivar <- checkpoint$bivar$analysis_status == "TESTED"
  numeric_bivar <- c("local_covariance", "rho", "rho.lower", "rho.upper", "r2", "r2.lower", "r2.upper", "p")
  if (nrow(checkpoint$bivar) && any(vapply(checkpoint$bivar[numeric_bivar], function(column) any(tested_bivar & !is.finite(column)), logical(1)))) stop("A TESTED bivariate row is numerically invalid")
  if (any(tested_bivar & (checkpoint$bivar$p < 0 | checkpoint$bivar$p > 1))) stop("A TESTED bivariate p-value is out of range")
  if (any(!tested_bivar & (!is.na(checkpoint$bivar$p) | is.na(checkpoint$bivar$error)))) stop("A failed bivariate row is not explicitly preserved")
  semantic_locus <- as.character(checkpoint$LOC)
  semantic_chromosome <- as.character(checkpoint$CHR)
  semantic_n_snps <- marker_number(checkpoint$status$n_snps[[1L]])
  semantic_pair <- if (nrow(checkpoint$bivar)) paste(checkpoint$bivar$pair_id, collapse = ";") else "NONE"
  semantic_qc <- status_value
  semantic_construction <- if (status_value == "PROCESS_FAILED") "FALSE" else "TRUE"
} else if (values$phase == "aggregate-discovery") {
  if (!identical(names(checkpoint), c("fingerprint", "phase", "status", "univ", "bivar"))) stop("Discovery aggregation list drifted")
  exact_names(checkpoint$status, c(status_pre_fields, "analysis_fingerprint", "lava_version", "reference_prefix"), "Aggregated status")
  exact_names(checkpoint$univ, c(univ_pre_fields, "univariate_test_family_n", "univariate_p_threshold", "p_bonferroni", "p_fdr"), "Aggregated univariate")
  exact_names(checkpoint$bivar, c(bivar_pre_fields, "bivariate_test_family_n", "p_fdr", "fdr_significant"), "Aggregated bivariate")
  if (nrow(checkpoint$status) != 2495L || !identical(as.character(checkpoint$status$LOC), loci$LOC)) stop("Aggregated status family drifted")
  expected_univ_loci <- rep(loci$LOC, each = 8L)
  expected_univ_traits <- rep(as.character(input_info$phenotype), times = 2495L)
  if (nrow(checkpoint$univ) != 2495L * 8L || !identical(as.character(checkpoint$univ$LOC), expected_univ_loci) || !identical(as.character(checkpoint$univ$phen), expected_univ_traits)) stop("Aggregated univariate family drifted")
  if (any(checkpoint$univ$univariate_test_family_n != 2495L * 8L) || any(!(checkpoint$bivar$pair_id %in% pairs$pair_id)) || anyDuplicated(paste(checkpoint$bivar$LOC, checkpoint$bivar$pair_id))) stop("Aggregated correction family drifted")
  if (nrow(checkpoint$bivar) && any(checkpoint$bivar$bivariate_test_family_n != nrow(checkpoint$bivar))) stop("Aggregated bivariate denominator drifted")
  if (any(checkpoint$status$analysis_fingerprint != values$fingerprint)) stop("Aggregated fingerprint column drifted")
  tested_univ <- checkpoint$univ$analysis_status == "TESTED"
  univ_family <- rep(1, nrow(checkpoint$univ))
  univ_family[tested_univ] <- checkpoint$univ$p[tested_univ]
  expected_univ_fdr <- p.adjust(univ_family, method = "BH")
  if (any(tested_univ & (!is.finite(checkpoint$univ$p) | checkpoint$univ$p < 0 | checkpoint$univ$p > 1)) ||
      !isTRUE(all.equal(checkpoint$univ$p_fdr[tested_univ], expected_univ_fdr[tested_univ], tolerance = 1e-12)) ||
      any(!tested_univ & !is.na(checkpoint$univ$p_fdr))) stop("Aggregated univariate BH family drifted")
  expected_pairs <- list()
  expected_index <- 0L
  for (locus_id in loci$LOC) {
    local_univ <- checkpoint$univ[as.character(checkpoint$univ$LOC) == locus_id, , drop = FALSE]
    eligible_traits <- local_univ$phen[local_univ$analysis_status == "TESTED" & is.finite(local_univ$p) & local_univ$p <= univariate_threshold]
    local_pairs <- pairs[pairs$trait1 %in% eligible_traits & pairs$trait2 %in% eligible_traits, , drop = FALSE]
    if (nrow(local_pairs)) for (pair_row in seq_len(nrow(local_pairs))) {
      expected_index <- expected_index + 1L
      expected_pairs[[expected_index]] <- c(locus_id, as.character(local_pairs$pair_id[[pair_row]]))
    }
  }
  expected_pair_keys <- if (length(expected_pairs)) vapply(expected_pairs, paste, character(1), collapse = "\t") else character(0)
  observed_pair_keys <- paste(as.character(checkpoint$bivar$LOC), checkpoint$bivar$pair_id, sep = "\t")
  if (!identical(observed_pair_keys, expected_pair_keys)) stop("Aggregated bivariate family differs from all locally eligible pair-locus rows")
  tested_bivar <- checkpoint$bivar$analysis_status == "TESTED"
  bivar_family <- rep(1, nrow(checkpoint$bivar))
  bivar_family[tested_bivar] <- checkpoint$bivar$p[tested_bivar]
  expected_bivar_fdr <- p.adjust(bivar_family, method = "BH")
  if ((nrow(checkpoint$bivar) && any(tested_bivar & (!is.finite(checkpoint$bivar$p) | checkpoint$bivar$p < 0 | checkpoint$bivar$p > 1))) ||
      !isTRUE(all.equal(checkpoint$bivar$p_fdr[tested_bivar], expected_bivar_fdr[tested_bivar], tolerance = 1e-12)) ||
      any(!tested_bivar & !is.na(checkpoint$bivar$p_fdr)) ||
      !identical(as.logical(checkpoint$bivar$fdr_significant), !is.na(checkpoint$bivar$p_fdr) & checkpoint$bivar$p_fdr <= bivariate_fdr_alpha)) stop("Aggregated bivariate BH family drifted")
  candidate_path <- file.path(dirname(values$rds), "conditional_candidates.tsv")
  if (!file.exists(candidate_path) || file.info(candidate_path)$size <= 0) stop("Aggregated checkpoint lacks conditional candidate attestation")
  candidates <- read.delim(candidate_path, stringsAsFactors = FALSE, check.names = FALSE)
  candidate_fields <- c("locus_index", "locus", "chromosome", "start", "stop", "n_snps", "eligible_models", "eligible_pairs")
  if (!identical(names(candidates), candidate_fields)) stop("Conditional candidate schema drifted")
  expected_candidates <- list()
  candidate_index <- 0L
  for (index in seq_len(nrow(loci))) {
    locus_id <- loci$LOC[[index]]
    targets <- checkpoint$bivar[as.character(checkpoint$bivar$LOC) == locus_id & checkpoint$bivar$fdr_significant & checkpoint$bivar$pair_id %in% c("A", "B"), , drop = FALSE]
    eligible_models <- 0L
    eligible_pairs <- character(0)
    if (nrow(targets)) for (target_index in seq_len(nrow(targets))) {
      target <- targets[target_index, , drop = FALSE]
      models <- conditioners[conditioners$pair_id == target$pair_id & conditioners$conditional_model_id != "NONE", , drop = FALSE]
      if (nrow(models)) for (model_index in seq_len(nrow(models))) {
        covars <- strsplit(models$covariates[[model_index]], ";", fixed = TRUE)[[1L]]
        local_univ <- checkpoint$univ[as.character(checkpoint$univ$LOC) == locus_id & checkpoint$univ$phen %in% covars, , drop = FALSE]
        eligible <- length(covars) > 0L && nrow(local_univ) == length(covars) && all(local_univ$analysis_status == "TESTED" & is.finite(local_univ$p) & local_univ$p <= univariate_threshold)
        if (eligible) {
          eligible_models <- eligible_models + 1L
          eligible_pairs <- unique(c(eligible_pairs, as.character(target$pair_id)))
        }
      }
    }
    if (eligible_models > 0L) {
      candidate_index <- candidate_index + 1L
      expected_candidates[[candidate_index]] <- c(as.character(index), locus_id, loci$CHR[[index]], loci$START[[index]], loci$STOP[[index]],
        as.character(checkpoint$status$n_snps[match(locus_id, as.character(checkpoint$status$LOC))]), as.character(eligible_models), paste(eligible_pairs, collapse = ";"))
    }
  }
  expected_candidate_matrix <- if (length(expected_candidates)) do.call(rbind, expected_candidates) else matrix(character(0), nrow = 0L, ncol = length(candidate_fields))
  observed_candidate_matrix <- as.matrix(data.frame(lapply(candidates, as.character), check.names = FALSE, stringsAsFactors = FALSE))
  if (!identical(unname(observed_candidate_matrix), unname(expected_candidate_matrix))) stop("Conditional candidate family differs from frozen full-family gates")
  semantic_locus <- "ALL"
  semantic_chromosome <- "ALL"
  semantic_n_snps <- marker_number(sum(checkpoint$status$n_snps, na.rm = TRUE))
  semantic_pair <- "A;B;CONTROL"
  semantic_qc <- "FULL_FAMILY_BH_COMPLETE"
  semantic_construction <- "TRUE"
} else if (values$phase == "conditional") {
  expected_names <- c("fingerprint", "phase", "locus_index", "LOC", "CHR", "n_snps", "pair", "state", "conditional")
  if (!identical(names(checkpoint), expected_names) || !identical(as.integer(checkpoint$locus_index), values$locus_index)) stop("Conditional checkpoint list drifted")
  official <- loci[values$locus_index, , drop = FALSE]
  if (!identical(as.character(checkpoint$LOC), official$LOC) || !identical(as.character(checkpoint$CHR), official$CHR)) stop("Conditional locus identity drifted")
  exact_names(checkpoint$conditional, conditional_pre_fields, "Conditional result")
  if (!same_coordinates(checkpoint$conditional, official) || any(!(checkpoint$conditional$pair_id %in% c("A", "B"))) || anyDuplicated(paste(checkpoint$conditional$pair_id, checkpoint$conditional$conditional_model_id))) stop("Conditional row identity drifted")
  if (nrow(checkpoint$conditional) && any(!(checkpoint$conditional$conditional_model_id %in% conditioners$conditional_model_id[conditioners$conditional_model_id != "NONE"]))) stop("Conditional model identity drifted")
  if (any(!(checkpoint$conditional$analysis_status %in% c("TESTED", "CONDITIONER_LOCAL_H2_INELIGIBLE", "CONDITIONAL_FAILED", "CONDITIONAL_UNSTABLE_MAX_R2")))) stop("Invalid conditional terminal status")
  if (!(checkpoint$state %in% c("NOT_BIVARIATE_FDR_ELIGIBLE", "NO_FROZEN_MODEL", "CONDITIONER_LOCAL_H2_INELIGIBLE", "CONDITIONAL_LOCUS_PROCESSED", "CONDITIONAL_LOCUS_PROCESS_FAILED"))) stop("Invalid conditional process state")
  family_path <- file.path("results/track_b/checkpoints/lava", values$fingerprint, "aggregate-discovery", "unit_all", "result.rds")
  if (!file.exists(family_path)) stop("Conditional checkpoint lacks its discovery-family dependency")
  family <- readRDS(family_path)
  if (!is.list(family) || !identical(family$fingerprint, values$fingerprint) || !identical(family$phase, "aggregate-discovery")) stop("Conditional discovery-family dependency drifted")
  targets <- family$bivar[as.character(family$bivar$LOC) == official$LOC & family$bivar$fdr_significant & family$bivar$pair_id %in% c("A", "B"), , drop = FALSE]
  expected_models <- character(0)
  expected_pairs <- character(0)
  if (nrow(targets)) for (target_index in seq_len(nrow(targets))) {
    models <- conditioners[conditioners$pair_id == targets$pair_id[[target_index]] & conditioners$conditional_model_id != "NONE", , drop = FALSE]
    expected_models <- c(expected_models, as.character(models$conditional_model_id))
    expected_pairs <- c(expected_pairs, rep(as.character(targets$pair_id[[target_index]]), nrow(models)))
  }
  if (!identical(as.character(checkpoint$conditional$conditional_model_id), expected_models) ||
      !identical(as.character(checkpoint$conditional$pair_id), expected_pairs)) stop("Conditional checkpoint omits or adds a frozen model")
  if (!nrow(targets) && !identical(checkpoint$state, "NOT_BIVARIATE_FDR_ELIGIBLE")) stop("Ineligible conditional locus has the wrong process state")
  tested_conditional <- checkpoint$conditional$analysis_status == "TESTED"
  numeric_conditional <- c("pcor", "ci.lower", "ci.upper", "p", "r2.trait1_z", "r2.trait2_z")
  if (nrow(checkpoint$conditional) && any(vapply(checkpoint$conditional[numeric_conditional], function(column) any(tested_conditional & !is.finite(column)), logical(1)))) stop("A TESTED conditional row is numerically invalid")
  if (any(tested_conditional & (checkpoint$conditional$p < 0 | checkpoint$conditional$p > 1))) stop("A TESTED conditional p-value is out of range")
  semantic_locus <- as.character(checkpoint$LOC)
  semantic_chromosome <- as.character(checkpoint$CHR)
  semantic_n_snps <- marker_number(checkpoint$n_snps)
  semantic_pair <- as.character(checkpoint$pair)
  semantic_qc <- as.character(checkpoint$state)
  semantic_construction <- if (identical(checkpoint$state, "CONDITIONAL_LOCUS_PROCESSED")) "TRUE" else "FALSE"
} else {
  if (!identical(names(checkpoint), c("fingerprint", "phase", "locus_index", "LOC", "CHR", "total_snps", "state")) ||
      !identical(checkpoint$LOC, "ALL") || !identical(checkpoint$CHR, "ALL") ||
      !identical(checkpoint$state, "STAGED_SEMANTIC_VALIDATION_COMPLETE")) stop("Finalization checkpoint drifted")
  staged_names <- c("lava_locus_status.tsv", "lava_univariate.tsv", "lava_bivariate.tsv", "lava_conditional.tsv",
    "04_lava_local_results.tsv", "05_local_conditional_results.tsv")
  staged_paths <- file.path(dirname(values$rds), staged_names)
  if (!all(file.exists(staged_paths)) || any(file.info(staged_paths)$size <= 0)) stop("Finalization checkpoint lacks its complete staged result family")
  semantic_locus <- "ALL"
  semantic_chromosome <- "ALL"
  semantic_n_snps <- marker_number(checkpoint$total_snps)
  semantic_pair <- "A;B;CONTROL"
  semantic_qc <- "STAGED_SEMANTIC_VALIDATION_COMPLETE"
  semantic_construction <- "TRUE"
}

cat(paste("TRACK_B_LAVA_CHECKPOINT_SEMANTIC", paste0("phase=", values$phase),
  paste0("index=", if (is.na(values$locus_index)) 0L else values$locus_index),
  paste0("locus=", semantic_locus), paste0("chromosome=", semantic_chromosome),
  paste0("n_snps=", semantic_n_snps), paste0("pair=", semantic_pair),
  paste0("qc=", semantic_qc), paste0("construction_complete=", semantic_construction),
  paste0("fingerprint=", values$fingerprint), "status=PASS", sep = "\t"), "\n", sep = "")
