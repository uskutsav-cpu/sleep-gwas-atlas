#!/usr/bin/env Rscript
# Diagnostic-only probe of pair cells admitted by the frozen trait-only LAVA gate.
# This script reads frozen inputs and writes only the explicitly supplied new output.
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(LAVA))

fail <- function(message) stop(message, call. = FALSE)
need <- function(test, message) if (!isTRUE(test)) fail(message)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) fail(
  "Usage: Rscript probe_lava_unique_gate_v1.R WORKER_CONFIG FAMILY_LOCK CANONICAL_RESULTS OUTPUT_TSV")
cfg <- fromJSON(args[[1L]], simplifyVector = FALSE)
family <- fromJSON(args[[2L]], simplifyVector = FALSE)
canonical_path <- args[[3L]]
output_path <- args[[4L]]
canonical <- fread(canonical_path)
need(as.character(packageVersion("LAVA")) == "0.1.5", "Pinned LAVA version mismatch")
need(!file.exists(output_path), "Refusing to overwrite diagnostic output")
set.seed(as.integer(cfg$random_seed))
need(cfg$analysis_id == "brain6-lava-local-rg-v2" &&
       family$analysis_id == "brain6-lava-canonical-v3" &&
       setequal(unlist(family$trait_ids), unique(canonical$phen)),
     "Unexpected pair-input or canonical-gate family")

# Reapply exactly the checksum-bound roundoff-only patch used by this source run.
ns <- asNamespace("LAVA")
original <- get("decompose.ld", envir = ns, inherits = FALSE)
body_lines <- deparse(body(original), width.cutoff = 500L)
target <- "M = t(R.base) %*% ld$ld %*% R.base"
hit <- which(trimws(body_lines) == target)
need(length(hit) == 1L, "Pinned LAVA decomposition source changed")
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
unlockBinding("decompose.ld", ns)
assign("decompose.ld", patched, envir = ns)
lockBinding("decompose.ld", ns)

need(nrow(canonical) == family$canonical_univariate$n_tests,
     "Canonical univariate table has unexpected row count")
threshold <- as.numeric(family$canonical_univariate$gate_p_strictly_less_than)
pair_ids <- vapply(cfg$pairs, `[[`, "", "pair_id")
traits_by_pair <- lapply(cfg$pairs, function(pair) as.character(unlist(pair$phenotypes, use.names = FALSE)))
canonical_key <- paste(as.character(canonical$phen), as.character(canonical$locus_id), sep = "\034")
canonical_status <- setNames(as.character(canonical$status), canonical_key)
canonical_p <- setNames(as.numeric(canonical$p), canonical_key)
eligible <- list()
for (j in seq_along(pair_ids)) {
  for (locus_id in as.character(seq_len(family$canonical_univariate$loci))) {
    p <- vapply(traits_by_pair[[j]], function(trait) {
      key <- paste(trait, locus_id, sep = "\034")
      if (is.na(canonical_status[[key]]) || canonical_status[[key]] != "TESTED") return(NA_real_)
      canonical_p[[key]]
    }, numeric(1))
    if (all(is.finite(p)) && all(p < threshold))
      eligible[[length(eligible) + 1L]] <- list(pair_index = j, locus_id = as.character(locus_id), canonical_p = p)
  }
}
eligible_keys <- sort(vapply(eligible, function(x) paste(pair_ids[[x$pair_index]], x$locus_id, sep = "/"), ""))
expected_keys <- sort(c("longsleep__scz/2207", "longsleep__parkinson/2207"))
need(identical(eligible_keys, expected_keys),
     paste("Frozen unique-trait gate candidate set changed; found", paste(eligible_keys, collapse = ",")))

loci <- LAVA::read.loci(cfg$loci_file)
locus <- loci[match("2207", as.character(loci$LOC)), , drop = FALSE]
out <- list()
for (candidate in eligible) {
  pair <- cfg$pairs[[candidate$pair_index]]
  phenotypes <- traits_by_pair[[candidate$pair_index]]
  input <- LAVA::process.input(
    input.info.file = pair$input_info,
    sample.overlap.file = pair$sample_overlap_file,
    ref.prefix = pair$reference_prefix,
    phenos = phenotypes,
    input.dir = dirname(pair$input_info)
  )
  policy <- cfg$execution_policy$locus_processing
  messages <- capture.output(loc <- LAVA::process.locus(
    locus, input, phenos = phenotypes,
    min.K = as.integer(policy$min_K), prune.thresh = as.numeric(policy$prune_thresh),
    max.prop.K = as.numeric(policy$max_prop_K), drop.failed = isTRUE(policy$drop_failed),
    max.block.size = as.integer(policy$max_block_size), cap.estimates = isTRUE(policy$cap_estimates)
  ))
  need(!is.null(loc), paste("Pair-context locus was not estimable:", pair$pair_id))
  context_univ <- as.data.table(LAVA::run.univ(
    loc, phenos = phenotypes, cap.estimates = isTRUE(cfg$execution_policy$univariate$cap_estimates)
  ))
  bivar <- as.data.table(LAVA::run.bivar(
    loc, phenos = phenotypes,
    adap.thresh = as.numeric(unlist(cfg$execution_policy$bivariate$adap_thresh)),
    p.values = isTRUE(cfg$execution_policy$bivariate$p_values),
    CIs = isTRUE(cfg$execution_policy$bivariate$confidence_intervals),
    param.lim = as.numeric(cfg$execution_policy$bivariate$parameter_limit),
    cap.estimates = isTRUE(cfg$execution_policy$bivariate$cap_estimates)
  ))
  need(nrow(context_univ) == 2L && nrow(bivar) == 1L,
       paste("Unexpected LAVA result shape:", pair$pair_id))
  out[[length(out) + 1L]] <- data.table(
    pair_id = pair$pair_id, locus_id = "2207", status = "DIAGNOSTIC_ONLY",
    n_snps = as.integer(loc$n.snps), K = as.integer(loc$K),
    trait1 = phenotypes[[1L]], canonical_trait1_p = candidate$canonical_p[[1L]],
    pair_context_trait1_p = context_univ[phen == phenotypes[[1L]], p],
    trait2 = phenotypes[[2L]], canonical_trait2_p = candidate$canonical_p[[2L]],
    pair_context_trait2_p = context_univ[phen == phenotypes[[2L]], p],
    local_rg = bivar$rho[[1L]], ci_lower = bivar$rho.lower[[1L]],
    ci_upper = bivar$rho.upper[[1L]], p = bivar$p[[1L]],
    reason = "Pair-independent frozen univariate gate; full-family QC and promotion not established"
  )
}
result <- rbindlist(out)
result[, q_bh_full_12475 := p.adjust(replace(rep(1, family$bivariate_family$n_slots),
  seq_len(.N), p), method = "BH")[seq_len(.N)]]
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
fwrite(result, output_path, sep = "\t", na = "NA")
print(result)
