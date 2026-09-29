#!/usr/bin/env Rscript
# Frozen trait-only comparison of two predeclared models for the Dashti release.
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(LAVA))

fail <- function(message) stop(message, call. = FALSE)
need <- function(ok, message) if (!isTRUE(ok)) fail(message)
args <- commandArgs(trailingOnly = TRUE)
need(length(args) == 1L, 'Usage: Rscript run_longsleep_linear_pilot_v1.R WORKER.json')
cfg <- fromJSON(args[[1L]], simplifyVector = FALSE)
need(cfg$analysis_id == 'brain6_longsleep_linear_88_pilot_v1', 'Unexpected analysis ID')
need(cfg$arm_id %in% c('linear', 'binary'), 'Unexpected model arm')
need(as.character(packageVersion('LAVA')) == '0.1.5', 'LAVA version mismatch')
need(!file.exists(cfg$output_tsv) && !file.exists(cfg$summary_json), 'Refusing to overwrite pilot outputs')
need(length(cfg$locus_ids) == as.integer(cfg$expected_loci), 'Worker locus count changed')
need(length(unique(unlist(cfg$locus_ids))) == as.integer(cfg$expected_loci), 'Duplicate worker locus')

# Exact roundoff-only numerical patch used by the immutable canonical v3 run.
ns <- asNamespace('LAVA')
original <- get('decompose.ld', envir = ns, inherits = FALSE)
body_lines <- deparse(body(original), width.cutoff = 500L)
hit <- which(trimws(body_lines) == 'M = t(R.base) %*% ld$ld %*% R.base')
need(length(hit) == 1L, 'Pinned LAVA block-decomposition source changed')
insert <- c('            asymmetry <- max(abs(M - t(M)))',
            '            tolerance <- 100 * .Machine$double.eps * max(1, nrow(M)) * max(1, max(abs(M)))',
            "            if (!is.finite(asymmetry) || asymmetry > tolerance) stop('Material asymmetry in block-reduced LD matrix; refusing numerical symmetrization', call. = FALSE)",
            '            if (asymmetry > 0) M <- (M + t(M)) / 2')
body_lines <- append(body_lines, insert, after = hit)
patched <- original
body(patched) <- parse(text = paste(body_lines, collapse = '\n'))[[1L]]
environment(patched) <- ns
unlockBinding('decompose.ld', ns)
assign('decompose.ld', patched, envir = ns)
lockBinding('decompose.ld', ns)

set.seed(as.integer(cfg$random_seed))
loci <- LAVA::read.loci(cfg$loci_file)
need(length(unique(loci$LOC)) == 88L, 'Prospective locus file changed')
policy <- cfg$execution_policy$locus_processing
univariate <- cfg$execution_policy$univariate
checkref <- get('check.reference', envir = ns, inherits = FALSE)
loadref <- get('load.reference', envir = ns, inherits = FALSE)
readss <- get('read.sumstats.file', envir = ns, inherits = FALSE)
harm <- get('harmonize.snps', envir = ns, inherits = FALSE)
align <- get('align', envir = ns, inherits = FALSE)
threshold <- as.numeric(cfg$strict_gate_p)
rows <- vector('list', length(cfg$locus_ids))
current_chromosome <- NA_integer_
reference <- NULL

for (index in seq_along(cfg$locus_ids)) {
  locus_id <- as.character(cfg$locus_ids[[index]])
  loc_index <- match(locus_id, as.character(loci$LOC))
  need(!is.na(loc_index), paste('Locus missing:', locus_id))
  input_count <- as.integer(cfg$input_rows[[locus_id]])
  need(length(input_count) == 1L && !is.na(input_count) && input_count >= 0L,
       paste('Missing or invalid input row count:', locus_id))
  if (input_count == 0L) {
    rows[[index]] <- data.table(trait_id = 'longsleep', locus_id = locus_id,
      status = 'NOT_RUN', n_snps = NA_integer_, n_components = NA_integer_,
      h2_obs = NA_real_, p = NA_real_, strict_gate_pass = FALSE,
      reason = 'NO_SOURCE_VERIFIED_SNPS')
    if (index %% 5L == 0L || index == length(rows))
      cat(sprintf('PILOT_PROGRESS arm=%s worker=%d done=%d total=%d\n',
                  cfg$arm_id, as.integer(cfg$worker_id), index, length(rows)))
    next
  }
  chromosome <- as.integer(loci$CHR[[loc_index]])
  if (is.na(current_chromosome) || chromosome != current_chromosome) {
    reference <- loadref(checkref(paste0(cfg$reference_prefix_stem, chromosome)))
    current_chromosome <- chromosome
  }
  input_dir <- file.path(cfg$input_root, paste0('locus_', locus_id))
  info <- fread(file.path(input_dir, 'input_info.tsv'))
  need(nrow(info) == 1L && info$phenotype[[1L]] == 'longsleep', 'Input-info identity mismatch')
  if (cfg$arm_id == 'linear') {
    need(is.na(info$cases[[1L]]) && is.na(info$controls[[1L]]), 'Linear arm is not declared continuous')
  } else {
    need(info$cases[[1L]] == 34184 && info$controls[[1L]] == 305742, 'Binary arm counts changed')
  }
  info$filename <- file.path(input_dir, info$filename)
  info$N <- info$cases + info$controls
  info$prop_cases <- info$cases / info$N
  info$binary <- !is.na(info$prop_cases) & info$prop_cases != 1
  result <- tryCatch({
    input <- new.env(parent = globalenv())
    input$info <- copy(info)
    input$P <- 1L
    input$reference <- reference
    input$sample.overlap <- NULL
    input$sum.stats <- list(longsleep = readss(info$filename[[1L]], 'longsleep'))
    need(all(input$sum.stats$longsleep$N == 339926), 'Pilot N proxy changed')
    harm(input)
    align(input)
    messages <- capture.output(locus <- LAVA::process.locus(
      loci[loc_index, , drop = FALSE], input, phenos = 'longsleep',
      min.K = as.integer(policy$min_K),
      prune.thresh = as.numeric(policy$prune_thresh),
      max.prop.K = as.numeric(policy$max_prop_K),
      drop.failed = isTRUE(policy$drop_failed),
      max.block.size = as.integer(policy$max_block_size),
      cap.estimates = isTRUE(policy$cap_estimates)))
    if (is.null(locus)) {
      reason <- if (any(grepl('Fewer than', messages, fixed = TRUE)))
        'FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS' else if
        (any(grepl('Negative variance estimate', messages, fixed = TRUE)))
        'LOW_LOCAL_H2_UNDERPOWERED' else 'LAVA_PROCESS_LOCUS_RETURNED_NULL'
      data.table(trait_id = 'longsleep', locus_id = locus_id, status = 'NOT_RUN',
                 n_snps = NA_integer_, n_components = NA_integer_, h2_obs = NA_real_,
                 p = NA_real_, strict_gate_pass = FALSE, reason = reason)
    } else {
      uni <- as.data.table(LAVA::run.univ(locus, phenos = 'longsleep',
                                         cap.estimates = isTRUE(univariate$cap_estimates)))
      need(nrow(uni) == 1L && is.finite(uni$p[[1L]]) &&
           uni$p[[1L]] >= 0 && uni$p[[1L]] <= 1, 'Invalid univariate result')
      data.table(trait_id = 'longsleep', locus_id = locus_id, status = 'TESTED',
                 n_snps = as.integer(locus$n.snps), n_components = as.integer(locus$K),
                 h2_obs = as.numeric(uni$h2.obs[[1L]]), p = as.numeric(uni$p[[1L]]),
                 strict_gate_pass = as.numeric(uni$p[[1L]]) < threshold, reason = '')
    }
  }, error = function(error) data.table(
    trait_id = 'longsleep', locus_id = locus_id, status = 'FAILED',
    n_snps = NA_integer_, n_components = NA_integer_, h2_obs = NA_real_,
    p = NA_real_, strict_gate_pass = FALSE,
    reason = substr(conditionMessage(error), 1L, 500L)))
  rows[[index]] <- result
  if (index %% 5L == 0L || index == length(rows))
    cat(sprintf('PILOT_PROGRESS arm=%s worker=%d done=%d total=%d\n',
                cfg$arm_id, as.integer(cfg$worker_id), index, length(rows)))
}

output <- rbindlist(rows, use.names = TRUE)
need(nrow(output) == as.integer(cfg$expected_loci) &&
     uniqueN(output$locus_id) == as.integer(cfg$expected_loci), 'Incomplete worker result')
write.table(output, cfg$output_tsv, sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
summary <- list(analysis_id = cfg$analysis_id, arm_id = cfg$arm_id,
                worker_id = cfg$worker_id, rows = nrow(output),
                tested = sum(output$status == 'TESTED'),
                not_run = sum(output$status == 'NOT_RUN'),
                failed = sum(output$status == 'FAILED'),
                strict_gate_pass = sum(output$strict_gate_pass), threshold = threshold)
writeLines(toJSON(summary, auto_unbox = TRUE, pretty = TRUE), cfg$summary_json)
cat(sprintf('PILOT_COMPLETE arm=%s worker=%d tested=%d not_run=%d failed=%d\n',
            cfg$arm_id, as.integer(cfg$worker_id), summary$tested, summary$not_run, summary$failed))
if (summary$failed > 0L) fail('Worker has failed loci; inspect result before adjudication')
