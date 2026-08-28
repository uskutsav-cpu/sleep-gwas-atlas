#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || args[[3L]] != "--execute") {
  stop("usage: 46_run_placo_pair.R TASK.tsv TASK.lock.tsv --execute")
}
if (!requireNamespace("data.table", quietly = TRUE) ||
    !requireNamespace("R.utils", quietly = TRUE)) {
  stop("pinned R environment must contain data.table and R.utils")
}

task_path <- normalizePath(args[[1L]], mustWork = TRUE)
lock_path <- normalizePath(args[[2L]], mustWork = TRUE)
root <- normalizePath(".", mustWork = TRUE)

sha256 <- function(path) {
  output <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("shasum failed for ", path)
  value <- strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", value)) stop("invalid SHA-256 output for ", path)
  value
}

lock <- read.delim(lock_path, stringsAsFactors = FALSE, check.names = FALSE)
task <- read.delim(task_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(lock) != 1L || nrow(task) != 1L) stop("task and task lock must each contain one row")
if (lock$task_sha256 != sha256(task_path)) stop("PLACO+ task differs from its lock")
if (lock$analysis_id != task$analysis_id || lock$pair_id != task$pair_id) {
  stop("PLACO+ task identity differs from its lock")
}

required <- c(
  "analysis_id", "pair_id", "pair_input", "pair_input_sha256", "pair_input_rows",
  "placo_source", "placo_source_sha256", "marginal_p_threshold",
  "z_squared_maximum", "abs_tolerance", "maximum_failure_fraction",
  "conventional_threshold", "family_threshold", "variant_hits_out", "summary_out"
)
if (!all(required %in% names(task))) stop("PLACO+ task lacks required columns")
resolve <- function(path) normalizePath(file.path(root, path), mustWork = TRUE)
pair_input <- resolve(task$pair_input)
placo_source <- resolve(task$placo_source)
if (sha256(pair_input) != task$pair_input_sha256) stop("materialized pair checksum drifted")
if (sha256(placo_source) != task$placo_source_sha256) stop("pinned PLACO+ source checksum drifted")

numeric_field <- function(name, positive = TRUE) {
  value <- suppressWarnings(as.numeric(task[[name]]))
  if (length(value) != 1L || !is.finite(value) || (positive && value <= 0)) {
    stop("invalid locked numeric field: ", name)
  }
  value
}
marginal_threshold <- numeric_field("marginal_p_threshold")
z_squared_maximum <- numeric_field("z_squared_maximum")
abs_tolerance <- numeric_field("abs_tolerance")
maximum_failure_fraction <- numeric_field("maximum_failure_fraction", positive = FALSE)
conventional_threshold <- numeric_field("conventional_threshold")
family_threshold <- numeric_field("family_threshold")
expected_rows <- as.integer(numeric_field("pair_input_rows"))
if (maximum_failure_fraction < 0 || maximum_failure_fraction >= 1) {
  stop("invalid locked numerical-failure fraction")
}
if (family_threshold >= conventional_threshold) {
  stop("family threshold must be stricter than the conventional threshold")
}

dat <- data.table::fread(
  pair_input,
  select = c("SNP", "CHR", "BP", "A1", "A2", "Z1", "Z2", "P1", "P2"),
  showProgress = FALSE
)
required_input <- c("SNP", "CHR", "BP", "A1", "A2", "Z1", "Z2", "P1", "P2")
if (!identical(names(dat), required_input)) stop("materialized pair schema differs from lock")
if (nrow(dat) != expected_rows || nrow(dat) == 0L) stop("materialized pair row count differs from lock")
if (anyDuplicated(dat$SNP)) stop("materialized pair contains duplicate SNPs")
if (any(!is.finite(dat$Z1)) || any(!is.finite(dat$Z2)) ||
    any(!is.finite(dat$P1)) || any(!is.finite(dat$P2)) ||
    any(dat$P1 <= 0 | dat$P1 > 1) || any(dat$P2 <= 0 | dat$P2 > 1)) {
  stop("materialized pair contains invalid Z/P values")
}
if (any(dat$Z1^2 > z_squared_maximum) || any(dat$Z2^2 > z_squared_maximum)) {
  stop("materialized pair retained a variant above the locked squared-Z limit")
}
if (nrow(dat) > 1L) {
  previous <- seq_len(nrow(dat) - 1L)
  following <- previous + 1L
  disordered <- dat$CHR[following] < dat$CHR[previous] |
    (dat$CHR[following] == dat$CHR[previous] & dat$BP[following] < dat$BP[previous]) |
    (dat$CHR[following] == dat$CHR[previous] & dat$BP[following] == dat$BP[previous] &
       dat$SNP[following] <= dat$SNP[previous])
  if (any(disordered)) stop("materialized pair is not in deterministic coordinate order")
}

source(placo_source, local = globalenv())
needed <- c("var.placo", "cor.pearson", "placo.plus")
if (!all(vapply(needed, exists, logical(1L), inherits = TRUE))) {
  stop("pinned PLACO+ source lacks required functions")
}
z_matrix <- as.matrix(dat[, .(Z1, Z2)])
p_matrix <- as.matrix(dat[, .(P1, P2)])
var_z <- var.placo(z_matrix, p_matrix, p.threshold = marginal_threshold)
cor_z <- cor.pearson(z_matrix, p_matrix, p.threshold = marginal_threshold, returnMatrix = FALSE)
if (length(var_z) != 2L || any(!is.finite(var_z)) || any(var_z <= 0)) {
  stop("invalid PLACO+ nuisance variance estimate")
}
if (length(cor_z) != 1L || !is.finite(cor_z) || cor_z <= -1 || cor_z >= 1) {
  stop("invalid PLACO+ nuisance correlation estimate")
}

hit_rows <- list()
hit_count <- 0L
failure_count <- 0L
minimum_p <- 1
chunk_size <- 10000L
for (start in seq.int(1L, nrow(dat), by = chunk_size)) {
  end <- min(nrow(dat), start + chunk_size - 1L)
  evaluated <- vapply(start:end, function(index) {
    tryCatch({
      result <- placo.plus(
        z_matrix[index, ], VarZ = var_z, CorZ = cor_z, AbsTol = abs_tolerance
      )
      p_value <- as.numeric(result$p.placo.plus)
      statistic <- as.numeric(result$T.placo.plus)
      if (!is.finite(p_value) || p_value < 0 || p_value > 1 || !is.finite(statistic)) {
        stop("non-finite PLACO+ result")
      }
      c(statistic = statistic, p_value = p_value, complete = 1)
    }, error = function(error) c(statistic = NA_real_, p_value = NA_real_, complete = 0))
  }, numeric(3L))
  failure_count <- failure_count + sum(evaluated["complete", ] == 0)
  complete_p <- evaluated["p_value", evaluated["complete", ] == 1]
  if (length(complete_p)) minimum_p <- min(minimum_p, complete_p)
  selected <- which(evaluated["complete", ] == 1 &
                    evaluated["p_value", ] <= conventional_threshold)
  if (length(selected)) {
    indices <- start:end
    for (offset in selected) {
      index <- indices[[offset]]
      hit_count <- hit_count + 1L
      hit_rows[[hit_count]] <- data.frame(
        analysis_id = task$analysis_id, pair_id = task$pair_id,
        sleep_trait = task$sleep_trait, non_sleep_trait = task$non_sleep_trait,
        analysis_tier = task$analysis_tier,
        SNP = dat$SNP[[index]], CHR = dat$CHR[[index]], BP = dat$BP[[index]],
        A1 = dat$A1[[index]], A2 = dat$A2[[index]],
        Z1 = dat$Z1[[index]], Z2 = dat$Z2[[index]],
        P1 = dat$P1[[index]], P2 = dat$P2[[index]],
        T_PLACO_PLUS = evaluated["statistic", offset],
        P_PLACO_PLUS = evaluated["p_value", offset],
        conventional_significant = TRUE,
        locked_family_significant = evaluated["p_value", offset] <= family_threshold,
        stringsAsFactors = FALSE
      )
    }
  }
  if (end %% 100000L == 0L || end == nrow(dat)) {
    message(sprintf("PLACO+ %s: evaluated %d/%d variants", task$pair_id, end, nrow(dat)))
  }
}

failure_fraction <- failure_count / nrow(dat)
if (failure_fraction > maximum_failure_fraction) {
  stop(sprintf(
    "PLACO+ numerical failure fraction %.8g exceeds locked maximum %.8g",
    failure_fraction, maximum_failure_fraction
  ))
}
hit_columns <- c(
  "analysis_id", "pair_id", "sleep_trait", "non_sleep_trait", "analysis_tier",
  "SNP", "CHR", "BP", "A1", "A2", "Z1", "Z2", "P1", "P2",
  "T_PLACO_PLUS", "P_PLACO_PLUS", "conventional_significant", "locked_family_significant"
)
if (length(hit_rows)) {
  hits <- do.call(rbind, hit_rows)
} else {
  hits <- as.data.frame(setNames(replicate(length(hit_columns), character(), simplify = FALSE), hit_columns))
}
hits <- hits[, hit_columns, drop = FALSE]
summary <- data.frame(
  analysis_id = task$analysis_id, pair_id = task$pair_id,
  sleep_trait = task$sleep_trait, non_sleep_trait = task$non_sleep_trait,
  analysis_tier = task$analysis_tier, evaluated_variant_count = nrow(dat),
  numerical_failure_count = failure_count, numerical_failure_fraction = failure_fraction,
  VarZ1 = var_z[[1L]], VarZ2 = var_z[[2L]], CorZ = cor_z,
  marginal_p_threshold = marginal_threshold, z_squared_maximum = z_squared_maximum,
  conventional_threshold = conventional_threshold, family_threshold = family_threshold,
  conventional_hit_count = sum(hits$conventional_significant == TRUE),
  locked_family_hit_count = sum(hits$locked_family_significant == TRUE),
  minimum_p_placo_plus = minimum_p, analysis_status = "PLACO_PLUS_COMPLETE",
  task_sha256 = sha256(task_path), input_sha256 = sha256(pair_input),
  placo_source_sha256 = sha256(placo_source), stringsAsFactors = FALSE
)

publish_table <- function(table, relative_path) {
  destination <- file.path(root, relative_path)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(destination, ".tmp")
  on.exit(unlink(temporary), add = TRUE)
  write.table(
    table, temporary, quote = FALSE, row.names = FALSE, col.names = TRUE,
    sep = "\t", na = "NA"
  )
  if (!file.rename(temporary, destination)) stop("could not atomically publish ", destination)
}
publish_table(hits, task$variant_hits_out)
publish_table(summary, task$summary_out)
message(sprintf(
  "PLACO_PLUS_PAIR_OK pair=%s variants=%d failures=%d conventional_hits=%d family_hits=%d",
  task$pair_id, nrow(dat), failure_count,
  summary$conventional_hit_count, summary$locked_family_hit_count
))
