#!/usr/bin/env Rscript

# Additive V2 execution copy of scripts/126_run_track_b_placo_pair.R.
# Scientific alignment, nuisance estimation, PLACO+, failure retention, and BH
# logic are intentionally unchanged.  V2 only binds the admitted worker count
# into every resumable family and requires Python bridge-controlled publication.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: 141_run_track_b_placo_pair_v2.R TASK.tsv --root ROOT --execute [--workers N] [--stage-only]")
task_argument <- args[[1L]]
root_argument <- "."
execute <- FALSE
stage_only <- FALSE
worker_override <- NA_integer_
benchmark_variants <- NA_integer_
i <- 2L
while (i <= length(args)) {
  if (args[[i]] == "--execute") {
    execute <- TRUE; i <- i + 1L; next
  }
  if (args[[i]] == "--stage-only") {
    stage_only <- TRUE; i <- i + 1L; next
  }
  if (args[[i]] %in% c("--root", "--workers", "--benchmark-variants") && i < length(args)) {
    if (args[[i]] == "--root") {
      root_argument <- args[[i + 1L]]
    } else if (args[[i]] == "--workers") {
      worker_override <- suppressWarnings(as.integer(args[[i + 1L]]))
    } else {
      benchmark_variants <- suppressWarnings(as.integer(args[[i + 1L]]))
    }
    i <- i + 2L; next
  }
  stop("unknown or incomplete PLACO+ runner argument: ", args[[i]])
}
if (!execute) stop("explicit --execute is required")
worker_lease_fd <- Sys.getenv("TRACK_B_PLACO_WORKER_LEASE_FD", unset = "")
worker_lease_payload <- Sys.getenv("TRACK_B_PLACO_WORKER_LEASE_PAYLOAD_B64", unset = "")
worker_lease_digest <- Sys.getenv("TRACK_B_PLACO_WORKER_LEASE_PAYLOAD_SHA256", unset = "")
if (!grepl("^[0-9]+$", worker_lease_fd)) {
  stop("PLACO+ V2 runner requires the inherited pair worker lease from script 143")
}
if (!nzchar(worker_lease_payload) || !grepl("^[0-9a-f]{64}$", worker_lease_digest)) {
  stop("PLACO+ V2 runner requires the coordinator-authenticated worker lease payload")
}
if (!requireNamespace("data.table", quietly = TRUE)) stop("pinned R environment must contain data.table")

root <- normalizePath(root_argument, mustWork = TRUE)
task_path <- normalizePath(if (grepl("^/", task_argument)) task_argument else file.path(root, task_argument), mustWork = TRUE)
if (!startsWith(task_path, paste0(root, .Platform$file.sep))) stop("PLACO+ task must be inside the execution root")
options(digits = 17)
runner_started <- proc.time()[["elapsed"]]

sha256 <- function(path) {
  output <- suppressWarnings(system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("shasum failed for ", path)
  if (!length(output)) stop("shasum returned no output for ", path)
  value <- strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", value)) stop("invalid SHA-256 output for ", path)
  value
}

safe_path <- function(relative, must_work = FALSE) {
  if (length(relative) != 1L || is.na(relative) || grepl("^/", relative) || any(strsplit(relative, "/", fixed = TRUE)[[1L]] == "..")) {
    stop("unsafe repository-relative task path: ", relative)
  }
  current <- root
  for (part in strsplit(relative, "/", fixed = TRUE)[[1L]]) {
    current <- file.path(current, part)
    link_target <- Sys.readlink(current)
    if (length(link_target) == 1L && !is.na(link_target) && nzchar(link_target)) {
      stop("repository-relative task path contains a symbolic link: ", relative)
    }
  }
  value <- file.path(root, relative)
  resolved <- normalizePath(value, mustWork = must_work)
  if (!(identical(resolved, root) || startsWith(resolved, paste0(root, .Platform$file.sep)))) {
    stop("repository-relative task path escapes through a symlink: ", relative)
  }
  resolved
}

inode_identity <- function(path) {
  code <- paste0(
    "import os,stat,sys;", "s=os.lstat(sys.argv[1]);",
    "assert stat.S_ISREG(s.st_mode);", "print(f'{s.st_dev}:{s.st_ino}')"
  )
  output <- suppressWarnings(system2(
    "python3", c("-c", shQuote(code), shQuote(path)), stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || length(output) != 1L ||
      !grepl("^[0-9]+:[0-9]+$", output[[1L]])) {
    stop("could not obtain regular-file inode identity for ", path)
  }
  output[[1L]]
}

validate_worker_lease <- function(fd, path, pair_id, run_fingerprint, payload_b64, payload_sha256, coordinator_path) {
  code <- paste0(
    "import base64,fcntl,hashlib,hmac,json,os,re,stat,sys;",
    "fd=int(sys.argv[1]);p,pair,run,b64,want,coord=sys.argv[2:];",
    "raw=base64.b64decode(b64,validate=True);",
    "assert raw.endswith(b'\\n') and hmac.compare_digest(hashlib.sha256(raw).hexdigest(),want);",
    "v=json.loads(raw);",
    "assert set(v)=={'schema_version','pair_id','run_fingerprint','device','inode','authorization_nonce','coordinator_sha256'};",
    "assert v['schema_version']=='sleep-atlas-track-b-placo-worker-lease.1';",
    "assert v['pair_id']==pair and v['run_fingerprint']==run;",
    "assert re.fullmatch(r'[0-9a-f]{64}',v['authorization_nonce']);",
    "assert re.fullmatch(r'[0-9a-f]{64}',v['coordinator_sha256']);",
    "a=os.fstat(fd);b=os.stat(p,follow_symlinks=False);",
    "assert stat.S_ISREG(a.st_mode) and stat.S_ISREG(b.st_mode) and b.st_nlink==1;",
    "assert (a.st_dev,a.st_ino)==(b.st_dev,b.st_ino)==(v['device'],v['inode']);",
    "assert os.pread(fd,a.st_size,0)==raw;",
    "assert hmac.compare_digest(hashlib.sha256(open(coord,'rb').read()).hexdigest(),v['coordinator_sha256']);",
    "probe=os.open(p,os.O_RDWR);",
    "\ntry:\n fcntl.flock(probe,fcntl.LOCK_EX|fcntl.LOCK_NB);raise RuntimeError('unlocked')",
    "\nexcept BlockingIOError:\n print('LOCKED_AUTHORIZED')\nfinally:\n os.close(probe)"
  )
  output <- suppressWarnings(system2(
    "python3", c(
      "-c", shQuote(code), shQuote(fd), shQuote(path), shQuote(pair_id),
      shQuote(run_fingerprint), shQuote(payload_b64), shQuote(payload_sha256),
      shQuote(coordinator_path)
    ), stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || length(output) != 1L || output[[1L]] != "LOCKED_AUTHORIZED") {
    stop("PLACO+ inherited worker lease is not authorized for this exact pair/run and inode")
  }
}

unlink_if_identity <- function(path, expected) {
  observed <- tryCatch(inode_identity(path), error = function(error) NULL)
  if (!is.null(observed) && identical(observed, expected) && !unlink(path)) {
    stop("could not remove exact V2 publication inode: ", path)
  }
}

link_without_replacement <- function(source, destination, expected_source) {
  if (!file.link(source, destination)) return(FALSE)
  observed <- tryCatch(inode_identity(destination), error = function(error) NULL)
  if (is.null(observed) || !identical(observed, expected_source)) {
    unlink_if_identity(destination, expected_source)
    return(FALSE)
  }
  TRUE
}

path_lexists <- function(path) {
  link_target <- Sys.readlink(path)
  file.exists(path) || (length(link_target) == 1L && !is.na(link_target) && nzchar(link_target))
}

require_absent <- function(paths) {
  if (any(vapply(paths, path_lexists, logical(1L)))) {
    stop("refusing to replace an existing V2 temporary artifact")
  }
}

free_bytes <- function(path) {
  output <- system2("df", c("-Pk", path), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || length(output) < 2L) stop("could not determine free storage")
  fields <- strsplit(trimws(output[[length(output)]]), "[[:space:]]+")[[1L]]
  if (length(fields) < 4L) stop("invalid df output")
  value <- suppressWarnings(as.numeric(fields[[4L]]) * 1024)
  if (!is.finite(value) || value < 0) stop("invalid free-storage value")
  value
}

task <- read.delim(task_path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = character())
if (nrow(task) != 1L) stop("PLACO+ task must contain exactly one row")
required_task <- c(
  "analysis_id", "pair_id", "trait1", "trait2", "family_role", "run_fingerprint",
  "aligned_input", "aligned_input_sha256", "aligned_input_rows", "aligned_input_schema",
  "materialization_provenance", "materialization_provenance_sha256", "policy", "policy_sha256",
  "contract_lock", "contract_lock_sha256", "input_gate_lock", "input_gate_lock_sha256",
  "placo_source", "placo_source_sha256", "marginal_p_threshold", "z_squared_maximum",
  "absolute_tolerance", "maximum_numerical_failure_fraction", "primary_headline_threshold",
  "pairwise_gws_threshold", "within_pair_bh_alpha", "minimum_free_bytes",
  "maximum_workers", "pair_concurrency_limit", "shard_size", "checkpoint_dir", "staged_ledger",
  "run_summary", "staged_provenance", "canonical_ledger", "canonical_provenance", "benchmark_raw", "benchmark_out",
  "materializer_sha256", "runner_sha256"
)
if (!identical(names(task), required_task)) stop("PLACO+ task schema differs from the locked execution contract")
expected_pairs <- list(
  A = c("snoring", "parental_lifespan", "PRIMARY_DISCOVERY"),
  B = c("insomnia", "adhd", "PRIMARY_DISCOVERY"),
  CONTROL = c("insomnia", "frailty", "POSITIVE_CONTROL")
)
if (!(task$pair_id %in% names(expected_pairs)) ||
    !identical(c(task$trait1, task$trait2, task$family_role), unname(expected_pairs[[task$pair_id]]))) {
  stop("PLACO+ task is not exactly frozen pair A, B, or CONTROL")
}
if (!grepl("^[0-9a-f]{64}$", task$run_fingerprint)) stop("invalid PLACO+ run fingerprint")
worker_lease_path <- safe_path(file.path(
  "results/track_b/pleiotropy/continuations/post_lava_terminal_v2/sequential_worker_leases",
  paste0(task$pair_id, ".lock")
), TRUE)
worker_coordinator_path <- safe_path("scripts/143_run_track_b_placo_sequential_v2.py", TRUE)
validate_worker_lease(
  worker_lease_fd, worker_lease_path, task$pair_id, task$run_fingerprint,
  worker_lease_payload, worker_lease_digest, worker_coordinator_path
)
if (as.integer(task$pair_concurrency_limit) != 1L) stop("PLACO+ resource envelope must run one pair at a time")

policy_path <- safe_path(task$policy, TRUE)
contract_lock <- safe_path(task$contract_lock, TRUE)
input_gate_lock <- safe_path(task$input_gate_lock, TRUE)
aligned_input <- safe_path(task$aligned_input, TRUE)
materialization_provenance <- safe_path(task$materialization_provenance, TRUE)
placo_source <- safe_path(task$placo_source, TRUE)
runner_path <- safe_path("scripts/141_run_track_b_placo_pair_v2.R", TRUE)
materializer_path <- safe_path("scripts/125_materialize_track_b_placo_pair.py", TRUE)
if (sha256(policy_path) != task$policy_sha256 || sha256(contract_lock) != task$contract_lock_sha256 ||
    sha256(input_gate_lock) != task$input_gate_lock_sha256 ||
    sha256(materialization_provenance) != task$materialization_provenance_sha256 ||
    sha256(aligned_input) != task$aligned_input_sha256 || sha256(placo_source) != task$placo_source_sha256 ||
    sha256(materializer_path) != task$materializer_sha256 || sha256(runner_path) != task$runner_sha256) {
  stop("PLACO+ task input, software, or upstream lock identity drifted")
}
task_sha256 <- sha256(task_path)
input_sha256 <- sha256(aligned_input)
placo_source_sha256 <- sha256(placo_source)
runner_sha256 <- sha256(runner_path)

numeric_field <- function(name, positive = TRUE) {
  value <- suppressWarnings(as.numeric(task[[name]]))
  if (length(value) != 1L || !is.finite(value) || (positive && value <= 0)) stop("invalid locked numeric field: ", name)
  value
}
expected_rows <- as.integer(numeric_field("aligned_input_rows"))
marginal_threshold <- numeric_field("marginal_p_threshold")
z_squared_maximum <- numeric_field("z_squared_maximum")
absolute_tolerance <- numeric_field("absolute_tolerance")
maximum_failure_fraction <- numeric_field("maximum_numerical_failure_fraction", FALSE)
primary_threshold <- numeric_field("primary_headline_threshold")
pairwise_threshold <- numeric_field("pairwise_gws_threshold")
bh_alpha <- numeric_field("within_pair_bh_alpha")
minimum_free <- numeric_field("minimum_free_bytes")
maximum_workers <- as.integer(numeric_field("maximum_workers"))
shard_size <- as.integer(numeric_field("shard_size"))
if (expected_rows <= 0L || maximum_failure_fraction < 0 || maximum_failure_fraction >= 1 ||
    primary_threshold > pairwise_threshold || bh_alpha >= 1 || maximum_workers < 1L || maximum_workers > 4L ||
    shard_size < 1L || shard_size > 100000L) stop("invalid locked PLACO+ execution envelope")
workers <- if (is.na(worker_override)) maximum_workers else worker_override
if (is.na(workers) || workers < 1L || workers > maximum_workers) stop("--workers must be between 1 and the locked maximum")
if (free_bytes(root) < minimum_free) stop("BLOCKED_BY_COMPUTE: free storage is below the locked PLACO+ minimum")

dat <- data.table::fread(aligned_input, showProgress = FALSE, nThread = workers)
required_input <- c("SNP", "CHR", "BP", "A1", "A2", "FRQ1", "FRQ2", "Z1", "Z2", "P1", "P2")
if (!identical(names(dat), required_input) || task$aligned_input_schema != paste(required_input, collapse = ",")) {
  stop("materialized PLACO+ pair schema differs from its lock")
}
if (nrow(dat) != expected_rows || anyDuplicated(dat$SNP) ||
    any(!is.finite(dat$CHR)) || any(dat$CHR < 1 | dat$CHR > 22) ||
    !identical(sort(unique(as.integer(dat$CHR))), 1:22) ||
    any(!is.finite(dat$BP)) || any(dat$BP <= 0) ||
    any(!is.finite(dat$FRQ1)) || any(!is.finite(dat$FRQ2)) ||
    any(dat$FRQ1 <= 0 | dat$FRQ1 >= 1 | dat$FRQ2 <= 0 | dat$FRQ2 >= 1) ||
    any(!is.finite(dat$Z1)) || any(!is.finite(dat$Z2)) ||
    any(dat$Z1^2 > z_squared_maximum) || any(dat$Z2^2 > z_squared_maximum) ||
    any(!is.finite(dat$P1)) || any(!is.finite(dat$P2)) ||
    any(dat$P1 < 0 | dat$P1 > 1 | dat$P2 < 0 | dat$P2 > 1)) {
  stop("materialized PLACO+ pair row family is incomplete, non-autosomal, duplicated, or invalid")
}
if (nrow(dat) > 1L) {
  previous <- seq_len(nrow(dat) - 1L); following <- previous + 1L
  disordered <- dat$CHR[following] < dat$CHR[previous] |
    (dat$CHR[following] == dat$CHR[previous] & dat$BP[following] < dat$BP[previous]) |
    (dat$CHR[following] == dat$CHR[previous] & dat$BP[following] == dat$BP[previous] & dat$SNP[following] <= dat$SNP[previous])
  if (any(disordered)) stop("materialized PLACO+ pair is not in deterministic coordinate order")
}

source(placo_source, local = globalenv())
needed <- c("var.placo", "cor.pearson", "placo.plus")
if (!all(vapply(needed, exists, logical(1L), inherits = TRUE))) stop("pinned PLACO+ source lacks required correlated-test functions")
checkpoint_dir <- safe_path(task$checkpoint_dir, FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
nuisance_path <- file.path(checkpoint_dir, "nuisance.rds")
nuisance_lock_path <- file.path(checkpoint_dir, "nuisance.sha256.tsv")
z_matrix <- as.matrix(dat[, .(Z1, Z2)])
p_matrix <- as.matrix(dat[, .(P1, P2)])
variance_null_rows <- nrow(dat) - sum(dat$P1 < marginal_threshold & dat$P2 < marginal_threshold)
correlation_null_rows <- nrow(dat) - sum(dat$P1 < marginal_threshold | dat$P2 < marginal_threshold)
nuisance_started <- proc.time()[["elapsed"]]
nuisance_reused <- FALSE
if (file.exists(nuisance_path) || file.exists(nuisance_lock_path)) {
  if (!(file.exists(nuisance_path) && file.exists(nuisance_lock_path))) stop("partial PLACO+ nuisance checkpoint family")
  nuisance_lock <- read.delim(nuisance_lock_path, stringsAsFactors = FALSE, check.names = FALSE)
  required_nuisance_lock <- c(
    "run_fingerprint", "input_sha256", "task_sha256", "placo_source_sha256", "runner_sha256",
    "workers", "marginal_p_threshold", "variance_null_rows", "correlation_null_rows", "nuisance_rds_sha256"
  )
  if (nrow(nuisance_lock) != 1L || !identical(names(nuisance_lock), required_nuisance_lock) ||
      nuisance_lock$run_fingerprint != task$run_fingerprint || nuisance_lock$input_sha256 != input_sha256 ||
      nuisance_lock$task_sha256 != task_sha256 || nuisance_lock$placo_source_sha256 != placo_source_sha256 ||
      nuisance_lock$runner_sha256 != runner_sha256 || nuisance_lock$nuisance_rds_sha256 != sha256(nuisance_path) ||
      as.integer(nuisance_lock$workers) != workers ||
      as.numeric(nuisance_lock$marginal_p_threshold) != marginal_threshold ||
      as.integer(nuisance_lock$variance_null_rows) != variance_null_rows ||
      as.integer(nuisance_lock$correlation_null_rows) != correlation_null_rows) {
    stop("stale or tampered PLACO+ nuisance checkpoint family")
  }
  nuisance <- readRDS(nuisance_path)
  required_nuisance <- c(
    "fingerprint", "input_sha256", "task_sha256", "placo_source_sha256", "runner_sha256",
    "workers", "marginal_p_threshold", "var_z", "cor_z", "variance_null_rows", "correlation_null_rows"
  )
  if (!identical(names(nuisance), required_nuisance) ||
      !identical(nuisance$fingerprint, task$run_fingerprint) ||
      nuisance$input_sha256 != input_sha256 || nuisance$task_sha256 != task_sha256 ||
      nuisance$placo_source_sha256 != placo_source_sha256 || nuisance$runner_sha256 != runner_sha256 ||
      !identical(nuisance$workers, workers) ||
      nuisance$marginal_p_threshold != marginal_threshold || nuisance$variance_null_rows != variance_null_rows ||
      nuisance$correlation_null_rows != correlation_null_rows) {
    stop("stale PLACO+ nuisance checkpoint")
  }
  var_z <- nuisance$var_z; cor_z <- nuisance$cor_z
  nuisance_reused <- TRUE
} else {
  var_z <- var.placo(z_matrix, p_matrix, p.threshold = marginal_threshold)
  cor_z <- cor.pearson(z_matrix, p_matrix, p.threshold = marginal_threshold, returnMatrix = FALSE)
  if (length(var_z) != 2L || any(!is.finite(var_z)) || any(var_z <= 0) ||
      length(cor_z) != 1L || !is.finite(cor_z) || cor_z <= -1 || cor_z >= 1) {
    stop("invalid official PLACO+ nuisance estimate")
  }
  nuisance <- list(
    fingerprint = task$run_fingerprint, input_sha256 = input_sha256, task_sha256 = task_sha256,
    placo_source_sha256 = placo_source_sha256, runner_sha256 = runner_sha256,
    workers = workers, marginal_p_threshold = marginal_threshold,
    var_z = as.numeric(var_z), cor_z = as.numeric(cor_z),
    variance_null_rows = variance_null_rows, correlation_null_rows = correlation_null_rows
  )
  temporary <- paste0(nuisance_path, ".", Sys.getpid(), ".tmp")
  temporary_lock <- paste0(nuisance_lock_path, ".", Sys.getpid(), ".tmp")
  require_absent(c(temporary, temporary_lock))
  saveRDS(nuisance, temporary, version = 3)
  nuisance_lock <- data.frame(
    run_fingerprint = task$run_fingerprint, input_sha256 = input_sha256, task_sha256 = task_sha256,
    placo_source_sha256 = placo_source_sha256, runner_sha256 = runner_sha256,
    workers = workers, marginal_p_threshold = marginal_threshold, variance_null_rows = variance_null_rows,
    correlation_null_rows = correlation_null_rows, nuisance_rds_sha256 = sha256(temporary),
    stringsAsFactors = FALSE
  )
  write.table(nuisance_lock, temporary_lock, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t", na = "NA")
  temporary_identity <- inode_identity(temporary)
  temporary_lock_identity <- inode_identity(temporary_lock)
  nuisance_linked <- link_without_replacement(temporary, nuisance_path, temporary_identity)
  lock_linked <- nuisance_linked && link_without_replacement(
    temporary_lock, nuisance_lock_path, temporary_lock_identity
  )
  if (!lock_linked) {
    if (nuisance_linked) unlink_if_identity(nuisance_path, temporary_identity)
    unlink_if_identity(temporary, temporary_identity)
    unlink_if_identity(temporary_lock, temporary_lock_identity)
    stop("could not publish PLACO+ nuisance checkpoint family without replacement")
  }
  unlink_if_identity(temporary, temporary_identity)
  unlink_if_identity(temporary_lock, temporary_lock_identity)
}
nuisance_elapsed <- proc.time()[["elapsed"]] - nuisance_started
if (length(var_z) != 2L || any(!is.finite(var_z)) || any(var_z <= 0) ||
    length(cor_z) != 1L || !is.finite(cor_z) || cor_z <= -1 || cor_z >= 1) {
  stop("invalid official PLACO+ nuisance estimate")
}
rm(z_matrix, p_matrix); invisible(gc())

sanitize_error <- function(value) {
  value <- gsub("[\t\r\n]+", " ", conditionMessage(value))
  if (!nzchar(value)) value <- "UNSPECIFIED_PLACO_PLUS_NUMERICAL_ERROR"
  substr(value, 1L, 500L)
}
evaluate_index <- function(index) {
  tryCatch({
    result <- withCallingHandlers(
      placo.plus(c(dat$Z1[[index]], dat$Z2[[index]]), VarZ = var_z, CorZ = cor_z, AbsTol = absolute_tolerance),
      warning = function(value) stop("PLACO+ numerical warning: ", conditionMessage(value), call. = FALSE)
    )
    statistic <- as.numeric(result$T.placo.plus); p_value <- as.numeric(result$p.placo.plus)
    if (length(statistic) != 1L || length(p_value) != 1L || !is.finite(statistic) || !is.finite(p_value) || p_value < 0 || p_value > 1) {
      stop("non-finite or out-of-range PLACO+ result")
    }
    list(statistic = statistic, p_value = p_value, status = "TESTED", error = NA_character_)
  }, error = function(error) list(
    statistic = NA_real_, p_value = 1, status = "NUMERICAL_FAILURE_P_SET_TO_ONE", error = sanitize_error(error)
  ))
}
nuisance_sha256 <- sha256(nuisance_path)

read_shard <- function(shard) {
  start <- (shard - 1L) * shard_size + 1L; end <- min(nrow(dat), shard * shard_size)
  checkpoint <- file.path(checkpoint_dir, sprintf("shard_%06d.rds", shard))
  checkpoint_lock <- paste0(checkpoint, ".sha256.tsv")
  if (!(file.exists(checkpoint) && file.exists(checkpoint_lock))) {
    if (file.exists(checkpoint) || file.exists(checkpoint_lock)) stop("partial PLACO+ shard checkpoint family: ", checkpoint)
    return(NULL)
  }
  lock <- read.delim(checkpoint_lock, stringsAsFactors = FALSE, check.names = FALSE)
  required_lock <- c(
    "run_fingerprint", "start", "end", "rows", "input_sha256", "task_sha256",
    "nuisance_checkpoint_sha256", "placo_source_sha256", "runner_sha256", "workers", "shard_rds_sha256"
  )
  if (nrow(lock) != 1L || !identical(names(lock), required_lock) ||
      lock$run_fingerprint != task$run_fingerprint || as.integer(lock$start) != start ||
      as.integer(lock$end) != end || as.integer(lock$rows) != end - start + 1L ||
      lock$input_sha256 != input_sha256 || lock$task_sha256 != task_sha256 ||
      lock$nuisance_checkpoint_sha256 != nuisance_sha256 || lock$placo_source_sha256 != placo_source_sha256 ||
      lock$runner_sha256 != runner_sha256 || as.integer(lock$workers) != workers ||
      lock$shard_rds_sha256 != sha256(checkpoint)) {
    stop("stale or tampered PLACO+ shard checkpoint lock: ", checkpoint)
  }
  payload <- readRDS(checkpoint)
  if (!identical(names(payload), c("fingerprint", "input_sha256", "task_sha256", "nuisance_sha256", "workers", "start", "end", "rows")) ||
      !identical(payload$fingerprint, task$run_fingerprint) || payload$input_sha256 != input_sha256 ||
      payload$task_sha256 != task_sha256 || payload$nuisance_sha256 != nuisance_sha256 ||
      !identical(payload$workers, workers) || payload$start != start || payload$end != end ||
      nrow(payload$rows) != end - start + 1L ||
      !identical(names(payload$rows), c("T_PLACO_PLUS", "P_PLACO_PLUS", "analysis_status", "numerical_error")) ||
      any(!is.finite(payload$rows$P_PLACO_PLUS)) || any(payload$rows$P_PLACO_PLUS < 0 | payload$rows$P_PLACO_PLUS > 1) ||
      any(!(payload$rows$analysis_status %in% c("TESTED", "NUMERICAL_FAILURE_P_SET_TO_ONE")))) {
    stop("stale or incomplete PLACO+ shard checkpoint: ", checkpoint)
  }
  tested <- payload$rows$analysis_status == "TESTED"
  failed <- !tested
  if (any(!is.finite(payload$rows$T_PLACO_PLUS[tested])) || any(!is.na(payload$rows$numerical_error[tested])) ||
      any(!is.na(payload$rows$T_PLACO_PLUS[failed])) || any(payload$rows$P_PLACO_PLUS[failed] != 1) ||
      any(is.na(payload$rows$numerical_error[failed])) || any(!nzchar(payload$rows$numerical_error[failed]))) {
    stop("PLACO+ shard checkpoint violates numerical-failure retention semantics: ", checkpoint)
  }
  payload
}

if (!is.na(benchmark_variants)) {
  if (benchmark_variants < 1L) stop("--benchmark-variants must be positive")
  benchmark_n <- min(nrow(dat), benchmark_variants)
  benchmark_indices <- unique(as.integer(round(seq.int(1L, nrow(dat), length.out = benchmark_n))))
  benchmark_started <- proc.time()[["elapsed"]]
  evaluated <- if (workers == 1L) lapply(benchmark_indices, evaluate_index) else
    parallel::mclapply(benchmark_indices, evaluate_index, mc.cores = workers, mc.preschedule = TRUE)
  benchmark_elapsed <- proc.time()[["elapsed"]] - benchmark_started
  if (!is.finite(benchmark_elapsed) || benchmark_elapsed <= 0) stop("PLACO+ benchmark timer resolution was insufficient")
  benchmark_failures <- sum(vapply(evaluated, `[[`, character(1L), "status") == "NUMERICAL_FAILURE_P_SET_TO_ONE")
  rate <- length(benchmark_indices) / benchmark_elapsed
  benchmark_path <- safe_path(task$benchmark_raw, FALSE)
  if (file.exists(benchmark_path)) stop("immutable PLACO+ raw benchmark already exists")
  benchmark <- data.frame(
    analysis_id = task$analysis_id, pair_id = task$pair_id, run_fingerprint = task$run_fingerprint,
    input_rows = nrow(dat), input_sha256 = sha256(aligned_input), task_sha256 = sha256(task_path),
    policy_sha256 = task$policy_sha256, contract_lock_sha256 = task$contract_lock_sha256,
    input_gate_lock_sha256 = task$input_gate_lock_sha256, materializer_sha256 = sha256(materializer_path),
    placo_source_sha256 = sha256(placo_source), nuisance_checkpoint_sha256 = sha256(nuisance_path),
    nuisance_null_rows_variance = variance_null_rows, nuisance_null_rows_correlation = correlation_null_rows,
    VarZ1 = var_z[[1L]], VarZ2 = var_z[[2L]], CorZ = cor_z,
    global_nuisance_input_rows = nrow(dat), global_nuisance_elapsed_seconds = nuisance_elapsed,
    global_nuisance_checkpoint_reused = nuisance_reused,
    benchmark_variants = length(benchmark_indices), workers = workers, numerical_failures = benchmark_failures,
    single_variant_testing_elapsed_seconds = benchmark_elapsed,
    measured_runner_elapsed_seconds = proc.time()[["elapsed"]] - runner_started, variants_per_second = rate,
    projected_full_family_seconds = nrow(dat) / rate, projected_full_family_hours = nrow(dat) / rate / 3600,
    scientific_equivalence = "GLOBAL_OFFICIAL_NUISANCE_ESTIMATED_ON_ALL_VALID_VARIANTS_THEN_IDENTICAL_SINGLE_VARIANT_PLACO_PLUS_CALLS_SAMPLED",
    runner_sha256 = sha256(runner_path), r_version = R.version.string,
    data_table_version = as.character(utils::packageVersion("data.table")), runtime_platform = R.version$platform,
    stringsAsFactors = FALSE
  )
  dir.create(dirname(benchmark_path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(benchmark_path, ".", Sys.getpid(), ".tmp")
  require_absent(temporary)
  write.table(benchmark, temporary, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t", na = "NA")
  temporary_identity <- inode_identity(temporary)
  if (!link_without_replacement(temporary, benchmark_path, temporary_identity)) {
    unlink_if_identity(temporary, temporary_identity)
    stop("could not publish raw PLACO+ benchmark without replacement")
  }
  unlink_if_identity(temporary, temporary_identity)
  message(sprintf(
    "TRACK_B_PLACO_RAW_BENCHMARK pair=%s input=%d benchmark=%d workers=%d rate=%.6g/s projected=%.6gh",
    task$pair_id, nrow(dat), length(benchmark_indices), workers, rate, nrow(dat) / rate / 3600
  ))
  quit(save = "no", status = 0L)
}

started <- proc.time()[["elapsed"]]
shard_count <- ceiling(nrow(dat) / shard_size)
for (shard in seq_len(shard_count)) {
  start <- (shard - 1L) * shard_size + 1L; end <- min(nrow(dat), shard * shard_size)
  checkpoint <- file.path(checkpoint_dir, sprintf("shard_%06d.rds", shard))
  checkpoint_lock <- paste0(checkpoint, ".sha256.tsv")
  observed <- read_shard(shard)
  if (!is.null(observed)) next
  indices <- start:end
  evaluated <- if (workers == 1L) lapply(indices, evaluate_index) else
    parallel::mclapply(indices, evaluate_index, mc.cores = workers, mc.preschedule = TRUE)
  rows <- data.frame(
    T_PLACO_PLUS = vapply(evaluated, `[[`, numeric(1L), "statistic"),
    P_PLACO_PLUS = vapply(evaluated, `[[`, numeric(1L), "p_value"),
    analysis_status = vapply(evaluated, `[[`, character(1L), "status"),
    numerical_error = vapply(evaluated, function(value) if (is.na(value$error)) NA_character_ else value$error, character(1L)),
    stringsAsFactors = FALSE
  )
  payload <- list(
    fingerprint = task$run_fingerprint, input_sha256 = input_sha256, task_sha256 = task_sha256,
    nuisance_sha256 = nuisance_sha256, workers = workers, start = start, end = end, rows = rows
  )
  temporary <- paste0(checkpoint, ".", Sys.getpid(), ".tmp")
  temporary_lock <- paste0(checkpoint_lock, ".", Sys.getpid(), ".tmp")
  require_absent(c(temporary, temporary_lock))
  saveRDS(payload, temporary, version = 3)
  lock <- data.frame(
    run_fingerprint = task$run_fingerprint, start = start, end = end, rows = end - start + 1L,
    input_sha256 = input_sha256, task_sha256 = task_sha256, nuisance_checkpoint_sha256 = nuisance_sha256,
    placo_source_sha256 = placo_source_sha256, runner_sha256 = runner_sha256,
    workers = workers, shard_rds_sha256 = sha256(temporary), stringsAsFactors = FALSE
  )
  write.table(lock, temporary_lock, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t", na = "NA")
  temporary_identity <- inode_identity(temporary)
  temporary_lock_identity <- inode_identity(temporary_lock)
  checkpoint_linked <- link_without_replacement(temporary, checkpoint, temporary_identity)
  lock_linked <- checkpoint_linked && link_without_replacement(
    temporary_lock, checkpoint_lock, temporary_lock_identity
  )
  if (!lock_linked) {
    if (checkpoint_linked) unlink_if_identity(checkpoint, temporary_identity)
    unlink_if_identity(temporary, temporary_identity)
    unlink_if_identity(temporary_lock, temporary_lock_identity)
    stop("could not publish PLACO+ shard checkpoint family without replacement")
  }
  unlink_if_identity(temporary, temporary_identity)
  unlink_if_identity(temporary_lock, temporary_lock_identity)
  message(sprintf("PLACO+ %s: checkpointed shard %d/%d (%d variants)", task$pair_id, shard, shard_count, end))
}

p_values <- numeric(nrow(dat)); failure_count <- 0L
for (shard in seq_len(shard_count)) {
  payload <- read_shard(shard)
  start <- (shard - 1L) * shard_size + 1L; end <- min(nrow(dat), shard * shard_size)
  if (is.null(payload)) stop("PLACO+ shard family changed during finalization")
  p_values[start:end] <- payload$rows$P_PLACO_PLUS
  failure_count <- failure_count + sum(payload$rows$analysis_status == "NUMERICAL_FAILURE_P_SET_TO_ONE")
}
failure_fraction <- failure_count / nrow(dat)
if (failure_fraction > maximum_failure_fraction) {
  stop(sprintf("PLACO+ numerical failure fraction %.12g exceeds locked maximum %.12g", failure_fraction, maximum_failure_fraction))
}
bh_q <- p.adjust(p_values, method = "BH", n = nrow(dat))
primary_count <- if (task$pair_id == "CONTROL") 0L else sum(p_values <= primary_threshold)
pairwise_count <- sum(p_values <= pairwise_threshold)
bh_count <- sum(bh_q <= bh_alpha)
terminal_status <- if (primary_count > 0L || pairwise_count > 0L || bh_count > 0L) "COMPLETE_WITH_HITS" else "TESTED_NO_HIT"

ledger <- safe_path(task$staged_ledger, FALSE)
summary_path <- safe_path(task$run_summary, FALSE)
ledger_fields <- c(
  "analysis_id", "pair_id", "SNP", "CHR", "BP", "A1", "A2", "Z1", "Z2", "P1", "P2",
  "T_PLACO_PLUS", "P_PLACO_PLUS", "PLACO_BH_Q", "within_pair_family_n", "analysis_status", "numerical_error"
)
summary_fields <- c(
  "analysis_id", "pair_id", "run_fingerprint", "aligned_rows", "within_pair_family_n",
  "nuisance_null_rows_variance", "nuisance_null_rows_correlation", "VarZ1", "VarZ2", "CorZ",
  "nuisance_checkpoint_sha256",
  "numerical_failure_count", "numerical_failure_fraction", "minimum_p_placo_plus",
  "primary_headline_count", "pairwise_gws_count", "within_pair_bh_count", "terminal_status",
  "shard_count", "shard_size", "workers", "elapsed_seconds", "ledger", "ledger_bytes",
  "ledger_sha256", "ledger_schema", "task_sha256", "aligned_sha256", "placo_source_sha256",
  "materializer_sha256", "runner_sha256", "r_version", "data_table_version", "runtime_platform"
)
if (file.exists(ledger) || file.exists(summary_path)) {
  if (!(file.exists(ledger) && file.exists(summary_path))) stop("partial staged PLACO+ result family exists")
  existing_summary <- read.delim(summary_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(existing_summary) != 1L || !identical(names(existing_summary), summary_fields) ||
      existing_summary$analysis_id != task$analysis_id || existing_summary$pair_id != task$pair_id ||
      existing_summary$run_fingerprint != task$run_fingerprint || existing_summary$task_sha256 != task_sha256 ||
      existing_summary$aligned_sha256 != input_sha256 || existing_summary$placo_source_sha256 != placo_source_sha256 ||
      existing_summary$runner_sha256 != runner_sha256 || existing_summary$ledger != task$staged_ledger ||
      existing_summary$ledger_schema != paste(ledger_fields, collapse = ",") ||
      as.numeric(existing_summary$ledger_bytes) != file.info(ledger)$size ||
      existing_summary$ledger_sha256 != sha256(ledger) || as.integer(existing_summary$aligned_rows) != nrow(dat) ||
      as.integer(existing_summary$within_pair_family_n) != nrow(dat) ||
      as.integer(existing_summary$numerical_failure_count) != failure_count ||
      as.integer(existing_summary$primary_headline_count) != primary_count ||
      as.integer(existing_summary$pairwise_gws_count) != pairwise_count ||
      as.integer(existing_summary$within_pair_bh_count) != bh_count ||
      as.integer(existing_summary$workers) != workers ||
      existing_summary$terminal_status != terminal_status) stop("stale staged PLACO+ result family")
} else {
  dir.create(dirname(ledger), recursive = TRUE, showWarnings = FALSE)
  temporary_ledger <- paste0(ledger, ".", Sys.getpid(), ".tmp")
  temporary_summary <- paste0(summary_path, ".", Sys.getpid(), ".tmp")
  require_absent(c(temporary_ledger, temporary_summary))
  connection <- gzfile(temporary_ledger, open = "wt", compression = 6)
  writeLines(paste(ledger_fields, collapse = "\t"), connection)
  for (shard in seq_len(shard_count)) {
    payload <- read_shard(shard)
    if (is.null(payload)) stop("PLACO+ shard disappeared during staged-ledger construction")
    start <- payload$start; end <- payload$end
    out <- data.frame(
      analysis_id = rep(task$analysis_id, end - start + 1L), pair_id = rep(task$pair_id, end - start + 1L),
      SNP = dat$SNP[start:end], CHR = dat$CHR[start:end], BP = dat$BP[start:end],
      A1 = dat$A1[start:end], A2 = dat$A2[start:end], Z1 = dat$Z1[start:end], Z2 = dat$Z2[start:end],
      P1 = dat$P1[start:end], P2 = dat$P2[start:end], T_PLACO_PLUS = payload$rows$T_PLACO_PLUS,
      P_PLACO_PLUS = payload$rows$P_PLACO_PLUS, PLACO_BH_Q = bh_q[start:end],
      within_pair_family_n = rep(nrow(dat), end - start + 1L),
      analysis_status = payload$rows$analysis_status, numerical_error = payload$rows$numerical_error,
      stringsAsFactors = FALSE
    )
    write.table(out, connection, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA")
  }
  close(connection)
  elapsed <- proc.time()[["elapsed"]] - started
  ledger_bytes <- file.info(temporary_ledger)$size
  ledger_sha256 <- sha256(temporary_ledger)
  summary <- data.frame(
    analysis_id = task$analysis_id, pair_id = task$pair_id, run_fingerprint = task$run_fingerprint,
    aligned_rows = nrow(dat), within_pair_family_n = nrow(dat),
    nuisance_null_rows_variance = variance_null_rows, nuisance_null_rows_correlation = correlation_null_rows,
    VarZ1 = var_z[[1L]], VarZ2 = var_z[[2L]], CorZ = cor_z, nuisance_checkpoint_sha256 = nuisance_sha256,
    numerical_failure_count = failure_count, numerical_failure_fraction = failure_fraction,
    minimum_p_placo_plus = min(p_values), primary_headline_count = primary_count,
    pairwise_gws_count = pairwise_count, within_pair_bh_count = bh_count, terminal_status = terminal_status,
    shard_count = shard_count, shard_size = shard_size, workers = workers, elapsed_seconds = elapsed,
    ledger = task$staged_ledger, ledger_bytes = ledger_bytes, ledger_sha256 = ledger_sha256,
    ledger_schema = paste(ledger_fields, collapse = ","),
    task_sha256 = task_sha256, aligned_sha256 = input_sha256,
    placo_source_sha256 = placo_source_sha256, materializer_sha256 = sha256(materializer_path),
    runner_sha256 = runner_sha256, r_version = R.version.string,
    data_table_version = as.character(utils::packageVersion("data.table")), runtime_platform = R.version$platform,
    stringsAsFactors = FALSE
  )
  if (!identical(names(summary), summary_fields)) stop("internal PLACO+ summary schema drifted")
  write.table(summary, temporary_summary, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t", na = "NA")
  temporary_ledger_identity <- inode_identity(temporary_ledger)
  temporary_summary_identity <- inode_identity(temporary_summary)
  ledger_linked <- link_without_replacement(temporary_ledger, ledger, temporary_ledger_identity)
  summary_linked <- ledger_linked && link_without_replacement(
    temporary_summary, summary_path, temporary_summary_identity
  )
  if (!summary_linked) {
    if (ledger_linked) unlink_if_identity(ledger, temporary_ledger_identity)
    unlink_if_identity(temporary_ledger, temporary_ledger_identity)
    unlink_if_identity(temporary_summary, temporary_summary_identity)
    stop("could not publish staged PLACO+ result family without replacement")
  }
  unlink_if_identity(temporary_ledger, temporary_ledger_identity)
  unlink_if_identity(temporary_summary, temporary_summary_identity)
}

if (stage_only) {
  message(sprintf("TRACK_B_PLACO_STAGE_COMPLETE pair=%s rows=%d failures=%d status=%s", task$pair_id, nrow(dat), failure_count, terminal_status))
  quit(save = "no", status = 0L)
}
stop("V2 PLACO+ runner is stage-only; scripts/140_materialize_track_b_placo_pair_v2.py must revalidate and publish")
