#!/usr/bin/env Rscript

# Extract one complete official-block signed Pearson correlation matrix from
# the pinned LAVA UKB v1.1 EUR GRCh37 reference.  This helper never squares,
# thresholds, thins, partitions, or lead-centres the LD matrix.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 14L) {
  stop(paste(
    "usage: 149_extract_track_b_signed_ld.R REFERENCE_PREFIX CHROMOSOME",
    paste(
      "VARIANT_ORDER SIGNED_LD_OUT QC_OUT MIN_REFERENCE_N REFERENCE_ID EXTRACTION_MODE",
      "RUNTIME_LOCK RUNTIME_LOCK_SHA256 RUNTIME_ATTESTATION_OUT ANALYSIS_ID",
      "LOCUS_ENTRY_ID EXECUTION_AMENDMENT_SHA256"
    )
  ))
}

reference_prefix <- args[[1L]]
chromosome <- suppressWarnings(as.integer(args[[2L]]))
order_path <- args[[3L]]
ld_out <- args[[4L]]
qc_out <- args[[5L]]
minimum_reference_n <- suppressWarnings(as.integer(args[[6L]]))
reference_id <- args[[7L]]
extraction_mode <- args[[8L]]
runtime_lock_path <- args[[9L]]
runtime_lock_sha256 <- args[[10L]]
runtime_attestation_out <- args[[11L]]
analysis_id <- args[[12L]]
locus_entry_id <- args[[13L]]
execution_amendment_sha256 <- args[[14L]]
if (!identical(extraction_mode, "SIGNED_PEARSON_FROM_LAVA_READ_LD")) {
  stop("signed-LD extraction mode is not the frozen Pearson-correlation mode")
}
if (
  is.na(chromosome) || chromosome < 1L || chromosome > 22L ||
  is.na(minimum_reference_n) || minimum_reference_n != 10000L ||
  !identical(reference_id, "LAVA_UKB_v1.1_EUR_GRCh37") ||
  !identical(analysis_id, "track-b-v1.0-finemapping-trait-coloc") ||
  !grepl("^[A-Za-z0-9_.-]+$", locus_entry_id) ||
  !grepl("^[0-9a-f]{64}$", execution_amendment_sha256) ||
  !grepl("^[0-9a-f]{64}$", runtime_lock_sha256)
) {
  stop("chromosome, analysis/reference identity, or immutable hash drifted")
}
if (
  !file.exists(order_path) || !file.exists(runtime_lock_path) ||
  file.exists(ld_out) || file.exists(qc_out) || file.exists(runtime_attestation_out)
) {
  stop("variant order is absent or a no-replace output already exists")
}

sha256_file <- function(path) {
  command <- Sys.which("sha256sum")
  command_args <- path
  if (!nzchar(command)) {
    command <- Sys.which("shasum")
    command_args <- c("-a", "256", path)
  }
  if (!nzchar(command)) stop("no SHA-256 utility is available")
  output <- system2(command, command_args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("SHA-256 utility failed")
  token <- strsplit(trimws(output[[1L]]), "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", token)) stop("SHA-256 utility returned an invalid digest")
  token
}

if (sha256_file(runtime_lock_path) != runtime_lock_sha256) {
  stop("signed-LD runtime lock hash differs before extraction")
}

pinned_library <- normalizePath(".r-env/lib/R/library", winslash = "/", mustWork = TRUE)
.libPaths(pinned_library)
if (!identical(normalizePath(.libPaths(), winslash = "/", mustWork = TRUE), pinned_library)) {
  stop("R library search path is not the one pinned repository library")
}
suppressPackageStartupMessages(loadNamespace("data.table", lib.loc = pinned_library))
suppressPackageStartupMessages(library(LAVA, lib.loc = pinned_library))
if (as.character(packageVersion("LAVA", lib.loc = pinned_library)) != "0.1.5") {
  stop("LAVA runtime differs from locked version 0.1.5")
}

runtime_lock_fields <- c(
  "record_type", "name", "version", "path", "file_count", "bytes", "sha256"
)
runtime_namespace_closure <- c(
  "base", "compiler", "data.table", "datasets", "graphics", "grDevices",
  "LAVA", "methods", "stats", "tools", "utils"
)

package_tree_identity <- function(package_path) {
  relative <- sort(list.files(
    package_path, all.files = TRUE, recursive = TRUE, full.names = FALSE,
    include.dirs = FALSE, no.. = TRUE
  ), method = "radix")
  full <- file.path(package_path, relative)
  info <- file.info(full)
  keep <- !is.na(info$isdir) & !info$isdir
  relative <- relative[keep]
  full <- full[keep]
  info <- info[keep, , drop = FALSE]
  if (!length(relative) || any(grepl("[\t\r\n]", relative))) {
    stop("signed-LD runtime package tree is empty or has an unsafe path")
  }
  links <- Sys.readlink(full)
  if (any(!is.na(links) & nzchar(links))) {
    stop("signed-LD runtime package tree contains a symbolic link")
  }
  hashes <- vapply(full, sha256_file, character(1L))
  sizes <- format(info$size, scientific = FALSE, trim = TRUE)
  manifest <- paste0(paste(relative, sizes, hashes, sep = "\t", collapse = "\n"), "\n")
  temporary <- tempfile(pattern = ".ld-runtime-tree.", tmpdir = dirname(ld_out))
  on.exit(unlink(temporary), add = TRUE)
  connection <- file(temporary, open = "wb")
  writeBin(charToRaw(enc2utf8(manifest)), connection)
  close(connection)
  list(
    file_count = length(relative), bytes = sum(info$size), sha256 = sha256_file(temporary)
  )
}

verify_runtime_lock <- function(phase) {
  if (sha256_file(runtime_lock_path) != runtime_lock_sha256) {
    stop(sprintf("%s signed-LD runtime lock hash differs", phase))
  }
  lock <- read.delim(
    runtime_lock_path, stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = "character"
  )
  if (
    !identical(names(lock), runtime_lock_fields) ||
    nrow(lock) != length(runtime_namespace_closure) + 1L
  ) stop("signed-LD runtime lock schema/closure count drifted")
  executable <- lock[lock$record_type == "R_EXECUTABLE", , drop = FALSE]
  packages <- lock[lock$record_type == "PACKAGE_TREE", , drop = FALSE]
  if (
    nrow(executable) != 1L || executable$name != "Rscript" ||
    executable$version != "R-4.3.3" || executable$file_count != "1" ||
    normalizePath(executable$path, winslash = "/", mustWork = TRUE) !=
      normalizePath(".r-env/bin/Rscript", winslash = "/", mustWork = TRUE) ||
    as.character(file.info(executable$path)$size) != executable$bytes ||
    sha256_file(executable$path) != executable$sha256 ||
    !identical(packages$name, runtime_namespace_closure)
  ) stop(sprintf("%s signed-LD executable or namespace family differs from lock", phase))
  if (!identical(sort(loadedNamespaces()), sort(runtime_namespace_closure))) {
    stop(sprintf("%s actually loaded signed-LD namespace closure differs from lock", phase))
  }
  for (index in seq_len(nrow(packages))) {
    row <- packages[index, , drop = FALSE]
    expected_path <- file.path(pinned_library, row$name)
    located <- find.package(row$name, lib.loc = pinned_library, quiet = TRUE)
    loaded_path <- if (row$name == "base") {
      file.path(R.home("library"), "base")
    } else {
      getNamespaceInfo(asNamespace(row$name), "path")
    }
    if (
      !nzchar(located) || normalizePath(located, winslash = "/", mustWork = TRUE) !=
        normalizePath(expected_path, winslash = "/", mustWork = TRUE) ||
      normalizePath(loaded_path, winslash = "/", mustWork = TRUE) !=
        normalizePath(row$path, winslash = "/", mustWork = TRUE) ||
      normalizePath(row$path, winslash = "/", mustWork = TRUE) !=
        normalizePath(expected_path, winslash = "/", mustWork = TRUE) ||
      as.character(utils::packageDescription(
        row$name, lib.loc = pinned_library, fields = "Version"
      )) != row$version
    ) stop(sprintf("%s loaded-package path/version mismatch for %s", phase, row$name))
    tree <- package_tree_identity(row$path)
    if (
      as.character(tree$file_count) != row$file_count ||
      as.character(tree$bytes) != row$bytes || tree$sha256 != row$sha256
    ) stop(sprintf("%s installed-package byte mismatch for %s", phase, row$name))
  }
  invisible(lock)
}

verify_runtime_lock("PRE")

order <- read.delim(order_path, stringsAsFactors = FALSE, check.names = FALSE)
required_order <- c("SNP", "CHR", "BP", "A1", "A2")
if (!identical(names(order), required_order) || nrow(order) < 50L) {
  stop("variant-order schema/count drifted")
}
order$SNP <- tolower(as.character(order$SNP))
order$A1 <- toupper(as.character(order$A1))
order$A2 <- toupper(as.character(order$A2))
if (
  anyDuplicated(order$SNP) || any(!grepl("^rs[0-9]+$", order$SNP)) ||
  any(order$CHR != chromosome) || any(order$BP < 1L) ||
  any(!order$A1 %in% c("A", "C", "G", "T")) ||
  any(!order$A2 %in% c("A", "C", "G", "T")) ||
  any(order$A1 == order$A2) ||
  any(paste0(order$A1, order$A2) %in% c("AT", "TA", "CG", "GC"))
) {
  stop("variant order contains duplicate, mismatched, or palindromic identities")
}

chromosome_prefix <- paste0(reference_prefix, "_chr", chromosome)
reference <- LAVA:::load.reference(LAVA:::check.reference(chromosome_prefix))
if (!is.finite(reference$sample.size) || reference$sample.size < minimum_reference_n) {
  stop("LAVA reference sample size is below 10,000")
}
loaded <- LAVA:::read.ld(reference, order$SNP, require.freq = TRUE)
if (
  !identical(loaded$mode, "ld") || !is.matrix(loaded$ld) ||
  !identical(rownames(loaded$ld), order$SNP) ||
  !identical(colnames(loaded$ld), order$SNP) ||
  nrow(loaded$info) != nrow(order) ||
  !identical(as.character(loaded$info$SNP), order$SNP) ||
  !identical(as.integer(loaded$info$CHR), as.integer(order$CHR)) ||
  !identical(as.integer(loaded$info$POS), as.integer(order$BP)) ||
  !identical(toupper(as.character(loaded$info$A1)), order$A1) ||
  !identical(toupper(as.character(loaded$info$A2)), order$A2) ||
  any(!is.finite(loaded$info$NOBS))
) {
  stop("LAVA signed-LD SNP/order/coordinate/allele family drifted")
}
reference_sample_size <- min(as.numeric(loaded$info$NOBS))
if (!is.finite(reference_sample_size) || reference_sample_size < minimum_reference_n) {
  stop("included LAVA reference variants have sample size below 10,000")
}
R <- loaded$ld
storage.mode(R) <- "double"
if (nrow(R) != nrow(order) || ncol(R) != nrow(order) || any(!is.finite(R))) {
  stop("LAVA signed-LD matrix is nonfinite or not square/full-sized")
}
symmetry <- max(abs(R - t(R)))
diagonal <- max(abs(diag(R) - 1))
maximum <- max(abs(R))
if (symmetry > 1e-8 || diagonal > 1e-6 || maximum > 1 + 1e-6) {
  stop("LAVA signed-LD finite/symmetry/diagonal/range QC failed")
}
minimum_eigenvalue <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
if (!is.finite(minimum_eigenvalue) || minimum_eigenvalue < -1e-6) {
  stop("LAVA signed-LD matrix fails exact complete-matrix PSD tolerance")
}

dir.create(dirname(ld_out), recursive = TRUE, showWarnings = FALSE)
ld_temporary <- paste0(ld_out, ".", Sys.getpid(), ".tmp")
qc_temporary <- paste0(qc_out, ".", Sys.getpid(), ".tmp")
runtime_attestation_temporary <- paste0(runtime_attestation_out, ".", Sys.getpid(), ".tmp")
if (
  file.exists(ld_temporary) || file.exists(qc_temporary) ||
  file.exists(runtime_attestation_temporary)
) stop("temporary output collision")
connection <- gzfile(ld_temporary, open = "wt", compression = 6)
tryCatch(
  write.table(
    data.frame(SNP = rownames(R), R, check.names = FALSE), connection,
    quote = FALSE, row.names = FALSE, sep = "\t", na = "NA", digits = 17
  ),
  finally = close(connection)
)
qc <- data.frame(
  reference_id = reference_id,
  reference_sample_size = format(reference_sample_size, scientific = FALSE, trim = TRUE),
  variant_count = nrow(R),
  symmetry_max_abs = format(symmetry, digits = 17, scientific = TRUE),
  diagonal_max_abs = format(diagonal, digits = 17, scientific = TRUE),
  minimum_eigenvalue = format(minimum_eigenvalue, digits = 17, scientific = TRUE),
  maximum_abs_correlation = format(maximum, digits = 17, scientific = TRUE),
  matrix_kind = "SIGNED_PEARSON_CORRELATION",
  variant_order_sha256 = sha256_file(order_path),
  status = "PASS",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write.table(
  qc, qc_temporary, quote = FALSE, row.names = FALSE, sep = "\t", na = "NA",
  col.names = TRUE
)
verify_runtime_lock("POST")
runtime_attestation <- data.frame(
  schema_version = "sleep-atlas-track-b-signed-ld-runtime-attestation.1",
  analysis_id = analysis_id, locus_entry_id = locus_entry_id,
  reference_id = reference_id,
  execution_amendment_sha256 = execution_amendment_sha256,
  runtime_package_lock_sha256 = runtime_lock_sha256,
  precheck_status = "PASS_EXACT_LOADED_CLOSURE_PATHS_AND_BYTES",
  postcheck_status = "PASS_EXACT_LOADED_CLOSURE_PATHS_AND_BYTES",
  loaded_namespace_count = length(loadedNamespaces()),
  loaded_namespace_closure = paste(runtime_namespace_closure, collapse = ","),
  lava_version = as.character(packageVersion("LAVA", lib.loc = pinned_library)),
  stringsAsFactors = FALSE, check.names = FALSE
)
write.table(
  runtime_attestation, runtime_attestation_temporary, quote = FALSE,
  row.names = FALSE, sep = "\t", na = "NA", col.names = TRUE
)
if (!file.link(ld_temporary, ld_out)) stop("could not exclusively install signed-LD output")
unlink(ld_temporary)
if (!file.link(qc_temporary, qc_out)) stop("could not exclusively install signed-LD QC output")
unlink(qc_temporary)
if (!file.link(runtime_attestation_temporary, runtime_attestation_out)) {
  stop("could not exclusively install signed-LD runtime attestation")
}
unlink(runtime_attestation_temporary)

cat(sprintf(
  paste0(
    "TRACK_B_SIGNED_PEARSON_LD_OK chromosome=%d variants=%d ",
    "reference_n=%s min_eigenvalue=%.17g\n"
  ), chromosome, nrow(R), format(reference_sample_size, scientific = FALSE),
  minimum_eigenvalue
))
rm(R, loaded, reference)
invisible(gc(full = TRUE))
