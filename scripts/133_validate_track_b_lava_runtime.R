#!/usr/bin/env Rscript

# Fail closed on the exact scientific implementation used by Track B.  The
# package version alone is insufficient: the preregistered official commit has
# a process.locus implementation that fixes the reference-sample correction
# internally, and another implementation carrying the same version string
# would change local h2, covariance, eligibility, and downstream test families.

expected_versions <- c(
  LAVA = "0.1.5",
  matrixsampling = "2.0.0",
  cpp11 = "0.5.2",
  keep = "1.0"
)

if (!identical(as.character(getRversion()), "4.3.3")) {
  stop(paste("Track B LAVA requires R 4.3.3; found", getRversion()))
}
for (package in names(expected_versions)) {
  if (!requireNamespace(package, quietly = TRUE)) stop(paste("Missing pinned package:", package))
  observed <- as.character(utils::packageVersion(package))
  if (!identical(observed, unname(expected_versions[[package]]))) {
    stop(paste("Pinned package version drifted:", package, observed))
  }
}

required_functions <- c("process.input", "read.loci", "process.locus", "run.univ", "run.bivar", "run.pcor")
namespace <- asNamespace("LAVA")
if (any(!vapply(required_functions, exists, logical(1), envir = namespace, inherits = FALSE))) {
  stop("The pinned LAVA namespace lacks a required Track B function")
}

process_locus <- get("process.locus", envir = namespace, inherits = FALSE)
expected_formals <- c(
  "locus", "input", "phenos", "min.K", "prune.thresh", "max.prop.K",
  "drop.failed", "max.block.size", "cap.estimates"
)
if (!identical(names(formals(process_locus)), expected_formals)) {
  stop("LAVA process.locus formals differ from official commit e729a245")
}
body_text <- paste(deparse(body(process_locus), width.cutoff = 500L), collapse = "\n")
normalized_body <- gsub("[[:space:]]+", "", body_text)
if (!startsWith(normalized_body, "{nref.correction=F") ||
    !grepl("loc\\$nref.scale=ifelse\\(nref.correction,", normalized_body) ||
    !grepl("loc\\$omega=\\(t\\(loc\\$delta\\)%\\*%loc\\$delta/loc\\$K-loc\\$sigma\\)\\*loc\\$nref.scale", normalized_body)) {
  stop("LAVA process.locus reference-sample correction semantics drifted from the pinned commit")
}

cat(paste(
  "TRACK_B_LAVA_RUNTIME",
  "status=PASS",
  "r_version=4.3.3",
  paste0("platform=", R.version$platform),
  "lava_version=0.1.5",
  "lava_commit=e729a245f7b6923967a96804fbf5246eadf2d6c6",
  "process_locus_nref_correction=FALSE",
  paste0("functions=", paste(required_functions, collapse = ";")),
  sep = "\t"
), "\n", sep = "")
