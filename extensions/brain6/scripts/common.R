# Shared native-adapter contracts. No automatic package installation here.
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(data.table))
need <- function(ok, message) if (!isTRUE(ok)) stop(message, call.=FALSE)
need_package <- function(pkg) {
  need(requireNamespace(pkg, quietly=TRUE), paste("BLOCKED_BY_SOFTWARE:", pkg))
}
args <- commandArgs(trailingOnly=TRUE)
need(length(args) == 1L, "Usage: Rscript adapter.R settings.json")
cfg <- fromJSON(args[[1]], simplifyVector=TRUE)
out <- Sys.getenv("BRAIN6_JOB_DIR", unset="")
need(nchar(out)>0L && dir.exists(out), "Use brain6 run-job: BRAIN6_JOB_DIR is required")
if (!is.null(cfg$expected_packages)) {
  for (pkg in names(cfg$expected_packages)) {
    need_package(pkg)
    need(as.character(packageVersion(pkg)) == cfg$expected_packages[[pkg]],
         paste("Package version drift:", pkg))
  }
}
set.seed(if (is.null(cfg$seed)) 20260908L else cfg$seed)
# All adapter outputs stay inside the executor's transaction directory.
opath <- function(name) {
  need(!grepl("[/\\\\]", name) && !name %in% c(".", ".."), "Unsafe adapter output filename")
  file.path(out, name)
}
finish <- function(status="PASS", extra=list()) {
  write_json(c(list(status=status),extra), opath("status.json"), auto_unbox=TRUE,
             pretty=TRUE, na="null", digits=NA)
  capture.output(sessionInfo(), file=opath("sessionInfo.txt"))
}
read_dt <- function(path, select=NULL) {
  need(file.exists(path), paste("Missing data:", path))
  # fread gzip support needs R.utils; uncompressed chunks need no extra package.
  if (grepl("\\.(gz|bgz)$",path)) need_package("R.utils")
  fread(path, select=select, na.strings=c("NA","NaN",""), showProgress=FALSE)
}
finite_cols <- function(d, cols) {
  need(all(cols %in% names(d)), "Missing required numeric columns")
  need(all(vapply(d[,..cols], function(x) all(is.finite(x)), logical(1))),
       "Nonfinite adapter input; do not silently na.omit")
}
