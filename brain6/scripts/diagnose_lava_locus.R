suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(LAVA))
if (as.character(utils::packageVersion("LAVA")) != "0.1.5") {
  stop("diagnostic requires the frozen LAVA 0.1.5 runtime")
}
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1L) stop("usage: diagnose_lava_locus.R worker_config.json")
cfg <- jsonlite::fromJSON(args[[1]], simplifyVector=FALSE)
loci <- LAVA::read.loci(cfg$loci_file)
locus <- loci[match(as.character(cfg$locus_id), as.character(loci$LOC)), , drop=FALSE]
if (nrow(locus) != 1L || is.na(locus$LOC[[1]])) stop("configured locus not found")
policy <- cfg$execution_policy$locus_processing
for (pair in cfg$pairs) {
  phenos <- as.character(unlist(pair$phenotypes, use.names=FALSE))
  cat("\n=== PAIR", pair$pair_id, "phenotypes", paste(phenos, collapse=","), "===\n")
  input <- tryCatch(LAVA::process.input(input.info.file=pair$input_info,
    sample.overlap.file=pair$sample_overlap_file, ref.prefix=pair$reference_prefix,
    phenos=phenos, input.dir=dirname(pair$input_info)), error=function(e) e)
  if (inherits(input, "error")) {
    cat("process.input error:", conditionMessage(input), "\n")
    next
  }
  captured <- character()
  printed <- capture.output(loc <- tryCatch(withCallingHandlers(
      LAVA::process.locus(locus, input, phenos=phenos,
        min.K=as.integer(policy$min_K), prune.thresh=as.numeric(policy$prune_thresh),
        max.prop.K=as.numeric(policy$max_prop_K), drop.failed=isTRUE(policy$drop_failed),
        max.block.size=as.integer(policy$max_block_size),
        cap.estimates=isTRUE(policy$cap_estimates)),
      message=function(m) { captured <<- c(captured, conditionMessage(m)); invokeRestart("muffleMessage") },
      warning=function(w) { captured <<- c(captured, paste("WARNING:", conditionMessage(w))); invokeRestart("muffleWarning") }),
    error=function(e) e))
  if (length(printed)) cat(paste(printed, collapse="\n"), "\n", sep="")
  if (length(captured)) cat(paste(captured, collapse="\n"), "\n", sep="")
  if (inherits(loc, "error")) {
    cat("process.locus error:", conditionMessage(loc), "\n")
  } else {
    cat("process.locus returned NULL:", is.null(loc), "\n")
    if (!is.null(loc)) cat("retained phenotypes:", paste(loc$phenos, collapse=","), "\n")
  }
}
