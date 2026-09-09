#!/usr/bin/env Rscript
# One locus per job. Low local h2 is NOT an execution failure.
script <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
need_package("LAVA")
need(cfg$mode=="locus", "Only locus mode is supported")
if (is.null(cfg$sample_overlap_file)) {
  need(isTRUE(cfg$no_overlap_reviewed), "Missing overlap matrix is not evidence of no sample overlap")
}
input <- LAVA::process.input(input.info.file=cfg$input_info,
          sample.overlap.file=cfg$sample_overlap_file,ref.prefix=cfg$reference_prefix,phenos=cfg$phenotypes)
loci <- LAVA::read.loci(cfg$loci_file)
row <- which(as.character(loci$LOC)==as.character(cfg$locus_id))
need(length(row)==1L, "Locus ID is absent or duplicated")
locus <- LAVA::process.locus(loci[row,],input)
need(!is.null(locus), "FAILED_QC: locus processing returned NULL")
u <- as.data.table(LAVA::run.univ(locus))
need(nrow(u)==length(cfg$phenotypes) && all(is.finite(u$p)), "Invalid/missing local univariate tests")
fwrite(u,opath("univariate.tsv"),sep="\t")
pass <- all(u$p < cfg$univariate_p_threshold)
if (pass) {
  b <- as.data.table(LAVA::run.bivar(locus,phenos=cfg$phenotypes))
  need(nrow(b)>=1L && all(is.finite(b$p)), "Invalid bivariate result")
  fwrite(b,opath("bivariate.tsv"),sep="\t")
  finish(extra=list(locus=cfg$locus_id,bivariate_eligible=TRUE,n_univariate=nrow(u)))
} else {
  fwrite(data.table(status="NOT_TESTED_LOW_LOCAL_H2"),opath("bivariate.tsv"),sep="\t")
  finish("INSUFFICIENT_EVIDENCE",list(locus=cfg$locus_id,bivariate_eligible=FALSE,
           reason="LOW_LOCAL_H2_NOT_EXECUTION_FAILURE",n_univariate=nrow(u)))
}
