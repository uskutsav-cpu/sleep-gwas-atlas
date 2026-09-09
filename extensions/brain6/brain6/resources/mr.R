#!/usr/bin/env Rscript
# Input: local, exposure-selected, LD-clumped instruments. No network lookups.
script <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
need_package("TwoSampleMR")
need(isTRUE(cfg$instruments_ld_clumped) && isTRUE(cfg$instrument_selection_exposure_only),
     "MR instruments must be exposure-selected and LD-clumped, not chosen by PLACO/outcome P")
need(isTRUE(cfg$sample_overlap_reviewed), "Sample-overlap assumptions require review")
x <- read_dt(cfg$exposure); y <- read_dt(cfg$outcome)
cols <- c("SNP","BETA","SE","P","A1","A2","EAF","N")
need(all(cols %in% names(x)) && all(cols %in% names(y)), "MR inputs must use canonical columns")
finite_cols(x,c("BETA","SE","P","N")); finite_cols(y,c("BETA","SE","P","N"))
need(all(x$SE>0 & x$N>0 & x$P>=0 & x$P<=1) && all(y$SE>0 & y$N>0 & y$P>=0 & y$P<=1),
     "Invalid MR SE/N/P values")
need(!anyDuplicated(x$SNP) && !anyDuplicated(y$SNP), "Duplicate MR SNP")
x <- x[P<cfg$instrument_p & (BETA/SE)^2>=cfg$minimum_F]
if (nrow(x)<3L || nrow(y)<3L) {
  fwrite(data.table(SNP=character(),mr_keep=logical()),opath("harmonized.tsv"),sep="\t")
  fwrite(data.table(status="INSUFFICIENT_INSTRUMENTS"),opath("mr.tsv"),sep="\t")
  finish("INSUFFICIENT_EVIDENCE",list(n_instruments=0L,reason="Fewer than three candidate or outcome rows"))
  quit(status=0L)
}
ex <- TwoSampleMR::format_data(as.data.frame(x),type="exposure",snp_col="SNP",beta_col="BETA",
      se_col="SE",effect_allele_col="A1",other_allele_col="A2",eaf_col="EAF",pval_col="P",samplesize_col="N")
oy <- TwoSampleMR::format_data(as.data.frame(y),type="outcome",snp_col="SNP",beta_col="BETA",
      se_col="SE",effect_allele_col="A1",other_allele_col="A2",eaf_col="EAF",pval_col="P",samplesize_col="N")
d <- TwoSampleMR::harmonise_data(ex,oy,action=3)
fwrite(as.data.table(d),opath("harmonized.tsv"),sep="\t")
keep <- d[d$mr_keep,]
if (nrow(keep)<3L) {
  fwrite(data.table(status="INSUFFICIENT_INSTRUMENTS"),opath("mr.tsv"),sep="\t")
  finish("INSUFFICIENT_EVIDENCE",list(n_instruments=nrow(keep)))
} else {
  methods <- c("mr_ivw","mr_weighted_median","mr_egger_regression")
  results <- TwoSampleMR::mr(keep,method_list=methods)
  primary <- results[results$method=="Inverse variance weighted",]
  need(nrow(primary)==1L && is.finite(primary$b) && is.finite(primary$se) && primary$se>0 &&
       is.finite(primary$pval) && primary$pval>=0 && primary$pval<=1, "Invalid native primary IVW estimate")
  fwrite(as.data.table(results),opath("mr.tsv"),sep="\t")
  fwrite(as.data.table(TwoSampleMR::mr_heterogeneity(keep)),opath("heterogeneity.tsv"),sep="\t")
  fwrite(as.data.table(TwoSampleMR::mr_pleiotropy_test(keep)),opath("egger_intercept.tsv"),sep="\t")
  fwrite(as.data.table(TwoSampleMR::mr_leaveoneout(keep)),opath("leave_one_out.tsv"),sep="\t")
  finish(extra=list(n_instruments=nrow(keep),interpretation="DIRECTIONAL_EVIDENCE_CONDITIONAL_ON_IV_ASSUMPTIONS",
         note="No binary causal label from method agreement; family correction is performed after collecting all directions."))
}
