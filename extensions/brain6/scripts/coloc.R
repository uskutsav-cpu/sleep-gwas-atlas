#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
need_package("coloc")
a <- readRDS(cfg$fit1); b <- readRDS(cfg$fit2)
need(isTRUE(a$converged) && isTRUE(b$converged), "Both fine-maps must have converged")
need(all(c("p1","p2","p12") %in% names(cfg)), "Explicit coloc priors required")
need(cfg$p1>0 && cfg$p1<1 && cfg$p2>0 && cfg$p2<1 && cfg$p12>0 && cfg$p12<min(cfg$p1,cfg$p2), "Invalid prior probabilities")
if (!length(a$sets$cs) || !length(b$sets$cs)) {
  fwrite(data.table(message="One or both traits have no credible set"),opath("coloc.tsv"),sep="\t")
  saveRDS(NULL,opath("coloc.rds"))
  finish("INSUFFICIENT_EVIDENCE",list(reason="NO_CREDIBLE_SET_IN_ONE_OR_BOTH_TRAITS"))
} else {
  res <- coloc::coloc.susie(a,b,p1=cfg$p1,p2=cfg$p2,p12=cfg$p12)
  need(!is.null(res$summary), "No colocalization summary")
  summary <- as.data.table(res$summary)
  probability_cols <- intersect(paste0("PP.H",0:4,".abf"),names(summary))
  need(length(probability_cols)==5L, "Unexpected coloc summary schema")
  finite_cols(summary,probability_cols)
  need(all(as.matrix(summary[,..probability_cols])>=0 & as.matrix(summary[,..probability_cols])<=1),
       "Invalid colocalization probability")
  need(all(abs(rowSums(as.matrix(summary[,..probability_cols]))-1)<1e-6), "Coloc posterior probabilities must sum to one")
  fwrite(summary,opath("coloc.tsv"),sep="\t")
  saveRDS(res,opath("coloc.rds"))
  # Retain all H0-H4, not just max H4. Prior sensitivity is a separate named output.
  sensitivity <- list()
  for (p12 in cfg$p12_sensitivity) {
    need(is.finite(p12) && p12>0 && p12<min(cfg$p1,cfg$p2), "Invalid sensitivity prior")
    s <- coloc::coloc.susie(a,b,p1=cfg$p1,p2=cfg$p2,p12=p12)$summary
    s$prior_p12 <- p12
    sensitivity[[length(sensitivity)+1L]] <- as.data.table(s)
  }
  if (length(sensitivity)) fwrite(rbindlist(sensitivity,fill=TRUE),opath("prior_sensitivity.tsv"),sep="\t")
  finish(extra=list(n_signal_comparisons=nrow(summary),interpretation="STATISTICAL_SHARED_SIGNAL_NOT_PROVEN_MECHANISM"))
}
