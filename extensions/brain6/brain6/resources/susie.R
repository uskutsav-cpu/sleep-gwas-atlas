#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
need_package("susieR"); need_package("coloc")
d <- read_dt(cfg$locus_file)
R <- as.matrix(read.table(cfg$ld_file,header=FALSE,sep="\t"))
need(!anyDuplicated(d$SNP) && nrow(R)==nrow(d) && ncol(R)==nrow(d), "LD order/shape mismatch")
need(all(is.finite(R)) && max(abs(R-t(R)))<1e-6 && max(abs(diag(R)-1))<1e-6,
     "LD is not a finite symmetric correlation matrix")
need(min(eigen(R,symmetric=TRUE,only.values=TRUE)$values)>=-1e-6, "LD is not PSD")
dimnames(R) <- list(d$SNP,d$SNP)
suffix <- as.character(cfg$trait_number)
need(suffix %in% c("1","2"), "trait_number must be 1 or 2")
b <- d[[paste0("BETA",suffix)]]; s <- d[[paste0("SE",suffix)]]
n <- d[[paste0("N",suffix)]]
need(all(is.finite(b)) && all(is.finite(s)) && all(s>0), "Invalid summary effects")
need(is.numeric(cfg$sample_size) && cfg$sample_size>0 && !is.null(cfg$sample_size_justification),
     "Provide a reviewed locus-level N; do not silently take a median")
need(all(is.finite(n)) && max(abs(n/cfg$sample_size-1))<=cfg$max_relative_n_deviation,
     "Per-SNP N variation exceeds reviewed locus-N approximation")
need(cfg$trait_type %in% c("cc","quant"), "Unknown trait type")
dataset <- list(snp=d$SNP,position=d$BP,beta=b,varbeta=s*s,N=cfg$sample_size,
                type=cfg$trait_type,LD=R)
if (cfg$trait_type=="cc") {
  need(is.numeric(cfg$case_fraction) && cfg$case_fraction>0 && cfg$case_fraction<1,
       "Case fraction required for case-control locus")
  dataset$s <- cfg$case_fraction
} else {
  need(is.numeric(cfg$sdY) && cfg$sdY>0, "Quantitative SD requires sourced value; do not invent sdY=1")
  dataset$sdY <- cfg$sdY
}
# Diagnose inconsistency. Do not flip z-scores or shrink LD to force a desired fit.
lambda <- susieR::estimate_s_rss(z=b/s,R=R,n=cfg$sample_size)
need(is.finite(lambda) && lambda<=cfg$max_ld_mismatch, "FAILED_QC: excessive z/LD mismatch")
fit <- coloc::runsusie(dataset,L=cfg$L,coverage=cfg$coverage,min_abs_corr=cfg$min_abs_corr,
                      estimate_residual_variance=FALSE,maxit=cfg$max_iterations,
                      repeat_until_convergence=FALSE)
need(isTRUE(fit$converged), "SuSiE did not converge")
saveRDS(fit,opath("fit.rds"))
fwrite(data.table(SNP=d$SNP,CHR=d$CHR,BP=d$BP,PIP=fit$pip),opath("pip.tsv"),sep="\t")
cs <- fit$sets$cs
cs_rows <- data.table(credible_set=character(),SNP=character(),PIP=numeric())
if (length(cs)) {
  cs_rows <- rbindlist(lapply(names(cs),function(name) {
    i<-cs[[name]];data.table(credible_set=name,SNP=d$SNP[i],PIP=fit$pip[i])
  }))
}
fwrite(cs_rows,opath("credible_sets.tsv"),sep="\t")
write_json(list(lambda=lambda,sets=fit$sets,n=nrow(d),settings=cfg),opath("diagnostics.json"),
           auto_unbox=TRUE,pretty=TRUE,na="null",digits=NA)
finish(if(length(cs)) "PASS" else "NO_SIGNAL",list(n_snps=nrow(d),n_credible_sets=length(cs),ld_mismatch=lambda))
