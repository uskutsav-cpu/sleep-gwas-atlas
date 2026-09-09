#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
need_package("GenomicSEM")
if (cfg$mode=="covariance") {
  need(length(cfg$traits)==length(cfg$trait_names) && !anyDuplicated(cfg$trait_names), "Trait naming mismatch")
  need(length(cfg$sample_prevalence)==length(cfg$traits) && length(cfg$population_prevalence)==length(cfg$traits),
       "Prevalence vector length mismatch; use JSON null for continuous traits")
  cov <- GenomicSEM::ldsc(traits=cfg$traits,sample.prev=cfg$sample_prevalence,
         population.prev=cfg$population_prevalence,ld=cfg$ld_directory,wld=cfg$weights_directory,
         trait.names=cfg$trait_names,n.blocks=cfg$n_blocks,ldsc.log=opath("ldsc.log"))
  need(all(c("S","V","I") %in% names(cov)), "Missing genuine multivariate S/V/I matrices")
  k <- length(cfg$traits); q <- k*(k+1)/2
  need(all(dim(cov$S)==c(k,k)) && all(dim(cov$V)==c(q,q)), "Covariance dimensions incorrect")
  for (name in c("S","V","I")) {
    need(all(is.finite(cov[[name]])), paste("Nonfinite",name))
    write.table(cov[[name]],opath(paste0(name,".tsv")),sep="\t",quote=FALSE,col.names=NA)
  }
  need(all(diag(cov$S)>0), "Nonpositive genetic variance")
  need(min(eigen(cov$S,symmetric=TRUE,only.values=TRUE)$values)>cfg$min_eigenvalue,
       "Genetic covariance is not usable; do not repair automatically")
  need(min(eigen(cov$V,symmetric=TRUE,only.values=TRUE)$values)>cfg$min_eigenvalue,
       "Sampling covariance is not usable")
  saveRDS(cov,opath("covariance.rds"))
  finish(extra=list(n_traits=k,n_sampling_parameters=q))
} else if (cfg$mode=="model") {
  need(isTRUE(cfg$model_reviewed), "Model specification requires review")
  cov <- readRDS(cfg$covariance_file)
  result <- GenomicSEM::usermodel(covstruc=cov,estimation="DWLS",model=cfg$model,
              std.lv=TRUE,imp_cov=TRUE,fix_resid=FALSE,toler=FALSE)
  saveRDS(result,opath("model.rds"))
  capture.output(str(result),file=opath("model_summary.txt"))
  # Fit interpretation is a scientific review, not inferred from exit code.
  finish("INSUFFICIENT_EVIDENCE",list(reason="MODEL_FIT_REQUIRES_SCIENTIFIC_REVIEW_NO_AUTOMATIC_CAUSAL_CLAIM"))
} else stop("Unknown GenomicSEM mode")
