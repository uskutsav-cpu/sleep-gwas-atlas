.libPaths(c('work/stage5_Rlib', .libPaths()))
library(susieR)
d <- read.csv('outputs/stage5_smad5_susie_ld_diagnostic_inputs.csv', stringsAsFactors=FALSE)
n <- nrow(d)
b <- readBin('work/stage5_susie_diagnostic/eur_ld_float64.bin', 'double', n=n*n, size=8, endian='little')
stopifnot(length(b)==n*n)
R <- matrix(b,n,n); rm(b)
stopifnot(max(abs(R-t(R)))<1e-6,max(abs(diag(R)-1))<1e-6)
attr(R,'eigen') <- eigen(R,symmetric=TRUE)
lambdas <- read.csv('outputs/stage5_smad5_susie_ld_consistency.csv')
for (i in seq_len(nrow(lambdas))) {
  trait <- lambdas$trait[i]
  z <- if (trait=='SMAD5_muscle_splice') d$qtl_beta/d$qtl_se else d$gwas_beta/d$gwas_se
  k <- kriging_rss(z,R,n=lambdas$n[i],s=lambdas$estimated_lambda[i])
  out <- data.frame(trait=trait,variant=d$variant,original_z=z,k$conditional_dist)
  out$possible_flip_flag <- out$logLR > 2 & abs(out$z) > 2
  out <- out[order(-out$logLR),]
  write.csv(out,paste0('outputs/stage5_smad5_',trait,'_variant_ld_diagnostic.csv'),row.names=FALSE)
  cat(trait,'flagged:',sum(out$possible_flip_flag),'max logLR:',max(out$logLR),'max abs standardized residual:',max(abs(out$z_std_diff)),'\n')
}
