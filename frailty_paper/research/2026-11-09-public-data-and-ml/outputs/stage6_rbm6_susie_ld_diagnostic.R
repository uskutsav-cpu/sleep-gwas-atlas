.libPaths(c('work/stage5_Rlib', .libPaths()))
library(susieR)
d <- read.csv('work/stage6_rbm6_ld_gwas_aligned.csv'); nvars <- nrow(d)
b <- readBin('work/stage6_rbm6_eur_ld_float64.bin','double',n=nvars*nvars,size=8,endian='little'); stopifnot(length(b)==nvars*nvars)
R <- matrix(b,nvars,nvars); rm(b); stopifnot(max(abs(R-t(R)))<1e-6,max(abs(diag(R)-1))<1e-5)
cat('Computing eigen decomposition for', nvars, 'variants\n'); flush.console()
e <- eigen(R,symmetric=TRUE); attr(R,'eigen') <- e
traits <- c('Jansen_insomnia','Mak_FI','FinnGen_HFRS'); ns <- c(313750,402910,500737)
res <- data.frame(trait=traits,working_n=ns,n_variants=nvars,estimated_lambda=NA_real_,reference_n=503,rank_bound=502,min_eigenvalue=min(e$values),package_version=as.character(packageVersion('susieR')))
for(i in seq_along(traits)) { z <- d[[paste0(traits[i], '_beta_alt')]] / d[[paste0(traits[i], '_se')]]; stopifnot(all(is.finite(z))); res$estimated_lambda[i] <- estimate_s_rss(z,R,n=ns[i]); cat(traits[i],res$estimated_lambda[i],'\n');flush.console() }
write.csv(res,'outputs/stage6_rbm6_susie_ld_consistency.csv',row.names=FALSE);print(res)
