.libPaths(c('work/stage5_Rlib', .libPaths()))
library(susieR)
d<-read.csv('work/stage5_susie_diagnostic/mak_fi_qtl_summary.csv',stringsAsFactors=FALSE)
n<-nrow(d);b<-readBin('work/stage5_susie_diagnostic/eur_ld_float64.bin','double',n=n*n,size=8,endian='little');stopifnot(length(b)==n*n);R<-matrix(b,n,n);rm(b);stopifnot(max(abs(R-t(R)))<1e-6, max(abs(diag(R)-1))<1e-6)
qz<-d$qtl_beta/d$qtl_se;gz<-d$gwas_beta/d$gwas_se;Nq<-702;Ng<-median(d$gwas_n,na.rm=TRUE);stopifnot(all(is.finite(qz)),all(is.finite(gz)),is.finite(Ng))
cat('Computing eigen decomposition for',n,'variants\n');flush.console();e<-eigen(R,symmetric=TRUE);attr(R,'eigen')<-e
sq<-estimate_s_rss(qz,R,n=Nq);sg<-estimate_s_rss(gz,R,n=Ng)
result<-data.frame(trait=c('SMAD5_muscle_splice','Mak_FI'),n=c(Nq,Ng),n_variants=n,estimated_lambda=c(sq,sg),reference_n=503,min_eigenvalue=min(e$values),max_eigenvalue=max(e$values),reference_rank_bound=502,package_version=as.character(packageVersion('susieR')))
write.csv(result,'outputs/stage5_smad5_susie_ld_consistency.csv',row.names=FALSE);print(result)
