.libPaths(c('work/stage5_Rlib',.libPaths()));library(susieR)
d<-read.csv('work/stage6_rbm6_ld_gwas_aligned.csv');p<-nrow(d); R<-matrix(readBin('work/stage6_rbm6_eur_ld_float64.bin','double',n=p*p,size=8,endian='little'),p,p);stopifnot(all(is.finite(R)))
traits<-c('Jansen_insomnia','FinnGen_HFRS'); ns<-c(313750,500737); outputs<-list()
for(i in seq_along(traits)){
 name<-traits[i];z<-d[[paste0(name,'_beta_alt')]]/d[[paste0(name,'_se')]];cat('Starting',name,'\n');flush.console();set.seed(42)
 f<-susie_rss(z=z,R=R,n=ns[i],L=10,estimate_residual_variance=FALSE,max_iter=500,coverage=.95,verbose=FALSE)
 cs<-f$sets$cs; sums<-data.frame(component=integer(),size=integer(),lead=character(),lead_pip=numeric(),min_abs_corr=numeric())
 if(length(cs)>0){for(k in seq_along(cs)){ids<-cs[[k]]; lead<-ids[which.max(f$pip[ids])]; pur<-f$sets$purity[k,'min.abs.corr'];sums<-rbind(sums,data.frame(component=f$sets$cs_index[k],size=length(ids),lead=d$rsid[lead],lead_pip=f$pip[lead],min_abs_corr=pur))}}
 write.csv(sums,paste0('outputs/stage6_rbm6_',name,'_susie_candidate_sets.csv'),row.names=FALSE)
 saveRDS(f,paste0('work/stage6_rbm6_',name,'_susie.rds')); outputs[[i]]<-data.frame(trait=name,converged=f$converged,iterations=f$niter,n_sets=length(cs),maximum_pip=max(f$pip),n_variants=p)
 cat('Finished',name,'converged',f$converged,'iterations',f$niter,'credible sets',length(cs),'\n');flush.console()
}
write.csv(do.call(rbind,outputs),'outputs/stage6_rbm6_susie_attempt_status.csv',row.names=FALSE)
