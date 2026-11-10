.libPaths(c('work/stage5_Rlib',.libPaths()));library(susieR);library(coloc)
d<-read.csv('work/meis1_gtex/meis1_hfrs_aligned.csv');p<-nrow(d);fits<-list();diags<-list();members<-list();traits<-c('qtl','HFRS');ns<-c(209,500737);pops<-c('nfe','fin')
for(i in 1:2){
 t<-traits[i];R<-matrix(readBin(paste0('work/meis1_gtex/meis1_hfrs_',pops[i],'_ld.bin'),'double',n=p*p,size=8,endian='little'),p,p);dimnames(R)<-list(d$id,d$id)
 z<-d[[paste0(t,'_beta')]]/d[[paste0(t,'_se')]];names(z)<-d$id;lambda<-tryCatch(estimate_s_rss(z,R,n=ns[i]),error=function(e) NA_real_)
 set.seed(42);f<-susie_rss(z=z,R=R,n=ns[i],L=10,coverage=.95,estimate_residual_variance=FALSE,max_iter=500);fits[[t]]<-f
 diags[[i]]<-data.frame(trait=t,n=ns[i],variants=p,reference=pops[i],ld_rank=qr(R,tol=1e-7)$rank,lambda=lambda,converged=f$converged,iterations=f$niter,credible_sets=length(f$sets$cs),max_pip=max(f$pip));cat(t,'lambda',lambda,'converged',f$converged,'sets',length(f$sets$cs),'\n');flush.console()
 if(length(f$sets$cs))for(k in seq_along(f$sets$cs)){j<-f$sets$cs[[k]];members[[length(members)+1]]<-data.frame(trait=t,component=f$sets$cs_index[k],variant=d$id[j],rsid=d$rsid[j],pip=f$pip[j])}
}
write.csv(do.call(rbind,diags),'outputs/stage6_meis1_hfrs_susie_diagnostics.csv',row.names=FALSE);if(length(members))write.csv(do.call(rbind,members),'outputs/stage6_meis1_hfrs_susie_members.csv',row.names=FALSE);saveRDS(fits,'work/meis1_gtex/meis1_hfrs_susie_fits.rds')
res<-list();if(all(sapply(fits,function(f) f$converged&&length(f$sets$cs)>0)))for(p12 in c(1e-6,1e-5,5e-5)){z<-as.data.frame(coloc.susie(fits$qtl,fits$HFRS,p1=1e-4,p2=1e-4,p12=p12)$summary);z$p12<-p12;res[[length(res)+1]]<-z};if(length(res))write.csv(do.call(rbind,res),'outputs/stage6_meis1_hfrs_susie_coloc.csv',row.names=FALSE)
