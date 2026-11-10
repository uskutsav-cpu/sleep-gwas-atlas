.libPaths(c('work/stage5_Rlib',.libPaths()));library(susieR);library(coloc)
d<-read.csv('work/meis1_gtex/meis1_pal_sensitivity_aligned.csv');p<-nrow(d);R<-matrix(readBin('work/meis1_gtex/meis1_pal_sensitivity_ld.bin','double',n=p*p,size=8,endian='little'),p,p);dimnames(R)<-list(d$id,d$id)
traits<-c('qtl','Jansen');ns<-c(209,313750);fits<-list();diag<-list();members<-list()
for(i in seq_along(traits)){
 t<-traits[i];z<-d[[paste0(t,'_beta')]]/d[[paste0(t,'_se')]];names(z)<-d$id
 lambda<-tryCatch(estimate_s_rss(z,R,n=ns[i]),error=function(e) NA_real_)
 set.seed(42);f<-susie_rss(z=z,R=R,n=ns[i],L=10,coverage=.95,estimate_residual_variance=FALSE,max_iter=500)
 fits[[t]]<-f;diag[[i]]<-data.frame(trait=t,n=ns[i],variants=p,ld_rank=qr(R,tol=1e-7)$rank,lambda=lambda,converged=f$converged,iterations=f$niter,credible_sets=length(f$sets$cs),max_pip=max(f$pip))
 if(length(f$sets$cs))for(k in seq_along(f$sets$cs)){j<-f$sets$cs[[k]];members[[length(members)+1]]<-data.frame(trait=t,component=f$sets$cs_index[k],variant=d$id[j],rsid=d$rsid[j],pip=f$pip[j])}
 cat(t,'lambda',lambda,'converged',f$converged,'sets',length(f$sets$cs),'\n');flush.console()
}
write.csv(do.call(rbind,diag),'outputs/stage6_meis1_susie_pal_diagnostics.csv',row.names=FALSE)
if(length(members))write.csv(do.call(rbind,members),'outputs/stage6_meis1_susie_pal_members.csv',row.names=FALSE)
saveRDS(fits,'work/meis1_gtex/meis1_susie_pal_fits.rds')
res<-list()
for(t in traits[-1])if(fits$qtl$converged&&fits[[t]]$converged&&length(fits$qtl$sets$cs)&&length(fits[[t]]$sets$cs))for(p12 in c(1e-6,1e-5,5e-5)){
 z<-as.data.frame(coloc.susie(fits$qtl,fits[[t]],p1=1e-4,p2=1e-4,p12=p12)$summary);z$trait<-t;z$p12<-p12;res[[length(res)+1]]<-z
}
if(length(res))write.csv(do.call(rbind,res),'outputs/stage6_meis1_susie_pal_coloc.csv',row.names=FALSE)
