.libPaths(c('work/stage5_Rlib',.libPaths()));library(susieR);library(coloc)
ind<-read.csv('outputs/stage6_foxp2_gnomad_narrow_variants.csv');all<-read.csv('work/stage6_foxp2_susie_aligned.csv');d<-all[match(ind$rs,all$rs),];p<-nrow(d);stopifnot(p==1257,all(d$rs==ind$rs));g<-matrix(readBin('work/stage6_foxp2_gnomad_narrow_1000g_genotypes.bin','double',n=p*503,size=8,endian='little'),503,p);fi<-scan('work/stage6_foxp2_fin_reference_indices.txt',quiet=TRUE);af<-colMeans(g[fi,,drop=FALSE])/2;keep<-af>=.01 & af<=.99;g<-g[,keep,drop=FALSE];d<-d[keep,];p<-nrow(d);out<-list();diagout<-list()
for(pop in c('EUR503','FIN99')){
 mat<-if(pop=='EUR503')g else g[fi,,drop=FALSE];stopifnot(all(apply(mat,2,sd)>0));R<-cor(mat);dimnames(R)<-list(d$rs,d$rs);e<-eigen(R,symmetric=TRUE);attr(R,'eigen')<-e;fits<-list()
 for(tr in c('MVP','HFRS')){
  z<-if(tr=='MVP')d$mvp_z_anchor else d$hfrs_z_anchor;names(z)<-d$rs;n<-if(tr=='MVP')253769 else 500737;lam<-estimate_s_rss(z,R,n=n);set.seed(42);f<-susie_rss(z=z,R=R,n=n,L=10,coverage=.95,estimate_residual_variance=FALSE,max_iter=500);fits[[tr]]<-f;diagout[[length(diagout)+1]]<-data.frame(reference=pop,trait=tr,n=p,lambda=lam,converged=f$converged,n_sets=length(f$sets$cs),max_pip=max(f$pip))
 }
 if(length(fits$MVP$sets$cs)>0&&length(fits$HFRS$sets$cs)>0)for(p12 in c(1e-6,1e-5,5e-5)){s<-as.data.frame(coloc.susie(fits$MVP,fits$HFRS,p1=1e-4,p2=1e-4,p12=p12)$summary);s$reference<-pop;s$p12<-p12;out[[length(out)+1]]<-s}
 cat(pop,'finished\n');flush.console()
}
write.csv(do.call(rbind,diagout),'outputs/stage6_foxp2_gnomad_narrow_1000g_diagnostics.csv',row.names=FALSE);if(length(out))write.csv(do.call(rbind,out),'outputs/stage6_foxp2_gnomad_narrow_1000g_coloc.csv',row.names=FALSE)
