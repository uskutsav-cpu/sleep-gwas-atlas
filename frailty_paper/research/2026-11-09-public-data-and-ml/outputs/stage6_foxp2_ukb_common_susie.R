.libPaths(c('work/stage5_Rlib',.libPaths()));library(susieR);library(coloc)
d<-read.csv('work/stage6_foxp2_ukb_common_aligned.csv');p<-nrow(d);stopifnot(p==553);diagout<-list();pairs<-list();k<-0
refs<-c('UKB_EUR_PSD','1000G_EUR503','1000G_FIN99');paths<-c('work/stage6_foxp2_ukb_ld_psd.bin','work/stage6_foxp2_eur503_ukbcommon_ld.bin','work/stage6_foxp2_fin99_ukbcommon_ld.bin');refn<-c(362063,503,99)
for(i in seq_along(refs)){
 R<-matrix(readBin(paths[i],'double',n=p*p,size=8,endian='little'),p,p);dimnames(R)<-list(d$rs,d$rs);e<-eigen(R,symmetric=TRUE);attr(R,'eigen')<-e;fits<-list()
 for(trait in c('MVP','HFRS')){
  z<-if(trait=='MVP')d$mvp_z_anchor else d$hfrs_z_anchor;names(z)<-d$rs;n<-if(trait=='MVP')253769 else 500737;lam<-estimate_s_rss(z,R,n=n);set.seed(42);f<-susie_rss(z=z,R=R,n=n,L=10,estimate_residual_variance=FALSE,max_iter=500,coverage=.95);fits[[trait]]<-f
  diagout[[length(diagout)+1]]<-data.frame(reference=refs[i],trait=trait,n_snps=p,reference_n=refn[i],working_n=n,lambda=lam,min_eigen=min(e$values),converged=f$converged,iterations=f$niter,n_sets=length(f$sets$cs),max_pip=max(f$pip))
 }
 if(length(fits$MVP$sets$cs)>0&&length(fits$HFRS$sets$cs)>0){for(p12 in c(1e-6,1e-5,5e-5)){s<-as.data.frame(coloc.susie(fits$MVP,fits$HFRS,p1=1e-4,p2=1e-4,p12=p12)$summary);s$reference<-refs[i];s$p12<-p12;pairs[[length(pairs)+1]]<-s}}
 saveRDS(fits,paste0('work/stage6_foxp2_',refs[i],'_ukbcommon_fits.rds'));cat(refs[i],'done\n');flush.console()
}
write.csv(do.call(rbind,diagout),'outputs/stage6_foxp2_ukbcommon_susie_diagnostics.csv',row.names=FALSE)
if(length(pairs))write.csv(do.call(rbind,pairs),'outputs/stage6_foxp2_ukbcommon_susie_coloc.csv',row.names=FALSE)
