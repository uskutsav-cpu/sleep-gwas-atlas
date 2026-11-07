.libPaths(c('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()));stopifnot(requireNamespace('susieR',quietly=TRUE))
p<-'outputs/shared-genetics-atlas/variant_context/rs62264764_brain_xqtl/ldz_readiness'; old<-'outputs/shared-genetics-atlas/variant_context/rs62264764_dense_qtl/multisignal_pilot';v<-read.delim(file.path(old,'variants.tsv'));R<-as.matrix(read.delim(gzfile(file.path(old,'ld_matrix.tsv.gz')),header=FALSE));dimnames(R)<-list(v$variant,v$variant)
stopifnot(nrow(v)==1548,nrow(R)==1548,max(abs(R-t(R)))<1e-9,max(abs(diag(R)-1))<1e-9)
result<-list();alls<-list()
for (ctx in c('oli','exc')){
 d<-read.delim(paste0('/Volumes/Brain6Work/xqtl-chr3/ldz_readiness/',ctx,'.tsv'));stopifnot(identical(d$variant,v$variant),all(d$num_samples==737),all(d$qtl_se>0));z<-d$qtl_beta/d$qtl_se;names(z)<-d$variant
 cat('starting',ctx,'\n');flush.console(); s<-susieR::estimate_s_rss(z,R,n=737); obj<-suppressWarnings(susieR::kriging_rss(z,R,n=737,s=s));x<-obj$conditional_dist;stopifnot(nrow(x)==nrow(d));x$variant<-d$variant;x$context<-ctx;alls[[ctx]]<-x
 result[[ctx]]<-data.frame(context=ctx,n=737,s=s,max_abs_z_std_diff=max(abs(x$z_std_diff)),count_abs_z_std_diff_gt3=sum(abs(x$z_std_diff)>3),count_abs_z_std_diff_gt5=sum(abs(x$z_std_diff)>5),max_logLR=max(x$logLR));print(result[[ctx]]);flush.console()
}
write.table(do.call(rbind,result),file.path(p,'ldz_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE);write.table(do.call(rbind,alls),file.path(p,'ldz_details.tsv'),sep='\t',quote=FALSE,row.names=FALSE);writeLines(capture.output(sessionInfo()),file.path(p,'runtime.txt'))
