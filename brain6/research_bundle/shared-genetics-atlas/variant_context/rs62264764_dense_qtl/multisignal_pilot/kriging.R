.libPaths(c('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()))
p<-'outputs/shared-genetics-atlas/variant_context/rs62264764_dense_qtl/multisignal_pilot';v<-read.delim(file.path(p,'variants.tsv'));R<-as.matrix(read.delim(gzfile(file.path(p,'ld_matrix.tsv.gz')),header=FALSE));dimnames(R)<-list(v$variant,v$variant)
d<-read.delim(file.path(p,'ld_diagnostics.tsv'));out<-list();sum<-list()
for(name in c('qtl','insomnia','adhd')){
 z<-v[[paste0(name,'_beta')]]/v[[paste0(name,'_se')]];names(z)<-v$variant;params<-d[d$fit==name,];stopifnot(nrow(params)==1)
 obj<-suppressWarnings(susieR::kriging_rss(z,R,n=params$n,s=params$estimated_s));x<-obj$conditional_dist;stopifnot(nrow(x)==nrow(v));x$variant<-v$variant;x$fit<-name;out[[length(out)+1L]]<-x
 sum[[length(sum)+1L]]<-data.frame(fit=name,lambda=params$estimated_s,max_abs_z_std_diff=max(abs(x$z_std_diff)),count_abs_z_std_diff_gt3=sum(abs(x$z_std_diff)>3),count_abs_z_std_diff_gt5=sum(abs(x$z_std_diff)>5),max_logLR=max(x$logLR),stringsAsFactors=FALSE)
}
write.table(do.call(rbind,out),file.path(p,'kriging_diagnostics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,sum),file.path(p,'kriging_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(do.call(rbind,sum))
