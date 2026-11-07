.libPaths(c('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()))
stopifnot(requireNamespace('susieR',quietly=TRUE))
p <- 'outputs/shared-genetics-atlas/variant_context/rs62264764_brain_xqtl/ldz_readiness/gctb_ukb_block'
old <- 'outputs/shared-genetics-atlas/variant_context/rs62264764_dense_qtl/multisignal_pilot'
variants<-readLines(file.path(p,'strand_safe_variants.txt')); v<-read.delim(file.path(old,'variants.tsv'))
ids<-match(variants,v$variant);stopifnot(!anyNA(ids),!anyDuplicated(ids))
Rold<-as.matrix(read.delim(gzfile(file.path(old,'ld_matrix.tsv.gz')),header=FALSE))[ids,ids]
z<-v$adhd_beta[ids]/v$adhd_se[ids];names(z)<-variants
summary<-list();details<-list()
for(ref in c('eigenclip','ridge_0p001','1000g_eur')) {
 R<-if(ref=='1000g_eur') Rold else as.matrix(read.delim(gzfile(file.path(p,paste0(ref,'.tsv.gz'))),header=FALSE))
 cat('starting',ref,'\n');flush.console()
 ans<-tryCatch({
  s<-susieR::estimate_s_rss(z,R,n=64699)
  x<-suppressWarnings(susieR::kriging_rss(z,R,n=64699,s=s))$conditional_dist
  x$variant<-variants;x$reference<-ref
  list(summary=data.frame(reference=ref,n=64699,markers=length(z),s=s,max_abs_z_std_diff=max(abs(x$z_std_diff)),count_gt3=sum(abs(x$z_std_diff)>3),count_gt5=sum(abs(x$z_std_diff)>5),error=''),detail=x)
 },error=function(e)list(summary=data.frame(reference=ref,n=64699,markers=length(z),s=NA,max_abs_z_std_diff=NA,count_gt3=NA,count_gt5=NA,error=conditionMessage(e)),detail=NULL))
 summary[[length(summary)+1]]<-ans$summary;if(!is.null(ans$detail))details[[length(details)+1]]<-ans$detail
 write.table(do.call(rbind,summary),file.path(p,'adhd_min_n_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 if(length(details))write.table(do.call(rbind,details),file.path(p,'adhd_min_n_details.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 print(ans$summary);flush.console()
}
writeLines(capture.output(sessionInfo()),file.path(p,'adhd_min_n_runtime.txt'))
