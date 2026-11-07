.libPaths(c('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()))
stopifnot(requireNamespace('susieR',quietly=TRUE))
p <- 'outputs/shared-genetics-atlas/variant_context/rs62264764_brain_xqtl/ldz_readiness/gctb_ukb_block'
old <- 'outputs/shared-genetics-atlas/variant_context/rs62264764_dense_qtl/multisignal_pilot'
stopifnot(file.exists(file.path(p,'repair_sensitivity_receipt.json')))
variants <- readLines(file.path(p,'strand_safe_variants.txt'))
v <- read.delim(file.path(old,'variants.tsv'))
ids <- match(variants,v$variant);stopifnot(!anyNA(ids),!anyDuplicated(ids))
Rold <- as.matrix(read.delim(gzfile(file.path(old,'ld_matrix.tsv.gz')),header=FALSE))[ids,ids]
sets <- list()
for(ctx in c('oli','exc')) {
 d <- read.delim(paste0('/Volumes/Brain6Work/xqtl-chr3/ldz_readiness/',ctx,'.tsv'))
 ix <- match(variants,d$variant);stopifnot(!anyNA(ix),all(d$num_samples[ix]==737))
 sets[[ctx]] <- list(z=d$qtl_beta[ix]/d$qtl_se[ix],n=737)
}
sets$insomnia <- list(z=v$insomnia_beta[ids]/v$insomnia_se[ids],n=313750)
sets$adhd <- list(z=v$adhd_beta[ids]/v$adhd_se[ids],n=128214)
summaries <- list(); details <- list()
for(ref in c('eigenclip','ridge_0p001','1000g_eur')) {
 R <- if(ref=='1000g_eur') Rold else as.matrix(read.delim(gzfile(file.path(p,paste0(ref,'.tsv.gz'))),header=FALSE))
 stopifnot(nrow(R)==length(variants),max(abs(diag(R)-1))<1e-6,max(abs(R-t(R)))<1e-6)
 for(trait in names(sets)) {
  z <- sets[[trait]]$z;names(z)<-variants;n<-sets[[trait]]$n
  cat('starting',ref,trait,'\n');flush.console()
  answer<-tryCatch({
   s <- susieR::estimate_s_rss(z,R,n=n)
   x <- suppressWarnings(susieR::kriging_rss(z,R,n=n,s=s))$conditional_dist
   x$variant<-variants;x$reference<-ref;x$trait<-trait
   list(summary=data.frame(reference=ref,trait=trait,n=n,markers=length(z),s=s,max_abs_z_std_diff=max(abs(x$z_std_diff)),count_gt3=sum(abs(x$z_std_diff)>3),count_gt5=sum(abs(x$z_std_diff)>5),error=''),detail=x)
  },error=function(e) list(summary=data.frame(reference=ref,trait=trait,n=n,markers=length(z),s=NA,max_abs_z_std_diff=NA,count_gt3=NA,count_gt5=NA,error=conditionMessage(e)),detail=NULL))
  summaries[[length(summaries)+1]]<-answer$summary
  if(!is.null(answer$detail)) details[[length(details)+1]]<-answer$detail
  write.table(do.call(rbind,summaries),file.path(p,'repair_ldz_comparison.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  if(length(details)) write.table(do.call(rbind,details),file.path(p,'repair_ldz_details.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  print(answer$summary);flush.console()
 }
}
writeLines(capture.output(sessionInfo()),file.path(p,'repair_runtime.txt'))
