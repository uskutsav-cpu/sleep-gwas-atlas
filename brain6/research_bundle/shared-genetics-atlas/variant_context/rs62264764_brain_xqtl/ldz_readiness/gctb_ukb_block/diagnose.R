.libPaths(c('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()))
stopifnot(requireNamespace('susieR',quietly=TRUE))
stop('The extracted GCTB matrix failed the pre-specified positive-semidefinite check. Do not run these diagnostics without a separately justified and documented repair protocol.')
p <- 'outputs/shared-genetics-atlas/variant_context/rs62264764_brain_xqtl/ldz_readiness/gctb_ukb_block'
old <- 'outputs/shared-genetics-atlas/variant_context/rs62264764_dense_qtl/multisignal_pilot'
names <- readLines(file.path(p,'strand_safe_variants.txt'))
v <- read.delim(file.path(old,'variants.tsv'))
ids <- match(names,v$variant);stopifnot(!anyNA(ids),!anyDuplicated(ids))
Rnew <- as.matrix(read.delim(gzfile(file.path(p,'strand_safe_ld.tsv.gz')),header=FALSE))
Rold <- as.matrix(read.delim(gzfile(file.path(old,'ld_matrix.tsv.gz')),header=FALSE))[ids,ids]
stopifnot(nrow(Rnew)==length(names),max(abs(diag(Rnew)-1))<1e-4,max(abs(Rnew-t(Rnew)))<1e-5)
sets <- list()
for(ctx in c('oli','exc')) {
 d <- read.delim(paste0('/Volumes/Brain6Work/xqtl-chr3/ldz_readiness/',ctx,'.tsv'))
 ix <- match(names,d$variant);stopifnot(!anyNA(ix),all(d$num_samples[ix]==737))
 sets[[ctx]] <- list(z=d$qtl_beta[ix]/d$qtl_se[ix],n=737)
}
sets$insomnia <- list(z=v$insomnia_beta[ids]/v$insomnia_se[ids],n=313750)
sets$adhd <- list(z=v$adhd_beta[ids]/v$adhd_se[ids],n=128214)
summary <- list();details <- list()
for(ref in c('gctb_ukb','1000g_eur')) for(trait in names(sets)) {
 R <- if(ref=='gctb_ukb') Rnew else Rold
 z <- sets[[trait]]$z;names(z)<-names;n<-sets[[trait]]$n
 cat('starting',ref,trait,'\n');flush.console()
 s <- susieR::estimate_s_rss(z,R,n=n)
 x <- suppressWarnings(susieR::kriging_rss(z,R,n=n,s=s))$conditional_dist
 x$variant<-names;x$reference<-ref;x$trait<-trait
 details[[length(details)+1]]<-x
 summary[[length(summary)+1]]<-data.frame(reference=ref,trait=trait,n=n,markers=length(z),s=s,max_abs_z_std_diff=max(abs(x$z_std_diff)),count_gt3=sum(abs(x$z_std_diff)>3),count_gt5=sum(abs(x$z_std_diff)>5))
 print(summary[[length(summary)]]);flush.console()
}
write.table(do.call(rbind,summary),file.path(p,'ldz_comparison.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,details),file.path(p,'ldz_details.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(p,'runtime.txt'))
