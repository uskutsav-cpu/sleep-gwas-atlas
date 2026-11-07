p<-'outputs/shared-genetics-atlas/frailty_coloc/chr7_three_trait/overlap_sensitivity'
r<-readRDS(file.path(p,'sensitivity_results.rds')); parent<-readRDS('outputs/shared-genetics-atlas/frailty_coloc/chr7_three_trait/pilot_results.rds')
stopifnot(length(r)==8,all(vapply(r,function(z)is.null(z$error),logical(1))))
stopifnot(identical(r[['0']]$results,parent$default$results))
stopifnot(abs(r[['-0.001']]$results$posterior_prob-r[['0']]$results$posterior_prob)<.01,abs(r[['0.001']]$results$posterior_prob-r[['0']]$results$posterior_prob)<.01)
tab<-do.call(rbind,lapply(names(r),function(k)cbind(effective_correlation=as.numeric(k),r[[k]]$results)))
write.table(tab,file.path(p,'clusters.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('Eight fits succeeded; baseline matches parent; near-zero values approach baseline.\n'); print(tab)
