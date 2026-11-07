p <- 'outputs/shared-genetics-atlas/variant_context/chr5_dense_qtl/hyprcoloc_pilot'
smoke<-readRDS(file.path(p,'smoke_test.rds')); stopifnot(nrow(smoke$results)==3)
stopifnot(identical(as.character(smoke$results$traits),c('T1, T2, T3, T4, T5','T6, T7, T8','T9, T10')))
stopifnot(identical(as.character(smoke$results$candidate_snp),c('rs11591147','rs12117612','rs7524677')))
r<-readRDS(file.path(p,'pilot_results.rds'));s<-read.delim(file.path(p,'settings.tsv'));stopifnot(identical(names(r),s$name))
rows<-list()
for(i in seq_along(r)){d<-r[[i]]$results;stopifnot(nrow(d)>0,all(d$posterior_prob>=0),all(d$posterior_prob<=1),all(d$regional_prob>=0),all(d$regional_prob<=1)); d$setting<-names(r)[i];rows[[i]]<-d}
res<-do.call(rbind,rows);write.table(res,file.path(p,'clusters.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('Documented example cluster names and lead markers match; 4 planned settings returned valid probability ranges. These checks do not calibrate the model.\n')
