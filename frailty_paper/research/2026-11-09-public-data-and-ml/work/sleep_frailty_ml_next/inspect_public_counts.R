# Independent base-R inspection; no original author scripts are executed.
p <- 'outputs/sleep_frailty_ai_v9_validation/next_source_screen'
summary_rows <- list()
for (stem in c('single_variants','diplotypes_1_ref_ref_vs_alt_ref','diplotypes_1_ref_ref_vs_ref_alt','diplotypes_1_ref_ref_vs_alt_alt')) {
  x <- readRDS(file.path(p,paste0('attributesData_',stem,'.rds')))
  stopifnot(is.data.frame(x))
  write.table(x,file.path(p,paste0('author_attributes_',stem,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='')
}
for (cell in c('K562','HL60','CHRF','THP1')) {
  x <- readRDS(file.path(p,paste0(cell,'_author_filtered_counts.rds')))
  stopifnot(is.data.frame(x))
  numcols <- names(x)[vapply(x,is.numeric,logical(1))]
  stopifnot(length(numcols)==14, !anyDuplicated(paste(x$Oligo,x$Barcode)))
  values <- as.matrix(x[,numcols])
  stopifnot(all(is.finite(values)),all(values>=0))
  summary_rows[[cell]] <- data.frame(cell=cell,rows=nrow(x),columns=ncol(x),numeric_columns=ncol(values),minimum=min(values),maximum=max(values),zeros=sum(values==0),fraction_noninteger=mean(abs(values-round(values))>1e-8))
  agg <- aggregate(x[,numcols],by=list(Oligo=x$Oligo),FUN=sum)
  nb <- aggregate(x$Barcode,by=list(Oligo=x$Oligo),FUN=function(z)length(unique(z)))
  names(nb)[2] <- 'distinct_barcodes';agg <- merge(agg,nb,by='Oligo')
  write.table(agg,file.path(p,paste0(cell,'_independent_oligo_count_sums.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
}
write.table(do.call(rbind,summary_rows),file.path(p,'independent_R_filtered_count_profile.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('Independent public count and allele-attribute audit PASS\n')
