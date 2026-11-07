#!/usr/bin/env Rscript
root<-'/Volumes/Brain6NHANES/source'
out<-'outputs/shared-genetics-atlas/component_exclusion/nhanes_public/item_mapping'
paths<-sort(list.files(root,pattern='\\.XPT$',full.names=TRUE));stopifnot(length(paths)==24)
files<-list();variables<-list()
for(path in paths){
 d<-foreign::read.xport(path);stopifnot('SEQN' %in% names(d),!anyNA(d$SEQN))
 name<-basename(path);files[[length(files)+1]]<-data.frame(file=name,rows=nrow(d),
  unique_participants=length(unique(d$SEQN)),duplicated_id_rows=sum(duplicated(d$SEQN)),columns=ncol(d))
 variables[[length(variables)+1]]<-data.frame(file=name,variable=names(d),
  nonmissing=vapply(d,function(x)sum(!is.na(x)),integer(1)),stringsAsFactors=FALSE)
}
write.table(do.call(rbind,files),file.path(out,'file_inventory.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,variables),file.path(out,'variable_inventory.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('PASS: parsed 24 public transport files; participant rows not exported.\n')
