args<-commandArgs(trailingOnly=TRUE)
d<-read.csv(args[1],check.names=FALSE)
models<-setdiff(names(d),c('row_id','effect','tile'))
out<-lapply(models,function(m){
  tiles<-split(seq_len(nrow(d)),d$tile)
  values<-sapply(tiles,function(ix){if(length(ix)>=10 && length(unique(d$effect[ix]))>1 && length(unique(d[[m]][ix]))>1) cor(d$effect[ix],d[[m]][ix],method='spearman') else NA_real_})
  data.frame(model=m,Spearman=cor(d$effect,d[[m]],method='spearman'),within_tile_Spearman=mean(values,na.rm=TRUE),tiles=sum(is.finite(values)))
})
options(digits=17);write.csv(do.call(rbind,out),args[2],row.names=FALSE)
