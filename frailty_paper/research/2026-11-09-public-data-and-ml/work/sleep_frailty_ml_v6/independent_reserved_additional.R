args <- commandArgs(trailingOnly=TRUE)
x <- read.csv(args[1],check.names=FALSE)
stopifnot(!anyDuplicated(x$row_id),all(is.finite(as.matrix(x[,-1]))))
correlations <- lapply(names(x)[-c(1,2)],function(nm) {
 data.frame(model=nm,Spearman=cor(x$effect,x[[nm]],method='spearman'))
})
write.csv(do.call(rbind,correlations),args[2],row.names=FALSE)
x <- read.csv(args[3],check.names=FALSE)
y<-x$label;ranked<-x$ranked==1;pos<-ranked&x$set1==1;neg<-ranked&x$set0==1
covered<-ifelse(y==1,x$set1==1,x$set0==1)
correct<-(pos&y==1)|(neg&y==0)
out<-data.frame(metric=c('coverage','detected_coverage','non_detected_coverage',
 'ranked_fraction','singleton_accuracy','positive_precision','positive_recall','Brier'),
 value=c(mean(covered),mean(covered[y==1]),mean(covered[y==0]),mean(ranked),
 mean(correct[ranked]),mean(y[pos]),sum(pos&y==1)/sum(y==1),mean((x$probability-y)^2)))
write.csv(out,args[4],row.names=FALSE)
