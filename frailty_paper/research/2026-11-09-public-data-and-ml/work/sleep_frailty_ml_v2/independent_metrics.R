args <- commandArgs(trailingOnly=TRUE)
x <- read.csv(args[1],check.names=FALSE)
stopifnot(!anyDuplicated(x$row_id),all(x$label %in% c(0,1)),all(is.finite(as.matrix(x[,-c(1,2)]))))
out <- lapply(names(x)[-c(1,2)],function(name) {
 y<-x$label;p<-x[[name]];n<-sum(y);m<-sum(y==0)
 auc<-(sum(rank(p,ties.method="average")[y==1])-n*(n+1)/2)/(n*m)
 o<-order(-p);ss<-p[o];yy<-y[o];ends<-c(which(diff(ss)!=0),length(ss));tp<-cumsum(yy)[ends]
 ap<-sum(diff(c(0,tp/n))*tp/ends)
 data.frame(model=name,n=length(y),positive=n,auroc=auc,average_precision=ap)
})
write.csv(do.call(rbind,out),args[2],row.names=FALSE)
