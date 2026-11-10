args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)) args[1] else "outputs/sleep_frailty_ai_v1"
x<-read.csv(file.path(root,"test_predictions.csv"),check.names=FALSE)
stopifnot(nrow(x)==3542,!anyDuplicated(x$row_index),all(x$chromosome %in% 19:22),sum(x$source_detected_call)==100)
y<-x$source_detected_call
models<-c("constant","kmer_logistic","frozen_hyenadna","finetuned_hyenadna","small_dna_cnn")
out<-lapply(models,function(name){
 p<-x[[name]];stopifnot(all(is.finite(p)),all(p>=0 & p<=1))
 # Mann–Whitney rank formula independently checks trapezoidal ROC calculation.
 auc<-(sum(rank(p,ties.method="average")[y==1])-sum(y)*(sum(y)+1)/2)/(sum(y)*sum(y==0))
 o<-order(-p);ss<-p[o];yy<-y[o];end<-c(which(diff(ss)!=0),length(ss));tp<-cumsum(yy)[end]
 ap<-sum(diff(c(0,tp/sum(y)))*tp/end)
 data.frame(model=name,n=length(y),detected=sum(y),auroc=auc,average_precision=ap)
})
out<-do.call(rbind,out)
write.csv(out,file.path(root,"independent_R_metrics.csv"),row.names=FALSE)
print(out,row.names=FALSE)
