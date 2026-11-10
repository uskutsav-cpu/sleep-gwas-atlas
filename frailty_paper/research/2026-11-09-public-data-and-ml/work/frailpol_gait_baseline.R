library(nnet)
d<-read.csv('work/frailpol_gait_three_class.csv',check.names=FALSE)
stopifnot(nrow(d)==668,!anyDuplicated(d$participant_id),!anyNA(d),all(d$Targets %in% 0:2))
X<-as.matrix(d[,setdiff(names(d),c('participant_id','Targets'))]); y<-as.integer(d$Targets); stopifnot(all(is.finite(X)))
set.seed(20261029); fold<-integer(length(y));for(k in 0:2){idx<-sample(which(y==k));fold[idx]<-rep(1:5,length.out=length(idx))}
pred<-rep(NA_integer_,length(y)); warnings<-character()
for (i in 1:5){tr<-fold!=i;te<-!tr;mu<-colMeans(X[tr,,drop=FALSE]);sd<-apply(X[tr,,drop=FALSE],2,sd);sd[sd==0]<-1
 z<-sweep(sweep(X,2,mu,'-'),2,sd,'/');train<-as.data.frame(z[tr,,drop=FALSE]); names(train)<-paste0('x',seq_len(ncol(train)));train$label<-factor(y[tr],levels=0:2)
 mod<-withCallingHandlers(multinom(label~.,data=train,decay=0.1,maxit=1000,trace=FALSE),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 test<-as.data.frame(z[te,,drop=FALSE]);names(test)<-paste0('x',seq_len(ncol(test)));pred[te]<-as.integer(as.character(predict(mod,test)))
 if(mod$convergence!=0)warnings<-c(warnings,paste('nonconvergence fold',i))
}
cm<-table(observed=factor(y,levels=0:2),predicted=factor(pred,levels=0:2));recall<-diag(cm)/rowSums(cm);acc<-sum(diag(cm))/sum(cm);baseline<-max(table(y))/length(y)
result<-data.frame(metric=c('N','accuracy','balanced_accuracy','majority_accuracy','recall_robust','recall_prefrail','recall_frail','warning_count'),value=c(length(y),acc,mean(recall),baseline,recall,length(warnings)))
write.csv(result,'outputs/stage6_frailpol_gait_baseline_metrics.csv',row.names=FALSE);write.csv(as.data.frame.matrix(cm),'outputs/stage6_frailpol_gait_baseline_confusion.csv');writeLines(warnings,'work/frailpol_baseline_warnings.txt');print(result);print(cm)
