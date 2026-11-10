library(nnet)
d<-read.csv('work/frailpol_gait_three_class.csv',check.names=FALSE)
flags<-read.csv(file('stdin'),colClasses=c('character','logical'))
stopifnot(nrow(d)==668,!anyDuplicated(d$participant_id),!anyDuplicated(flags$participant_id),nrow(flags)==668)
idx<-match(as.character(d$participant_id),flags$participant_id);stopifnot(!anyNA(idx));keep<-flags$eligible[idx]
X<-as.matrix(d[,setdiff(names(d),c('participant_id','Targets'))]);y<-as.integer(d$Targets)
set.seed(20261029);fold<-integer(length(y));for(k in 0:2){idx<-sample(which(y==k));fold[idx]<-rep(1:5,length.out=length(idx))}
fit<-function(traineligible){
 pred<-rep(NA_integer_,length(y));warnings<-0L
 for(i in 1:5){tr<-(fold!=i)&traineligible;te<-(fold==i)&keep
  mu<-colMeans(X[tr,,drop=FALSE]);sd<-apply(X[tr,,drop=FALSE],2,sd);sd[sd==0]<-1;z<-sweep(sweep(X,2,mu,'-'),2,sd,'/')
  train<-as.data.frame(z[tr,,drop=FALSE]);names(train)<-paste0('x',seq_len(ncol(train)));train$label<-factor(y[tr],levels=0:2)
  mod<-withCallingHandlers(multinom(label~.,data=train,decay=0.1,maxit=1000,trace=FALSE),warning=function(w){warnings<<-warnings+1L;invokeRestart('muffleWarning')})
  test<-as.data.frame(z[te,,drop=FALSE]);names(test)<-paste0('x',seq_len(ncol(test)));pred[te]<-as.integer(as.character(predict(mod,test)))
  if(mod$convergence!=0)warnings<-warnings+1L
 }
 cm<-table(observed=factor(y[keep],levels=0:2),predicted=factor(pred[keep],levels=0:2));list(cm=cm,pred=pred,warnings=warnings)
}
a<-fit(rep(TRUE,length(y)));b<-fit(keep)
metrics<-function(label,x){cm<-x$cm;rec<-diag(cm)/rowSums(cm);data.frame(model=label,n=sum(cm),accuracy=sum(diag(cm))/sum(cm),balanced_accuracy=mean(rec),recall_robust=rec[1],recall_prefrail=rec[2],recall_frail=rec[3],warnings=x$warnings)}
r<-rbind(metrics('all_training_common_evaluation',a),metrics('exclude_discrepant_or_incomplete_training',b));r$prediction_changes<-sum(a$pred[keep]!=b$pred[keep]);write.csv(r,'outputs/stage6_frailpol_label_sensitivity.csv',row.names=FALSE)
write.csv(as.data.frame.matrix(a$cm),'outputs/stage6_frailpol_label_sensitivity_original_cm.csv');write.csv(as.data.frame.matrix(b$cm),'outputs/stage6_frailpol_label_sensitivity_excluded_cm.csv');print(r);print(a$cm);print(b$cm)
