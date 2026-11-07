a<-commandArgs(TRUE)[1];d<-list();h<-list()
for(t in c('FI','FFS','HFRS')){
 b<-file.path(a,paste0(t,t,'.sumstats.gz_insomnia.sumstats.gz'));h[[t]]<-lapply(c('hsq1','hsq2','gencov'),function(k)scan(paste0(b,'.',k,'.delete'),quiet=TRUE));stopifnot(all(lengths(h[[t]])==200),all(h[[t]][[1]]>0),all(h[[t]][[2]]>0));d[[t]]<-h[[t]][[3]]/sqrt(h[[t]][[1]]*h[[t]][[2]])
}
stopifnot(identical(h$FI[[2]],h$FFS[[2]]),identical(h$FI[[2]],h$HFRS[[2]]));j<-function(x)sqrt(199/200*sum((x-mean(x))^2));names<-c('HFRS minus FI','HFRS minus FFS','FI minus FFS');se<-c(j(d$HFRS-d$FI),j(d$HFRS-d$FFS),j(d$FI-d$FFS));x<-read.delim(file.path(a,'contrasts.tsv'));stopifnot(identical(x$contrast,names),all(abs(x$new_se_paired-se)<1e-12));cat('independent paired SEs match:',se,'\n')
