a<-commandArgs(TRUE)[1]; stopifnot(!is.na(a)); d<-do.call(rbind,lapply(file.path(a,c('runs.tsv','resume_runs.tsv','resume2_runs.tsv')),read.delim)); stopifnot(nrow(d)==80,!anyDuplicated(d[c('seed','rho','method')]),setequal(d$seed,71001:71010),setequal(d$rho,c(0,.3)),setequal(d$method,c('smooth','unsmoothed','maxp_bonferroni','maxp_BY')))
old<-read.delim(file.path(a,'complete_summary.tsv')); out<-list(); pair<-list()
for(m in unique(d$method)){
 for(r in c(0,.3)){
 x<-subset(d,method==m & rho==r);stopifnot(nrow(x)==10,all(abs(x$fdp_markers-x$false_markers/pmax(x$markers,1))<1e-14))
 avg<-mean(x$fdp_markers);expected<-subset(old,method==m & rho==r)$mean_marker_fdp;stopifnot(abs(avg-expected)<1e-12)
 se<-sd(x$fdp_markers)/sqrt(10);ci<-if(se==0) c(NA_real_,NA_real_) else avg+c(-1,1)*qt(.975,9)*se;b<-binom.test(sum(x$false_markers>0),10)$conf.int
 out[[length(out)+1]]<-data.frame(method=m,rho=r,mean_fdp=avg,mcse=se,t_lower=ci[1],t_upper=ci[2],runs_with_false=sum(x$false_markers>0),run_false_lower=b[1],run_false_upper=b[2])
 }
 a0<-subset(d,method==m & rho==0); a3<-subset(d,method==m & rho==.3);a0<-a0[order(a0$seed),];a3<-a3[order(a3$seed),];stopifnot(identical(a0$seed,a3$seed));diff<-a3$fdp_markers-a0$fdp_markers;se<-sd(diff)/sqrt(10);ci<-if(se==0) c(NA_real_,NA_real_) else mean(diff)+c(-1,1)*qt(.975,9)*se
 pair[[length(pair)+1]]<-data.frame(method=m,mean_difference=mean(diff),paired_mcse=se,t_lower=ci[1],t_upper=ci[2])
}
write.table(do.call(rbind,out),file.path(a,'uncertainty','per_setting.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,pair),file.path(a,'uncertainty','paired.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
print(do.call(rbind,out));print(do.call(rbind,pair))
