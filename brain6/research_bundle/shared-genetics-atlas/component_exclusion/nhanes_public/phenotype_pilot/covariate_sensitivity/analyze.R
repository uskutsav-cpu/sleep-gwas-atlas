#!/usr/bin/env Rscript
source('outputs/shared-genetics-atlas/component_exclusion/nhanes_public/phenotype_pilot/score_and_analyze.R')
out <- 'outputs/shared-genetics-atlas/component_exclusion/nhanes_public/phenotype_pilot/covariate_sensitivity'
a <- do.call(rbind,data_all[grep(' primary$',names(data_all))]);stopifnot(!anyDuplicated(a$SEQN))
demo <- do.call(rbind,lapply(c('G','H'),function(s)readfile('DEMO',s)[,c('SEQN','RIAGENDR')]));stopifnot(!anyDuplicated(demo$SEQN))
a$female <- code(demo$RIAGENDR[match(a$SEQN,demo$SEQN)],c('1'=0,'2'=1))
a$eligible <- !is.na(a$age)&a$age>=60 & complete.cases(a[,c('fi','fi48','fi47','ffs','ffs4','sleep','fatigue')])
vars<-c('fi','ffs','fi43','ffs2','age','female');stopifnot(sum(a$eligible)==2433,all(complete.cases(a[a$eligible,vars])))
design<-svydesign(id=~SDMVPSU,strata=~STRATUM,weights=~weight,nest=TRUE,data=a);sub<-subset(design,eligible)
for(i in seq_along(vars))for(j in i:length(vars)){name<-paste0('m',i,'_',j);sub$variables[[name]]<-sub$variables[[vars[i]]]*sub$variables[[vars[j]]]} 
prods<-unlist(lapply(seq_along(vars),function(i)paste0('m',i,'_',i:length(vars))))
m<-svymean(as.formula(paste('~',paste(c(vars,prods),collapse='+'))),sub);v<-as.numeric(m);V<-vcov(m)
covariance<-function(t){u<-t[1:6];s<-matrix(0,6,6);k<-7L;for(i in 1:6)for(j in i:6){s[i,j]<-s[j,i]<-t[k]-u[i]*u[j];k<-k+1L};s}
calc<-function(t,i,j,adjust){s<-covariance(t);if(adjust){c<-s[c(i,j),c(i,j),drop=FALSE]-s[c(i,j),5:6,drop=FALSE]%*%solve(s[5:6,5:6,drop=FALSE],s[5:6,c(i,j),drop=FALSE])}else c<-s[c(i,j),c(i,j),drop=FALSE];c[1,2]/sqrt(c[1,1]*c[2,2])}
grad<-function(f,t,step=1e-5)vapply(seq_along(t),function(k){h<-rep(0,length(t));h[k]<-step*max(abs(t[k]),.01);(f(t+h)-f(t-h))/(2*h[k])},numeric(1))
q<-a[a$eligible,];crit<-qt(.975,df=degf(sub));rows<-list()
for(pair in list(c(1,2),c(3,4))){i<-pair[1];j<-pair[2]
 for(adjust in c(FALSE,TRUE)){
  f<-function(t)calc(t,i,j,adjust);est<-f(v);g<-grad(f,v);stopifnot(max(abs(g-grad(f,v,1e-6)))<1e-4);se<-sqrt(drop(t(g)%*%V%*%g))
  x<-q[[vars[i]]];y<-q[[vars[j]]];w<-q$weight
  if(adjust){x<-residuals(lm(x~q$age+q$female,weights=w));y<-residuals(lm(y~q$age+q$female,weights=w))}
  mx<-weighted.mean(x,w);my<-weighted.mean(y,w);direct<-sum(w*(x-mx)*(y-my))/sqrt(sum(w*(x-mx)^2)*sum(w*(y-my)^2));stopifnot(abs(est-direct)<1e-10)
  rows[[length(rows)+1]]<-data.frame(fi=vars[i],fried=vars[j],adjusted=adjust,n=nrow(q),estimate=est,se=se,low=est-crit*se,high=est+crit*se,df=degf(sub))
 }
 f<-function(t)calc(t,i,j,TRUE)-calc(t,i,j,FALSE);g<-grad(f,v);se<-sqrt(drop(t(g)%*%V%*%g));est<-f(v)
 rows[[length(rows)+1]]<-data.frame(fi=vars[i],fried=vars[j],adjusted='adjusted-minus-unadjusted',n=nrow(q),estimate=est,se=se,low=est-crit*se,high=est+crit*se,df=degf(sub))
}
write.table(do.call(rbind,rows),file.path(out,'results.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('PASS: same sample; independent weighted residual correlations and two-step gradients agree; no individual output.\n');print(sessionInfo())
