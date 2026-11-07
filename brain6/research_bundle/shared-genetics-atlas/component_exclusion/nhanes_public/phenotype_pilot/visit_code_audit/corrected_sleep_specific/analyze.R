#!/usr/bin/env Rscript
source('outputs/shared-genetics-atlas/component_exclusion/nhanes_public/phenotype_pilot/visit_code_audit/score_corrected.R')
out<-'outputs/shared-genetics-atlas/component_exclusion/nhanes_public/phenotype_pilot/visit_code_audit/corrected_sleep_specific'
a<-do.call(rbind,data_all[grep(' primary$',names(data_all))]);stopifnot(!anyDuplicated(a$SEQN))
slq<-rbind(readfile('SLQ','G')[,c('SEQN','SLQ050')],readfile('SLQ','H')[,c('SEQN','SLQ050')]);stopifnot(!anyDuplicated(slq$SEQN))
a$slq<-code(slq$SLQ050[match(a$SEQN,slq$SEQN)],c('1'=1,'2'=0))
a$base<-!is.na(a$age)&a$age>=60 & complete.cases(a[,c('fi','fi48','fi47','ffs','ffs4','sleep','fatigue')])
a$slq_domain<-a$base & !is.na(a$slq)
design<-svydesign(id=~SDMVPSU,strata=~STRATUM,weights=~weight,nest=TRUE,data=a)
paired_diff<-function(d,x,y,target){
 f<-as.formula(paste0('~',x,'+',y,'+',target,'+I(',x,'^2)+I(',y,'^2)+I(',target,'^2)+I(',x,'*',target,')+I(',y,'*',target,')'))
 m<-svymean(f,d);names(m)<-c('a','b','s','aa','bb','ss','u','v');dimnames(attr(m,'var'))<-list(names(m),names(m))
 z<-svycontrast(m,quote((u-a*s)/sqrt((aa-a*a)*(ss-s*s))-(v-b*s)/sqrt((bb-b*b)*(ss-s*s))))
 estimate<-unname(coef(z));se<-unname(SE(z));crit<-qt(.975,df=degf(d))
 x0<-as.numeric(m);fn<-function(t)(t[7]-t[1]*t[3])/sqrt((t[4]-t[1]^2)*(t[6]-t[3]^2))-(t[8]-t[2]*t[3])/sqrt((t[5]-t[2]^2)*(t[6]-t[3]^2))
 grad<-vapply(seq_along(x0),function(j){h<-rep(0,length(x0));h[j]<-1e-6;(fn(x0+h)-fn(x0-h))/(2e-6)},numeric(1))
 stopifnot(abs(sqrt(drop(t(grad)%*%vcov(m)%*%grad))-se)<1e-7)
 c(estimate=estimate,se=se,low=estimate-crit*se,high=estimate+crit*se)
}
allcorr<-list();alldiff<-list();counts<-list()
for(which in c('dpq','slq')){
 mask<-if(which=='dpq')a$base else a$slq_domain;sub<-design[mask,];target<-if(which=='dpq')'sleep' else 'slq'
 counts[[length(counts)+1]]<-data.frame(target=which,n=sum(mask),df=degf(sub))
 for(x in c('fi','fi48','fi47')){
  v<-correlate(sub,x,target);q<-a[mask,];wx<-q[[x]];wy<-q[[target]];wt<-q$weight
  mx<-sum(wt*wx)/sum(wt);my<-sum(wt*wy)/sum(wt)
  direct<-sum(wt*(wx-mx)*(wy-my))/sqrt(sum(wt*(wx-mx)^2)*sum(wt*(wy-my)^2))
  stopifnot(abs(v['estimate']-direct)<1e-12)
  allcorr[[length(allcorr)+1]]<-data.frame(target=which,score=x,t(v))
 }
 for(pair in list(c('fi','fi48'),c('fi48','fi47'),c('fi','fi47'))){
  v<-paired_diff(sub,pair[1],pair[2],target);alldiff[[length(alldiff)+1]]<-data.frame(target=which,score_a=pair[1],score_b=pair[2],t(v))
 }
}
write.table(do.call(rbind,counts),file.path(out,'counts.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,allcorr),file.path(out,'correlations.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,alldiff),file.path(out,'contrasts.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('Completed: direct weighted correlations agree; no participant rows exported.\n')
