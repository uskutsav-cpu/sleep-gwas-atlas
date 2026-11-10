.libPaths(c('work/Rlibs',.libPaths()));library(survey);library(foreign);options(survey.lonely.psu='fail')
result<-list()
for (year in c(2011,2013)) {
 suffix<-if(year==2011)'G' else 'H';path<-paste0('work/nhanes',year);readx<-function(n) read.xport(file.path(path,paste0(n,'_',suffix,'.XPT')))
 d<-readx('DEMO');g<-readx('MGX');p<-readx('PAXDAY');stopifnot(!anyDuplicated(d$SEQN),!anyDuplicated(g$SEQN),!anyDuplicated(p[,c('SEQN','PAXDAYD')]))
 p<-subset(p,!is.na(PAXTMD)&PAXTMD>=1439&!is.na(PAXWWMD)&!is.na(PAXSWMD)&PAXWWMD+PAXSWMD>=1200)
 stopifnot(all(p$PAXSWMD>=0 & p$PAXSWMD<=1440))
 count<-aggregate(PAXSWMD~SEQN,p,function(x)c(days=length(x),mean=mean(x)/60));names(count)[2]<-'summary';count$days<-count$summary[,'days'];count$sleep_hours<-count$summary[,'mean'];count$summary<-NULL
 a<-merge(merge(d,g[,c('SEQN','MGDCGSZ')],by='SEQN',all.x=T),count,by='SEQN',all.x=T);a<-subset(a,!is.na(WTMEC2YR)&WTMEC2YR>0&!is.na(SDMVSTRA)&!is.na(SDMVPSU));a$female<-as.numeric(a$RIAGENDR==2);a$centered_sleep<-a$sleep_hours-8
 for (min_days in c(4,5)) {
  a$domain<-!is.na(a$RIDAGEYR)&a$RIDAGEYR>=60&a$RIAGENDR%in%c(1,2)&!is.na(a$MGDCGSZ)&!is.na(a$days)&a$days>=min_days
  design<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC2YR,nest=T,data=a); sub<-subset(design,domain);fit<-svyglm(MGDCGSZ~centered_sleep+RIDAGEYR+female,design=sub);co<-summary(fit)$coefficients['centered_sleep',];ci<-confint(fit,parm='centered_sleep')
  result[[length(result)+1]]<-data.frame(cycle=paste0(year,'-',year+1),min_days=min_days,older_examined=sum(a$RIDAGEYR>=60,na.rm=T),n=sum(a$domain),design_df=degf(sub),residual_df=fit$df.residual,mean_sleep_unweighted=mean(a$sleep_hours[a$domain]),effect_kg_per_hour=unname(co['Estimate']),se=unname(co['Std. Error']),ci_low=ci[1,1],ci_high=ci[1,2],p_value=unname(co['Pr(>|t|)']))
 }
}
o<-do.call(rbind,result);write.csv(o,'outputs/stage6_nhanes_device_sleep_grip_results.csv',row.names=F);print(o)
