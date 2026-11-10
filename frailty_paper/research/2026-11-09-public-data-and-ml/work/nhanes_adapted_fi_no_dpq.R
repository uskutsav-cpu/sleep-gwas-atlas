.libPaths(c('work/Rlibs',.libPaths()));library(survey);library(foreign);options(survey.lonely.psu='fail')
yn<-c('CDQ010','CDQ001','BPQ020','MCQ160E','MCQ160F','MCQ160B','DIQ010','MCQ160A','MCQ160L','MCQ220','MCQ084')
dpq<-c('DPQ050','DPQ070','DPQ080','DPQ020','DPQ010','DPQ090');pfq<-paste0('PFQ061',c('F','G','K','J','A','B','R'));out<-list()
for (year in c(2011,2013)) {
 s<-if(year==2011)'G' else 'H';path<-paste0('work/nhanes',year);readx<-function(n)read.xport(file.path(path,paste0(n,'_',s,'.XPT')))
 a<-readx('DEMO');for (f in c('CDQ','BPQ','MCQ','DIQ','DPQ','PFQ','SLQ')) {b<-readx(f); vars<-intersect(c(yn,dpq,pfq,'SLQ050'),names(b)); stopifnot(!anyDuplicated(b$SEQN));a<-merge(a,b[,c('SEQN',vars),drop=FALSE],by='SEQN',all.x=TRUE)}
 stopifnot(all(c(yn,dpq,pfq,'SLQ050')%in%names(a)),!anyDuplicated(a$SEQN));a<-subset(a,!is.na(WTMEC2YR)&WTMEC2YR>0&!is.na(SDMVSTRA)&!is.na(SDMVPSU));a$female<-as.numeric(a$RIAGENDR==2);a$sleep_problem<-ifelse(a$SLQ050%in%c(1,2),as.numeric(a$SLQ050==1),NA_real_)
 for(v in yn)a[[paste0('i_',v)]]<-ifelse(a[[v]]%in%c(1,2),as.numeric(a[[v]]==1),NA_real_)
 for(v in dpq)a[[paste0('i_',v)]]<-ifelse(a[[v]]%in%0:3,as.numeric(a[[v]]>0),NA_real_)
 for(v in pfq)a[[paste0('i_',v)]]<-ifelse(a[[v]]%in%1:4,as.numeric(a[[v]]>1),NA_real_)
 primary_mat<-as.matrix(a[,paste0('i_',c(yn,dpq,pfq))]);a$primary_available<-rowSums(!is.na(primary_mat))>=20
 a$physical_available<-rowSums(!is.na(as.matrix(a[,paste0('i_',c(yn,pfq))])))>=15
 for(spec in c('primary','no_dpq')) {vars<-c(yn,dpq,pfq);if(spec=='no_dpq')vars<-setdiff(vars,dpq);mat<-as.matrix(a[,paste0('i_',vars)]); nobs<-rowSums(!is.na(mat));a$index<-rowMeans(mat,na.rm=TRUE);a$index[nobs<ceiling(.8*length(vars))]<-NA; a$domain<-a$RIDAGEYR>=60&a$RIAGENDR%in%c(1,2)&!is.na(a$sleep_problem)&!is.na(a$index)&a$primary_available&a$physical_available;a$domain[is.na(a$domain)]<-FALSE
 des<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC2YR,nest=TRUE,data=a);sub<-subset(des,domain);fit<-svyglm(index~sleep_problem+RIDAGEYR+female,design=sub);co<-summary(fit)$coefficients['sleep_problem',];ci<-confint(fit,parm='sleep_problem');out[[length(out)+1]]<-data.frame(cycle=paste0(year,'-',year+1),spec=spec,items=length(vars),eligible=sum(a$RIDAGEYR>=60,na.rm=TRUE),n=sum(a$domain),sleep_yes=sum(a$sleep_problem[a$domain]==1),design_df=degf(sub),estimate=unname(co['Estimate']),ci_low=ci[1,1],ci_high=ci[1,2],p=unname(co['Pr(>|t|)']))
 }
}
o<-do.call(rbind,out);write.csv(o,'outputs/stage6_nhanes_adapted_fi_no_dpq_results.csv',row.names=FALSE);print(o)
