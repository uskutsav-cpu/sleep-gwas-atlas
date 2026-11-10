.libPaths(c('work/Rlibs',.libPaths()));library(survey);library(foreign);options(survey.lonely.psu='fail')
yn<-c('CDQ010','CDQ001','BPQ020','MCQ160E','MCQ160F','MCQ160B','DIQ010','MCQ160A','MCQ160L','MCQ220','MCQ084')
dpq<-c('DPQ050','DPQ070','DPQ080','DPQ020','DPQ010','DPQ090');pfq<-paste0('PFQ061',c('F','G','K','J','A','B','R'));out<-list()
for (year in c(2011,2013)) {
 s<-if(year==2011)'G' else 'H';path<-paste0('work/nhanes',year);readx<-function(n)read.xport(file.path(path,paste0(n,'_',s,'.XPT')))
 a<-readx('DEMO');for (f in c('CDQ','BPQ','MCQ','DIQ','DPQ','PFQ','SLQ','PAQ','WHQ')) {b<-readx(f); vars<-intersect(c(yn,dpq,pfq,'SLQ050','PFQ061E','PFQ061H','DPQ040','PAD680','WHD020','WHD050','WHQ060'),names(b)); stopifnot(!anyDuplicated(b$SEQN));a<-merge(a,b[,c('SEQN',vars),drop=FALSE],by='SEQN',all.x=TRUE)}
 stopifnot(all(c(yn,dpq,pfq,'SLQ050')%in%names(a)),!anyDuplicated(a$SEQN));a<-subset(a,!is.na(WTMEC2YR)&WTMEC2YR>0&!is.na(SDMVSTRA)&!is.na(SDMVPSU));a$female<-as.numeric(a$RIAGENDR==2);a$sleep_problem<-ifelse(a$SLQ050%in%c(1,2),as.numeric(a$SLQ050==1),NA_real_)
 for(v in yn)a[[paste0('i_',v)]]<-ifelse(a[[v]]%in%c(1,2),as.numeric(a[[v]]==1),NA_real_)
 for(v in dpq)a[[paste0('i_',v)]]<-ifelse(a[[v]]%in%0:3,as.numeric(a[[v]]>0),NA_real_)
 for(v in pfq)a[[paste0('i_',v)]]<-ifelse(a[[v]]%in%1:4,as.numeric(a[[v]]>1),NA_real_)
 a$weak<-ifelse(a$PFQ061E%in%1:4,as.numeric(a$PFQ061E>1),NA_real_);a$slow<-ifelse(a$PFQ061H%in%1:4,as.numeric(a$PFQ061H>1),NA_real_);a$exhaust<-ifelse(a$DPQ040%in%0:3,as.numeric(a$DPQ040>=2),NA_real_)
 a$sedentary<-NA_real_;a$agegroup<-cut(a$RIDAGEYR,breaks=c(59,69,79,Inf),labels=c('60-69','70-79','80+'));for (ag in levels(a$agegroup)){idx<-!is.na(a$agegroup)&a$agegroup==ag&!is.na(a$PAD680)&a$PAD680>=0&a$PAD680<1440;cutoff<-quantile(a$PAD680[idx],.8,type=1,names=FALSE);a$sedentary[idx]<-as.numeric(a$PAD680[idx]>=cutoff);cat(year,ag,'sedentary cutoff',cutoff,'n',sum(idx),'flag',sum(a$sedentary[idx]),'\n')}
 validweight<-!is.na(a$WHD020)&!is.na(a$WHD050)&a$WHD020>=50&a$WHD020<=750&a$WHD050>=50&a$WHD050<=750;loss<-a$WHD050-a$WHD020;a$loss<-ifelse(validweight&loss<10,0,ifelse(validweight&loss>=10&a$WHQ060%in%c(1,2),as.numeric(a$WHQ060==2),NA_real_));a$count5<-a$weak+a$slow+a$exhaust+a$sedentary+a$loss;a$count4<-a$weak+a$slow+a$sedentary+a$loss
 for(spec in c('proportion_difference')) {vars<-c(yn,dpq,pfq);mat<-as.matrix(a[,paste0('i_',vars)]); nobs<-rowSums(!is.na(mat));a$index<-rowMeans(mat,na.rm=TRUE);a$index[nobs<ceiling(.8*length(vars))]<-NA; a$proportion_difference<-a$index-a$count5/5; a$domain<-a$RIDAGEYR>=60&a$RIAGENDR%in%c(1,2)&!is.na(a$sleep_problem)&!is.na(a$index)&!is.na(a$count5);a$domain[is.na(a$domain)]<-FALSE
 des<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC2YR,nest=TRUE,data=a);sub<-subset(des,domain);form<-proportion_difference~sleep_problem+RIDAGEYR+female;fit<-svyglm(form,design=sub);co<-summary(fit)$coefficients['sleep_problem',];ci<-confint(fit,parm='sleep_problem');out[[length(out)+1]]<-data.frame(cycle=paste0(year,'-',year+1),spec=spec,items=NA_integer_,eligible=sum(a$RIDAGEYR>=60,na.rm=TRUE),n=sum(a$domain),sleep_yes=sum(a$sleep_problem[a$domain]==1),design_df=degf(sub),estimate=unname(co['Estimate']),ci_low=ci[1,1],ci_high=ci[1,2],p=unname(co['Pr(>|t|)']))
 }
}
o<-do.call(rbind,out);write.csv(o,'outputs/stage6_nhanes_paired_proportion_contrast.csv',row.names=FALSE);print(o)
