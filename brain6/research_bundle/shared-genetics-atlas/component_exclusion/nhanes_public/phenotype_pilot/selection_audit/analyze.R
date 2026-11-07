#!/usr/bin/env Rscript
source('outputs/shared-genetics-atlas/component_exclusion/nhanes_public/phenotype_pilot/score_and_analyze.R')
out<-'outputs/shared-genetics-atlas/component_exclusion/nhanes_public/phenotype_pilot/selection_audit'
a<-do.call(rbind,data_all[grep(' primary$',names(data_all))]);stopifnot(!anyDuplicated(a$SEQN))
demo<-do.call(rbind,lapply(c('G','H'),function(s)readfile('DEMO',s)[,c('SEQN','RIAGENDR')]))
pfq<-do.call(rbind,lapply(c('G','H'),function(s)readfile('PFQ',s)[,c('SEQN','PFQ054')]))
slq<-do.call(rbind,lapply(c('G','H'),function(s)readfile('SLQ',s)[,c('SEQN','SLQ050')]))
stopifnot(!anyDuplicated(demo$SEQN),!anyDuplicated(pfq$SEQN),!anyDuplicated(slq$SEQN))
a$female<-code(demo$RIAGENDR[match(a$SEQN,demo$SEQN)],c('1'=0,'2'=1))
a$mobility<-code(pfq$PFQ054[match(a$SEQN,pfq$SEQN)],c('1'=1,'2'=0))
a$disclosure<-code(slq$SLQ050[match(a$SEQN,slq$SEQN)],c('1'=1,'2'=0))
a$sleep_any<-ifelse(is.na(a$sleep),NA,as.numeric(a$sleep>0))
a$older<-!is.na(a$age)&a$age>=60
a$eligible<-a$older & complete.cases(a[,c('fi','fi48','fi47','ffs','ffs4','sleep','fatigue')])
stopifnot(sum(a$eligible)==2433)
design<-svydesign(id=~SDMVPSU,strata=~STRATUM,weights=~weight,nest=TRUE,data=a)
older<-subset(design,older);m<-svymean(~eligible,older);v<-unname(coef(m)[2]);se<-unname(SE(m)[2]);crit<-qt(.975,df=degf(older))
direct<-sum(a$weight[a$older]*a$eligible[a$older])/sum(a$weight[a$older]);stopifnot(abs(v-direct)<1e-12)
write.table(data.frame(older_n=sum(a$older),included_n=sum(a$eligible),estimate=v,se=se,low=v-crit*se,high=v+crit*se,df=degf(older)),file.path(out,'inclusion.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
rows<-list()
for(metric in c('age','female','mobility','sleep_any','disclosure'))for(group in c(FALSE,TRUE)){
 mask<-a$older & (a$eligible==group) & !is.na(a[[metric]])
 sub<-design[mask,];z<-svymean(as.formula(paste0('~',metric)),sub);est<-unname(coef(z));se<-unname(SE(z));df<-degf(sub);crit<-qt(.975,df=df)
 direct<-sum(a$weight[mask]*a[[metric]][mask])/sum(a$weight[mask]);stopifnot(abs(est-direct)<1e-12)
 rows[[length(rows)+1]]<-data.frame(metric=metric,included=group,n=sum(mask),estimate=est,se=se,low=est-crit*se,high=est+crit*se,df=df)
}
write.table(do.call(rbind,rows),file.path(out,'groups.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('PASS: domain means match direct weighted calculations; no participant-level exports.\n')
