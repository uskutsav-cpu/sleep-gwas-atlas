source('outputs/shared-genetics-atlas/component_exclusion/nhanes_public/phenotype_pilot/visit_code_audit/sleep_covariance_decomposition/score_copy.R')
out<-'outputs/shared-genetics-atlas/component_exclusion/nhanes_public/phenotype_pilot/visit_code_audit/sleep_covariance_decomposition'
a<-do.call(rbind,data_all[grep(' primary$',names(data_all))]);stopifnot(!anyDuplicated(a$SEQN));a$eligible<-!is.na(a$age)&a$age>=60 & complete.cases(a[,c('fi','fi48','fi47','ffs','ffs4','sleep','fatigue')]);stopifnot(sum(a$eligible)==2433)
a$direct<-a$sleep/a$n_deficits;a$rest<-(1-1/a$n_deficits)*a$fi48
stopifnot(all(a$n_deficits[a$eligible]>=40),max(abs(a$fi[a$eligible]-a$direct[a$eligible]-a$rest[a$eligible]))<1e-12)
design<-svydesign(id=~SDMVPSU,strata=~STRATUM,weights=~weight,nest=TRUE,data=a);sub<-subset(design,eligible)
m<-svymean(~sleep+direct+rest+fi+I(sleep*direct)+I(sleep*rest)+I(sleep*fi),sub);names(m)<-c('s','d','r','f','sd','sr','sf');dimnames(attr(m,'var'))<-list(names(m),names(m))
exprs<-list(direct=quote(sd-s*d),rest=quote(sr-s*r),total=quote(sf-s*f));get<-function(e){x<-svycontrast(m,e);crit<-qt(.975,df=degf(sub));z<-unname(coef(x));se<-unname(SE(x));c(estimate=z,se=se,low=z-crit*se,high=z+crit*se)}
res<-do.call(rbind,lapply(exprs,get));stopifnot(abs(res['total','estimate']-sum(res[c('direct','rest'),'estimate']))<1e-12)
q<-a[a$eligible,];w<-q$weight;wcov<-function(x,y){mx<-sum(w*x)/sum(w);my<-sum(w*y)/sum(w);sum(w*(x-mx)*(y-my))/sum(w)}
stopifnot(abs(wcov(q$fi,q$sleep)-res['total','estimate'])<1e-12,abs(wcov(q$direct,q$sleep)-res['direct','estimate'])<1e-12,abs(wcov(q$rest,q$sleep)-res['rest','estimate'])<1e-12)
r<-correlate(sub,'fi','sleep');saved<-read.delim('outputs/shared-genetics-atlas/component_exclusion/nhanes_public/phenotype_pilot/visit_code_audit/corrected/correlations.tsv');target<-subset(saved,mode=='primary' & x=='fi' & y=='sleep')$estimate;stopifnot(length(target)==1,abs(r['estimate']-target)<1e-12)
write.table(data.frame(term=rownames(res),res,row.names=NULL),file.path(out,'covariances.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(n_people=nrow(q),n_min=min(q$n_deficits),n_median=median(q$n_deficits),n_max=max(q$n_deficits),df=degf(sub),fi_sleep_r=r['estimate']),file.path(out,'summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('Survey covariance terms and identity verified; no person-level data saved.\n');print(res);print(r)
