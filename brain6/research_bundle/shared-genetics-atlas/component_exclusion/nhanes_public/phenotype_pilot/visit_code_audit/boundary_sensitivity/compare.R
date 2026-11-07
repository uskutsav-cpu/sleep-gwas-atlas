#!/usr/bin/env Rscript
root<-'outputs/shared-genetics-atlas/component_exclusion/nhanes_public/phenotype_pilot/visit_code_audit';out<-file.path(root,'boundary_sensitivity')
e1<-new.env(parent=.GlobalEnv);e2<-new.env(parent=.GlobalEnv)
sys.source(file.path(root,'score_corrected.R'),envir=e1);sys.source(file.path(out,'score_high.R'),envir=e2)
a<-do.call(rbind,e1$data_all[grep(' primary$',names(e1$data_all))]);b<-do.call(rbind,e2$data_all[grep(' primary$',names(e2$data_all))]);b<-b[match(a$SEQN,b$SEQN),];stopifnot(identical(a$SEQN,b$SEQN))
eligible<-function(x)!is.na(x$age)&x$age>=60&complete.cases(x[,c('fi','fi48','fi47','ffs','ffs4','sleep','fatigue')]);ea<-eligible(a);eb<-eligible(b);stopifnot(identical(ea,eb),sum(ea)==2433)
h<-e1$readfile('HUQ','H');stopifnot(!anyDuplicated(h$SEQN));visit<-h$HUQ051[match(a$SEQN,h$SEQN)];is3<-!is.na(visit)&visit==3&a$cycle=='H'
changed<-ea & abs(a$fi-b$fi)>1e-12;stopifnot(all(!changed|is3),all(b$fi[changed]>a$fi[changed]))
a$eligible<-ea;a$fi_high<-b$fi;a$fi47_high<-b$fi47
design<-survey::svydesign(id=~SDMVPSU,strata=~STRATUM,weights=~weight,nest=TRUE,data=a);sub<-subset(design,eligible);q<-a[ea,];crit<-qt(.975,df=survey::degf(sub))
jobs<-list(c('fi','fi_high','ffs'),c('fi','fi_high','sleep'),c('fi47','fi47_high','sleep'));rows<-list()
for(job in jobs){x<-job[1];y<-job[2];target<-job[3];low<-e1$correlate(sub,x,target);high<-e1$correlate(sub,y,target)
 for(nm in c(x,y)){xx<-q[[nm]];yy<-q[[target]];w<-q$weight;mx<-weighted.mean(xx,w);my<-weighted.mean(yy,w);v<-sum(w*(xx-mx)*(yy-my))/sqrt(sum(w*(xx-mx)^2)*sum(w*(yy-my)^2));stopifnot(abs(v-if(nm==x)low['estimate'] else high['estimate'])<1e-12)}
 f<-as.formula(paste0('~',x,'+',y,'+',target,'+I(',x,'^2)+I(',y,'^2)+I(',target,'^2)+I(',x,'*',target,')+I(',y,'*',target,')'))
 m<-survey::svymean(f,sub);names(m)<-c('a','b','s','aa','bb','ss','u','v');dimnames(attr(m,'var'))<-list(names(m),names(m))
 d<-survey::svycontrast(m,quote((v-b*s)/sqrt((bb-b*b)*(ss-s*s))-(u-a*s)/sqrt((aa-a*a)*(ss-s*s))));est<-unname(coef(d));se<-unname(SE(d));stopifnot(abs(est-(high['estimate']-low['estimate']))<1e-12)
 rows[[length(rows)+1]]<-data.frame(score=x,target=target,low_r=low['estimate'],high_r=high['estimate'],difference=est,se=se,low=est-crit*se,high=est+crit*se,row.names=NULL)
}
write.table(do.call(rbind,rows),file.path(out,'comparison.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(eligible=sum(ea),older_code3=sum(is3 & !is.na(a$age)&a$age>=60),eligible_code3=sum(is3&ea),changed=sum(changed)),file.path(out,'counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('PASS: same sample, changes limited to code 3, weighted points and paired differences independently checked; no individual outputs.\n')
