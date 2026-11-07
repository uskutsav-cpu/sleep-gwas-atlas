#!/usr/bin/env Rscript
root<-'outputs/shared-genetics-atlas/component_exclusion/nhanes_public/phenotype_pilot/visit_code_audit'
e1<-new.env(parent=.GlobalEnv);e2<-new.env(parent=.GlobalEnv)
sys.source(file.path(root,'score_original.R'),envir=e1);sys.source(file.path(root,'score_corrected.R'),envir=e2)
a<-do.call(rbind,e1$data_all[grep(' primary$',names(e1$data_all))]);b<-do.call(rbind,e2$data_all[grep(' primary$',names(e2$data_all))]);b<-b[match(a$SEQN,b$SEQN),]
stopifnot(identical(a$SEQN,b$SEQN));eligible<-function(x)!is.na(x$age)&x$age>=60&complete.cases(x[,c('fi','fi48','fi47','ffs','ffs4','sleep','fatigue')])
ea<-eligible(a);eb<-eligible(b);stopifnot(identical(ea,eb),sum(ea)==2433)
h<-e1$readfile('HUQ','H');stopifnot(!anyDuplicated(h$SEQN));visit<-h$HUQ051[match(a$SEQN,h$SEQN)];code8<-!is.na(visit)&visit==8 & a$cycle=='H'
change<-ea & abs(a$fi-b$fi)>1e-12;stopifnot(all(!change|code8),all(!is.na(b$fi[ea])))
z<-a;z$eligible<-ea;z$fi_corrected<-b$fi
des<-survey::svydesign(id=~SDMVPSU,strata=~STRATUM,weights=~weight,nest=TRUE,data=z);sub<-subset(des,eligible)
m<-survey::svymean(~fi+fi_corrected+ffs+I(fi^2)+I(fi_corrected^2)+I(ffs^2)+I(fi*ffs)+I(fi_corrected*ffs),sub)
names(m)<-c('a','b','f','aa','bb','ff','af','bf');dimnames(attr(m,'var'))<-list(names(m),names(m))
diff<-survey::svycontrast(m,quote((bf-b*f)/sqrt((bb-b*b)*(ff-f*f))-(af-a*f)/sqrt((aa-a*a)*(ff-f*f))))
est<-unname(coef(diff));se<-unname(SE(diff));crit<-qt(.975,df=survey::degf(sub))
write.table(data.frame(older_code8=sum(code8 & !is.na(a$age)&a$age>=60),eligible_code8=sum(code8&ea),eligible_fi_changed=sum(change),same_eligible=sum(ea==eb),eligible=sum(ea),corr_difference=est,se=se,low=est-crit*se,high=est+crit*se),file.path(root,'comparison.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('PASS: eligibility identical and FI changed only in HUQ051 code 8 cases. No participant records exported.\n')
