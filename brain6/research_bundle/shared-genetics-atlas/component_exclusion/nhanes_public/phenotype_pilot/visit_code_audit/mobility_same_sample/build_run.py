#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent;base=p.parent/'score_corrected.R'
source=base.read_text();old="c('primary','mobility_proxy','visits_low','visits_high')";assert source.count(old)==2
source=source.replace(old,"c('mobility_zero','mobility_one')")
old="if(mode=='mobility_proxy')for(k in c('PFQB','PFQC')){sel<-is.na(mm[,k])&!is.na(d$PFQ054)&d$PFQ054==1;mm[sel,k]<-1}"
new="if(mode %in% c('mobility_zero','mobility_one'))for(k in c('PFQB','PFQC')){sel<-is.na(mm[,k])&!is.na(d$PFQ054)&d$PFQ054==1;mm[sel,k]<-if(mode=='mobility_zero') 0 else 1}"
assert source.count(old)==1;source=source.replace(old,new)
old="out<-'outputs/shared-genetics-atlas/component_exclusion/nhanes_public/phenotype_pilot/visit_code_audit/corrected'"
new="out<-'outputs/shared-genetics-atlas/component_exclusion/nhanes_public/phenotype_pilot/visit_code_audit/mobility_same_sample'"
assert source.count(old)==1;source=source.replace(old,new)
extra=r'''
# Evaluate paired correlation changes using joint survey moments on the same sample.
a0<-do.call(rbind,data_all[grep(' mobility_zero$',names(data_all))]);a1<-do.call(rbind,data_all[grep(' mobility_one$',names(data_all))])
stopifnot(!anyDuplicated(a0$SEQN),!anyDuplicated(a1$SEQN),setequal(a0$SEQN,a1$SEQN))
a1<-a1[match(a0$SEQN,a1$SEQN),];stopifnot(identical(a0$SEQN,a1$SEQN),identical(a0$STRATUM,a1$STRATUM),identical(a0$SDMVPSU,a1$SDMVPSU),identical(a0$weight,a1$weight))
elig<-function(a)!is.na(a$age)&a$age>=60&complete.cases(a[,c('fi','fi48','fi47','ffs','ffs4','sleep','fatigue')])
e0<-elig(a0);e1<-elig(a1);stopifnot(identical(e0,e1))
comparisons<-list(c('fi','ffs'),c('fi','sleep'),c('fi47','sleep'),c('fi47','ffs4'));result<-list()
for(pair in comparisons){
 data<-data.frame(STRATUM=a0$STRATUM,SDMVPSU=a0$SDMVPSU,weight=a0$weight,eligible=e0,x0=a0[[pair[1]]],y0=a0[[pair[2]]],x1=a1[[pair[1]]],y1=a1[[pair[2]]])
 d<-svydesign(id=~SDMVPSU,strata=~STRATUM,weights=~weight,nest=TRUE,data=data);d<-subset(d,eligible)
 m<-svymean(~x0+y0+x1+y1+I(x0^2)+I(y0^2)+I(x0*y0)+I(x1^2)+I(y1^2)+I(x1*y1),d)
 names(m)<-c('a','b','x','y','aa','bb','ab','xx','yy','xy');dimnames(attr(m,'var'))<-list(names(m),names(m))
 z<-svycontrast(m,quote((xy-x*y)/sqrt((xx-x*x)*(yy-y*y))-(ab-a*b)/sqrt((aa-a*a)*(bb-b*b))))
 est<-unname(coef(z));se<-unname(SE(z));crit<-qt(.975,df=degf(d))
 v<-as.numeric(m);fun<-function(t)(t[10]-t[3]*t[4])/sqrt((t[8]-t[3]^2)*(t[9]-t[4]^2))-(t[7]-t[1]*t[2])/sqrt((t[5]-t[1]^2)*(t[6]-t[2]^2))
 gradient<-vapply(seq_along(v),function(j){h<-rep(0,length(v));h[j]<-1e-6;(fun(v+h)-fun(v-h))/(2e-6)},numeric(1))
 stopifnot(abs(fun(v)-est)<1e-10,abs(sqrt(drop(t(gradient)%*%vcov(m)%*%gradient))-se)<1e-7)
 result[[length(result)+1]]<-data.frame(x=pair[1],y=pair[2],eligible_n=sum(e0),df=degf(d),difference_one_minus_zero=est,se=se,low=est-crit*se,high=est+crit*se)
}
write.table(do.call(rbind,result),file.path(out,'paired_mobility.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
# Post-hoc paired difference between the two FI-vs-FI47 sleep contrasts.
stopifnot(identical(a0$sleep,a1$sleep))
data<-data.frame(STRATUM=a0$STRATUM,SDMVPSU=a0$SDMVPSU,weight=a0$weight,eligible=e0,a=a0$fi,b=a0$fi47,x=a1$fi,y=a1$fi47,s=a0$sleep)
des<-subset(svydesign(id=~SDMVPSU,strata=~STRATUM,weights=~weight,nest=TRUE,data=data),eligible)
m<-svymean(~a+b+x+y+s+I(a^2)+I(b^2)+I(x^2)+I(y^2)+I(s^2)+I(a*s)+I(b*s)+I(x*s)+I(y*s),des)
names(m)<-c('a','b','x','y','s','aa','bb','xx','yy','ss','u','v','w','z');dimnames(attr(m,'var'))<-list(names(m),names(m))
change<-svycontrast(m,quote((w-x*s)/sqrt((xx-x*x)*(ss-s*s))-(z-y*s)/sqrt((yy-y*y)*(ss-s*s))-(u-a*s)/sqrt((aa-a*a)*(ss-s*s))+(v-b*s)/sqrt((bb-b*b)*(ss-s*s))))
est<-unname(coef(change));se<-unname(SE(change));crit<-qt(.975,df=degf(des))
v<-as.numeric(m);fun<-function(t)(t[13]-t[3]*t[5])/sqrt((t[8]-t[3]^2)*(t[10]-t[5]^2))-(t[14]-t[4]*t[5])/sqrt((t[9]-t[4]^2)*(t[10]-t[5]^2))-(t[11]-t[1]*t[5])/sqrt((t[6]-t[1]^2)*(t[10]-t[5]^2))+(t[12]-t[2]*t[5])/sqrt((t[7]-t[2]^2)*(t[10]-t[5]^2))
grad<-vapply(seq_along(v),function(j){h<-rep(0,length(v));h[j]<-1e-6;(fun(v+h)-fun(v-h))/(2e-6)},numeric(1))
stopifnot(abs(fun(v)-est)<1e-10,abs(sqrt(drop(t(grad)%*%vcov(m)%*%grad))-se)<1e-7)
write.table(data.frame(eligible_n=sum(e0),df=degf(des),difference_one_minus_zero=est,se=se,low=est-crit*se,high=est+crit*se),file.path(out,'paired_sleep_contrast_change.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
'''
marker="cat('Completed exploratory survey-weighted scoring; no participant-level outputs saved.\\n');print(sessionInfo())"
assert source.count(marker)==1;source=source.replace(marker,extra+'\n'+marker)
(p/'run_generated.R').write_text(source)
receipt={'corrected_source_sha256':hashlib.sha256(base.read_bytes()).hexdigest(),'generated_sha256':hashlib.sha256((p/'run_generated.R').read_bytes()).hexdigest(),'changes':'only mode list, mobility-imputation value, output directory and paired moment analysis; individual-level data remain in memory'}
(p/'build_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
